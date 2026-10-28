#!/usr/bin/env python3
"""CZ Aqr aperture photometry from reduced LCO SCI images and a fixed APASS table.

Save per-frame and comparison measurements, basic uncertainties, mid-exposure
HJD/BJD times, light curves, and optional source-position diagnostics.
"""

import argparse
import json
import math
import warnings
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from astropy import units as u
from astropy.coordinates import SkyCoord, EarthLocation
from astropy.io import fits
from astropy.stats import sigma_clip
from astropy.table import Table
from astropy.time import Time
from astropy.wcs import WCS
from astropy.wcs.utils import skycoord_to_pixel

from photutils.aperture import CircularAperture, CircularAnnulus
from photutils.centroids import centroid_sources, centroid_com

from plots import save_diagnostic, make_summary_plots, save_zero_point_plot, save_reference_plot


FITS_PATTERNS = ("*.fits", "*.fit", "*.fits.fz", "*.fit.fz")


@dataclass
class FrameData:
    path: Path
    data: np.ndarray
    error: np.ndarray
    bad_pixels: np.ndarray
    header: fits.Header
    wcs: WCS


def find_fits(input_dir: str | Path) -> list[Path]:
    input_dir = Path(input_dir)
    files = []
    for pattern in FITS_PATTERNS:
        files.extend(input_dir.rglob(pattern))
    # Resolve paths so symlinked images are measured only once.
    return sorted(set(p.resolve() for p in files))


def open_science_frame(path: str | Path) -> FrameData:
    path = Path(path)
    with fits.open(path, memmap=False) as hdul:
        data = np.array(hdul["SCI"].data, dtype=float, copy=True)
        error = np.array(hdul["ERR"].data, dtype=float, copy=True)
        bad_pixels = np.asarray(hdul["BPM"].data) != 0
        header = hdul["SCI"].header.copy()

    if data.ndim != 2 or not np.isfinite(data).any():
        raise ValueError(f"SCI must contain a 2-D image with finite pixels: {path.name}")
    if error.shape != data.shape or bad_pixels.shape != data.shape:
        raise ValueError(f"SCI, ERR and BPM must have matching shapes: {path.name}")
    bad_pixels |= ~np.isfinite(data) | ~np.isfinite(error) | (error <= 0)

    return FrameData(path=path, data=data, error=error, bad_pixels=bad_pixels,
                     header=header, wcs=WCS(header))


def target_coord(config):
    return SkyCoord(config["target"]["ra"], config["target"]["dec"],
                    unit=(u.hourangle, u.deg))


def frame_times(header: fits.Header, coord: SkyCoord) -> dict:
    exptime = float(header["EXPTIME"])
    location = EarthLocation.from_geocentric(
        *(float(header[key]) for key in ("OBSGEO-X", "OBSGEO-Y", "OBSGEO-Z")), unit=u.m
    )
    t0 = Time(header["DATE-OBS"], format="isot", scale="utc", location=location)
    tmid = t0 + 0.5 * exptime * u.s

    ltt_bary = tmid.light_travel_time(coord, kind="barycentric")
    bjd_tdb = (tmid.tdb + ltt_bary).jd

    # Historical eclipse ephemerides are often quoted as HJD in UTC.
    ltt_helio = tmid.light_travel_time(coord, kind="heliocentric")
    hjd_utc = (tmid.utc + ltt_helio).jd

    return {
        "date_obs": header["DATE-OBS"],
        "hjd_utc": float(hjd_utc),
        "bjd_tdb": float(bjd_tdb),
    }


def estimate_fwhm_px(header, config):
    """Use the image seeing estimate, or the configured fallback aperture scale."""
    pixscale = float(header.get("PIXSCALE", np.nan))
    fwhm_arcsec = float(header.get("L1FWHM", np.nan))
    lower, upper = config["photometry"]["valid_fwhm_arcsec"]
    if np.isfinite(pixscale) and pixscale > 0 and lower <= fwhm_arcsec <= upper:
        return fwhm_arcsec / pixscale, fwhm_arcsec
    fwhm_px = config["photometry"]["fallback_fwhm_px"]
    return fwhm_px, fwhm_px * pixscale


def centroid_positions(data, x, y, box_size, bad_pixels):
    """Refine WCS positions, keeping the originals when centroiding fails."""
    try:
        xc, yc = centroid_sources(data, x, y, box_size=box_size,
                                  mask=bad_pixels, centroid_func=centroid_com)
        valid = np.isfinite(xc) & np.isfinite(yc)
        return np.where(valid, xc, x), np.where(valid, yc, y)
    except ValueError:
        return x, y


def photometer_positions(data, wcs, coords, fwhm_px, exptime, config, error, bad_pixels):
    """Measure sky-subtracted flux, SNR and instrumental magnitude for each star."""
    settings = config["photometry"]
    aperture_radius = float(settings["aperture_radius_fwhm"]) * fwhm_px
    annulus_inner = float(settings["annulus_inner_fwhm"]) * fwhm_px
    annulus_outer = float(settings["annulus_outer_fwhm"]) * fwhm_px

    x, y = skycoord_to_pixel(coords, wcs, origin=0)
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    margin = annulus_outer + 2
    in_frame = (
        np.isfinite(x) & np.isfinite(y)
        & (x > margin) & (x < data.shape[1] - margin)
        & (y > margin) & (y < data.shape[0] - margin)
    )

    # Recenter only sources that are safely in frame.
    x_center, y_center = x.copy(), y.copy()
    if np.any(in_frame):
        box = int(settings["centroid_box_size"])
        if box % 2 == 0:
            box += 1
        xx, yy = centroid_positions(data, x[in_frame], y[in_frame], box, bad_pixels)
        x_center[in_frame] = xx
        y_center[in_frame] = yy

    positions = np.column_stack([x_center, y_center])
    apertures = CircularAperture(positions, r=aperture_radius)
    annuli = CircularAnnulus(positions, r_in=annulus_inner, r_out=annulus_outer)

    aperture_masks = apertures.to_mask(method="exact")
    annulus_masks = annuli.to_mask(method="center")
    aperture_sum = np.full(len(coords), np.nan)
    pixel_variance = np.full(len(coords), np.nan)
    aperture_bad = np.zeros(len(coords), dtype=bool)
    sky_mean = np.full(len(coords), np.nan)
    sky_mean_variance = np.full(len(coords), np.nan)
    n_sky = np.zeros(len(coords), dtype=int)

    for i, (ap_mask, ann_mask) in enumerate(zip(aperture_masks, annulus_masks)):
        if not in_frame[i]:
            continue
        inside = ap_mask.data > 0
        aperture_bad[i] = np.any(ap_mask.cutout(bad_pixels, fill_value=True)[inside])
        # Reject compromised apertures rather than silently losing stellar flux.
        if aperture_bad[i]:
            continue
        weights = ap_mask.data[inside]
        aperture_sum[i] = np.sum(weights * ap_mask.cutout(data)[inside])
        # Propagate the actual weighted sum: squared fractional-pixel weights.
        pixel_variance[i] = np.sum((weights * ap_mask.cutout(error)[inside])**2)

        sky_good = (ann_mask.data > 0) & ~ann_mask.cutout(bad_pixels, fill_value=True)
        sky_values = ann_mask.cutout(data)[sky_good]
        errors = ann_mask.cutout(error)[sky_good]
        if sky_values.size < 10:
            continue
        clipped = sigma_clip(sky_values, sigma=3.0, maxiters=5)
        keep = ~np.ma.getmaskarray(clipped)
        n_sky[i] = int(keep.sum())
        if n_sky[i] < 10:
            continue
        sky_mean[i] = np.mean(sky_values[keep])
        sky_mean_variance[i] = np.sum(errors[keep]**2) / n_sky[i]**2

    aperture_area = float(apertures.area)
    flux = aperture_sum - sky_mean * aperture_area
    # ERR already contains pixel shot/read/calibration noise. Add only the
    # uncertainty of the local sky level estimated from the annulus.
    sky_estimate_variance = aperture_area**2 * sky_mean_variance
    flux_err = np.sqrt(pixel_variance + sky_estimate_variance)

    snr = flux / flux_err
    inst_mag = np.where(flux > 0, -2.5 * np.log10(flux / exptime), np.nan)
    mag_err = np.where(snr > 0, (2.5 / np.log(10)) / snr, np.nan)

    # Peak in a small region around each recentered coordinate, useful for QC.
    peak = np.full(len(coords), np.nan)
    half = max(2, int(math.ceil(fwhm_px)))
    for i, (xx, yy) in enumerate(zip(x_center, y_center)):
        if not in_frame[i]:
            continue
        ix, iy = int(round(xx)), int(round(yy))
        region = data[iy-half:iy+half+1, ix-half:ix+half+1]
        peak[i] = np.nanmax(region)

    return {
        "x": x_center, "y": y_center, "in_frame": in_frame,
        "flux": flux, "flux_err": flux_err,
        "snr": snr, "inst_mag": inst_mag, "mag_err": mag_err,
        "n_sky": n_sky, "aperture_bad": aperture_bad, "peak": peak,
    }

def robust_zero_point(comparisons, catalog, config, saturation):
    """Weighted APASS zero point after quality cuts and sigma clipping."""
    catalog_mag = np.ma.filled(catalog["Vmag"], np.nan).astype(float)
    catalog_error = np.ma.filled(catalog["e_Vmag"], np.nan).astype(float)
    magnitude = comparisons["inst_mag"]
    error = comparisons["mag_err"]
    quality = config["quality"]
    good = (
        comparisons["in_frame"] & ~comparisons["aperture_bad"]
        & np.isfinite(magnitude) & np.isfinite(error) & (error > 0)
        & np.isfinite(catalog_mag) & np.isfinite(catalog_error) & (catalog_error >= 0)
        & (comparisons["snr"] >= quality["comparison_min_snr"])
        & (comparisons["peak"] < quality["saturation_fraction"] * saturation)
    )
    zero_points = catalog_mag - magnitude
    used = np.zeros(len(catalog), dtype=bool)
    if good.sum() < quality["min_comparison_stars"]:
        return np.nan, used
    clipped = sigma_clip(zero_points[good], sigma=config["comparison"]["zero_point_clip_sigma"],
                         maxiters=5)
    used[np.flatnonzero(good)[~np.ma.getmaskarray(clipped)]] = True
    if used.sum() < quality["min_comparison_stars"]:
        return np.nan, used
    weights = 1 / (error[used]**2 + catalog_error[used]**2)
    weights /= weights.sum()
    return float(np.sum(weights * zero_points[used])), used


def process_one(path, catalog, config, output_dir):
    """Measure CZ Aqr and the comparison stars in one image."""
    frame = open_science_frame(path)
    header = frame.header
    target = target_coord(config)
    exposure = float(header["EXPTIME"])
    times = frame_times(header, target)
    fwhm_px, fwhm_arcsec = estimate_fwhm_px(header, config)
    saturation = float(header.get("MAXLIN", header.get("SATURATE", np.inf)))
    comparison_coords = SkyCoord(catalog["RAJ2000"], catalog["DEJ2000"], unit=u.deg)

    target_phot = photometer_positions(
        frame.data, frame.wcs, SkyCoord([target]), fwhm_px, exposure, config,
        frame.error, frame.bad_pixels,
    )
    comparison_phot = photometer_positions(
        frame.data, frame.wcs, comparison_coords, fwhm_px, exposure, config,
        frame.error, frame.bad_pixels,
    )
    zero_point, used = robust_zero_point(comparison_phot, catalog, config, saturation)
    # The target measurement contains one value per field.
    target_phot = {name: values[0] for name, values in target_phot.items()}
    magnitude = target_phot["inst_mag"] + zero_point

    quality = config["quality"]
    flags = []
    if not target_phot["in_frame"]:
        flags.append("TARGET_OUT_OF_FRAME")
    if target_phot["aperture_bad"]:
        flags.append("BAD_TARGET_APERTURE")
    if not np.isfinite(target_phot["mag_err"]) or target_phot["mag_err"] <= 0:
        flags.append("INVALID_MAG_ERROR")
    if not np.isfinite(target_phot["flux"]) or target_phot["flux"] <= 0:
        flags.append("NONPOSITIVE_TARGET_FLUX")
    if target_phot["snr"] < quality["target_min_snr"]:
        flags.append("LOW_TARGET_SNR")
    if target_phot["peak"] >= quality["saturation_fraction"] * saturation:
        flags.append("TARGET_NEAR_SATURATION")
    if used.sum() < quality["min_comparison_stars"]:
        flags.append("TOO_FEW_COMPARISONS")
    if not np.isfinite(zero_point):
        flags.append("NO_ZERO_POINT")

    ephemeris = config["ephemeris"]
    phase = ((times["hjd_utc"] - ephemeris["t0_hjd"]) / ephemeris["period_days"]) % 1.0
    row = {
        "filename": path.name, **times, "phase": phase,
        "exptime_s": exposure,
        "filter": header.get("FILTER", ""),
        "instrument": header.get("INSTRUME", ""),
        "airmass": header.get("AIRMASS", np.nan),
        "fwhm_arcsec": fwhm_arcsec,
        "maxlin_counts": saturation,
        "target_flux_counts": target_phot["flux"],
        "target_flux_err_counts": target_phot["flux_err"],
        "target_snr": target_phot["snr"],
        "target_inst_mag": target_phot["inst_mag"],
        "target_inst_mag_err": target_phot["mag_err"],
        "zero_point": zero_point,
        "n_comp": int(used.sum()),
        "v_mag": magnitude,
        "quality_flag": ";".join(flags) if flags else "OK",
    }
    comparison_rows = []
    for i, coord in enumerate(comparison_coords):
        comparison_rows.append({
            "filename": path.name,
            "star_id": f"{coord.ra.deg:.7f}_{coord.dec.deg:.7f}",
            "Vmag_catalog": float(catalog["Vmag"][i]),
            "inst_mag": comparison_phot["inst_mag"][i],
            "inst_mag_err": comparison_phot["mag_err"][i],
            "snr": comparison_phot["snr"][i],
            "peak_counts": comparison_phot["peak"][i],
            "bad_aperture": comparison_phot["aperture_bad"][i],
            "used_for_zp": used[i],
        })

    if config["output"]["save_diagnostics"]:
        diagnostics = output_dir / "diagnostics"
        diagnostics.mkdir(parents=True, exist_ok=True)
        comparison_xy = np.column_stack([comparison_phot["x"], comparison_phot["y"]])
        save_diagnostic(frame, (target_phot["x"], target_phot["y"]), comparison_xy,
                        used, row, diagnostics / f"{path.name}.png")
    return row, comparison_rows


def apply_reference_errors(frames: pd.DataFrame, comparisons: pd.DataFrame, config: dict):
    """Class prescription: target SNR error plus the temporal two-star scatter.

    Use one fixed pair, matched by filename, and sample standard deviation.
    Zero-point clipping must not suppress the reference differences.
    """
    ref1, ref2 = config["uncertainty"]["reference_star_ids"]
    if ref1 == ref2:
        raise ValueError("Choose two distinct uncertainty reference stars")
    accepted = frames.loc[frames["quality_flag"].eq("OK"),
                          ["filename", "bjd_tdb", "instrument", "maxlin_counts"]]
    stars = comparisons.merge(accepted, on="filename", validate="many_to_one")
    valid = (
        ~stars["bad_aperture"] & np.isfinite(stars["inst_mag"])
        & (stars["snr"] >= config["quality"]["comparison_min_snr"])
        & (stars["peak_counts"] < config["quality"]["saturation_fraction"] * stars["maxlin_counts"])
    )
    stars = stars.loc[valid]
    stars = stars[stars["star_id"].isin([ref1, ref2])]
    pair = stars.pivot(index="filename", columns="star_id", values="inst_mag")
    pair = pair.reindex(columns=[ref1, ref2]).dropna()
    pair = pair.rename(columns={ref1: "ref1_inst_mag", ref2: "ref2_inst_mag"}).reset_index()
    pair = pair.merge(accepted[["filename", "bjd_tdb", "instrument"]], on="filename")
    pair = pair.sort_values("bjd_tdb").reset_index(drop=True)
    if len(pair) < 3:
        raise ValueError("Need at least three accepted frames with both reference stars")
    pair["ref1_star_id"], pair["ref2_star_id"] = ref1, ref2
    pair["difference_mag"] = pair["ref1_inst_mag"] - pair["ref2_inst_mag"]
    instrumental_err = float(pair["difference_mag"].std(ddof=1))
    result = frames.copy()
    result["instrumental_err"] = instrumental_err
    result["v_mag_err"] = np.hypot(result["target_inst_mag_err"], instrumental_err)
    return result, pair


def run_pipeline(input_dir: str | Path, output_dir: str | Path):
    output_dir = Path(output_dir)
    (output_dir / "tables").mkdir(parents=True, exist_ok=True)

    config = json.loads(Path("config.json").read_text())
    files = find_fits(input_dir)
    if not files:
        raise FileNotFoundError(f"No FITS/FITS.FZ files found under {input_dir}")

    catalog = Table.read(config["comparison"]["catalog_path"], format="ascii.ecsv")
    rows, comparison_rows = [], []

    for i, path in enumerate(files, start=1):
        print(f"[{i}/{len(files)}] {path.name}")
        try:
            row, star_rows = process_one(path, catalog, config, output_dir)
            rows.append(row)
            comparison_rows.extend(star_rows)
            print(f"    V={row['v_mag']:.4f}, Ncomp={row['n_comp']}, {row['quality_flag']}")
        except Exception as exc:
            warnings.warn(f"{path.name}: {exc}")
            rows.append({
                "filename": path.name,
                "quality_flag": f"PROCESSING_ERROR:{type(exc).__name__}:{exc}",
            })

    frames = pd.DataFrame(rows)
    if not comparison_rows:
        raise ValueError("No images could be measured; inspect the warnings above")
    frames = frames.sort_values("bjd_tdb", na_position="last")
    comparisons = pd.DataFrame(comparison_rows)
    frames, pair = apply_reference_errors(frames, comparisons, config)
    frames.to_csv(output_dir / "tables" / "lightcurve_all_frames.csv", index=False)
    comparisons.to_csv(
        output_dir / "tables" / "comparison_star_photometry.csv", index=False
    )
    pair.to_csv(output_dir / "tables" / "reference_pair.csv", index=False)
    print(f"\nReference-pair scatter: {frames['instrumental_err'].iloc[0]:.5f} mag "
          f"from {len(pair)} paired observations")

    phase_start = float(config["output"]["phase_start"])
    make_summary_plots(frames, output_dir, phase_start)
    save_zero_point_plot(frames, comparisons, output_dir)
    save_reference_plot(pair, output_dir)

    print(f"\nDone. Results written to: {output_dir.resolve()}")
    return frames


def main():
    parser = argparse.ArgumentParser(description="Batch aperture photometry for CZ Aqr.")
    parser.add_argument("input_dir", help="Directory containing reduced FITS/FITS.FZ images")
    parser.add_argument("--output", default="output", help="Output directory")
    args = parser.parse_args()

    run_pipeline(args.input_dir, args.output)


if __name__ == "__main__":
    main()
