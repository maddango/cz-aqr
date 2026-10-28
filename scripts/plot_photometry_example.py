"""Create two illustrative report figures from one image and saved photometry."""

from pathlib import Path
import json
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Wedge
import matplotlib.patheffects as pe
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.stats import sigma_clipped_stats
from astropy.table import Table
from astropy.visualization import ImageNormalize, AsinhStretch
from astropy.wcs.utils import skycoord_to_pixel

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from photometry import open_science_frame, target_coord, estimate_fwhm_px, centroid_positions

EXAMPLE = "cpt0m439-sq39-20261004-0106-e91.fits.fz"


def main():
    config = json.loads((ROOT / "config.json").read_text())
    frames = pd.read_csv(ROOT / "output/tables/lightcurve_all_frames.csv")
    comparisons = pd.read_csv(ROOT / "output/tables/comparison_star_photometry.csv")
    row = frames.loc[frames["filename"].eq(EXAMPLE)].iloc[0]
    if row["quality_flag"] != "OK":
        raise ValueError("The example must be an accepted observation")
    stars = comparisons.loc[comparisons["filename"].eq(EXAMPLE)].copy()
    used = stars["used_for_zp"].to_numpy(dtype=bool)
    finite = np.isfinite(stars["inst_mag"]) & np.isfinite(stars["Vmag_catalog"])
    omitted = finite & ~used
    zp = float(row["zero_point"])
    individual_zp = stars["Vmag_catalog"] - stars["inst_mag"]
    scatter = float(individual_zp[used].std(ddof=1))
    output = ROOT / "output/plots"
    output.mkdir(parents=True, exist_ok=True)
    subtitle = f"{row['date_obs'][:10]}  {row['date_obs'][11:19]} UTC · {row['instrument']} · {row['exptime_s']:g} s"
    green, blue, orange = "#25834b", "#2468ad", "#e66732"

    fig, (left, right) = plt.subplots(1, 2, figsize=(10, 4.4), layout="constrained")
    x = stars.loc[finite, "inst_mag"]
    line_x = np.linspace(x.min() - 0.12, x.max() + 0.12, 100)
    left.plot(line_x, line_x + zp, color=blue, lw=1.6,
              label=rf"$V = m_{{\rm inst}} + {zp:.4f}$")
    left.scatter(stars.loc[used, "inst_mag"], stars.loc[used, "Vmag_catalog"],
                 color=green, s=30, zorder=3, label=f"Comparisons used ({used.sum()})")
    right.axhspan(zp - scatter, zp + scatter, color=blue, alpha=0.10,
                  label=rf"Comparison scatter: $\pm s$ ({scatter:.4f} mag)")
    right.axhline(zp, color=blue, lw=1.6, label=f"Weighted mean ZP = {zp:.4f}")
    right.scatter(stars.loc[used, "Vmag_catalog"], individual_zp[used],
                  color=green, s=30, zorder=3)
    if omitted.any():
        left.scatter(stars.loc[omitted, "inst_mag"], stars.loc[omitted, "Vmag_catalog"],
                     color="0.55", marker="x", s=30, label="Not used")
        right.scatter(stars.loc[omitted, "Vmag_catalog"], individual_zp[omitted],
                      color="0.55", marker="x", s=30, label="Not used")
    left.set(xlabel="Instrumental magnitude", ylabel="APASS V magnitude",
             title="Adopted calibration · slope fixed at 1")
    right.set(xlabel="APASS V magnitude", ylabel="APASS V − instrumental magnitude",
              title="Individual zero points")
    right.margins(y=0.25)
    for ax in (left, right):
        ax.grid(alpha=0.18)
        ax.legend(fontsize=8, loc="best")
        ax.ticklabel_format(useOffset=False, axis="both")
    fig.suptitle("CZ Aqr — calibration example\n" + subtitle, fontsize=12)
    save(fig, output / "zero_point_example")

    frame = open_science_frame(ROOT / "inputs" / EXAMPLE)
    catalog = Table.read(ROOT / config["comparison"]["catalog_path"], format="ascii.ecsv")
    catalog_ids = [f"{ra:.7f}_{dec:.7f}" for ra, dec in zip(catalog["RAJ2000"], catalog["DEJ2000"])]
    used_ids = set(stars.loc[used, "star_id"])
    ref_ids = config["uncertainty"]["reference_star_ids"]
    coordinates = SkyCoord(catalog["RAJ2000"], catalog["DEJ2000"], unit="deg")
    x, y = skycoord_to_pixel(coordinates, frame.wcs, origin=0)
    box = int(config["photometry"]["centroid_box_size"])
    box += int(box % 2 == 0)
    x, y = centroid_positions(frame.data, x, y, box, frame.bad_pixels)
    tx, ty = skycoord_to_pixel(target_coord(config), frame.wcs, origin=0)
    tx, ty = centroid_positions(frame.data, np.array([tx]), np.array([ty]), box, frame.bad_pixels)
    tx, ty = float(tx[0]), float(ty[0])
    fwhm, _ = estimate_fwhm_px(frame.header, config)
    radii = [config["photometry"][key] * fwhm for key in
             ("aperture_radius_fwhm", "annulus_inner_fwhm", "annulus_outer_fwhm")]
    _, median, std = sigma_clipped_stats(frame.data, sigma=3.0, mask=frame.bad_pixels)
    norm = ImageNormalize(vmin=median - std, vmax=median + 35 * std, stretch=AsinhStretch())
    fig, ax = plt.subplots(figsize=(8, 8), layout="constrained")
    ax.imshow(frame.data, origin="lower", cmap="gray", norm=norm, rasterized=True)
    normal = np.array([sid in used_ids and sid not in ref_ids for sid in catalog_ids])
    ax.scatter(x[normal], y[normal], s=170, facecolors="none", edgecolors="#6bdd8a",
               linewidths=1.2, label="APASS comparisons")
    special = [(tx, ty, "CZ Aqr", orange)]
    for index, sid in enumerate(ref_ids, 1):
        j = catalog_ids.index(sid)
        special.append((float(x[j]), float(y[j]), f"REF{index}", "#6cd9ff"))
    for sx, sy, label, color in special:
        ax.scatter([sx], [sy], s=220, facecolors="none", edgecolors=color, linewidths=1.6)
        ax.annotate(label, (sx, sy), xytext=(11, 10), textcoords="offset points",
                    fontsize=10, color=color, weight="bold",
                    path_effects=[pe.withStroke(linewidth=2.5, foreground="black")])
    ax.set(xlabel="x [pixel]", ylabel="y [pixel]", title="CZ Aqr — field and aperture photometry\n" + subtitle)
    ax.legend(loc="upper left", fontsize=9, facecolor="white", framealpha=0.9)

    # Only the inset circles represent actual aperture/annulus radii in pixels.
    inset = ax.inset_axes([0.69, 0.69, 0.28, 0.28])
    half = int(np.ceil(radii[2] * 1.35))
    ix, iy = int(round(tx)), int(round(ty))
    patch = frame.data[iy-half:iy+half+1, ix-half:ix+half+1]
    extent = (ix-half-.5-tx, ix+half+.5-tx, iy-half-.5-ty, iy+half+.5-ty)
    patch_norm = ImageNormalize(vmin=median - std, vmax=patch.max(), stretch=AsinhStretch())
    inset.imshow(patch, origin="lower", cmap="gray", norm=patch_norm,
                  extent=extent, interpolation="nearest", rasterized=True)
    inset.add_patch(Wedge((0, 0), radii[2], 0, 360, width=radii[2]-radii[1],
                          facecolor="#6cd9ff", alpha=0.16, edgecolor="none"))
    inset.add_patch(Circle((0, 0), radii[0], fill=False, edgecolor=orange, lw=1.5))
    for radius in radii[1:]:
        inset.add_patch(Circle((0, 0), radius, fill=False, edgecolor="#6cd9ff", lw=1.2, ls="--"))
    inset.set(xticks=[], yticks=[], title="CZ Aqr · enlarged")
    inset.title.set(color="white", fontsize=10,
                    path_effects=[pe.withStroke(linewidth=2, foreground="black")])
    inset.text(0.5, -0.05, "Orange: stellar aperture\nBlue: sky annulus", transform=inset.transAxes,
               ha="center", va="top", fontsize=9, color="white",
               path_effects=[pe.withStroke(linewidth=2, foreground="black")])
    save(fig, output / "field_example")
    print(f"Exposure: {EXAMPLE}")
    print(f"ZP = {zp:.8f}; comparison scatter = {scatter:.8f} mag")
    print(f"Aperture / annulus radii: {radii} pixels")
    print(f"Figures saved as PNG and PDF in {output}")


def save(fig, stem):
    fig.savefig(stem.with_suffix(".png"), dpi=180)
    fig.savefig(stem.with_suffix(".pdf"))
    plt.close(fig)


if __name__ == "__main__":
    main()
