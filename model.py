#!/usr/bin/env python3
"""Fit the CZ Aqr V-band light curve with PHOEBE. See docs/guide.md."""

import argparse
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize

# The project uses the Johnson V tables already installed in .venv-phoebe.
os.environ.setdefault("PHOEBE_ENABLE_ONLINE_PASSBANDS", "FALSE")
import phoebe

from plots import save_modeling_plots

PARAMETERS = ("incl_deg", "r1_fraction", "temperature_ratio", "phase_shift")


def load_observations(path, ephemeris):
    """Read accepted V-band measurements and calculate their HJD phases."""
    data = pd.read_csv(path)
    data = data[data["quality_flag"].eq("OK")].copy()
    if len(data) <= 5:
        raise ValueError("Need more than five accepted observations")
    if not data["filter"].eq("V").all():
        raise ValueError("This model expects V-band measurements")
    values = data[["hjd_utc", "v_mag", "v_mag_err"]].to_numpy()
    if not np.isfinite(values).all() or (data["v_mag_err"] <= 0).any():
        raise ValueError("Measurements must be finite, with positive errors")
    data["phase"] = ((data["hjd_utc"] - ephemeris["t0_hjd"])
                     / ephemeris["period_days"]) % 1.0
    return data.sort_values("hjd_utc").reset_index(drop=True)


def build_binary(config, period):
    """Circular, synchronous binary with a Roche-filling secondary."""
    adopted = config["adopted"]
    binary = phoebe.default_binary()
    settings = {
        "period@binary": period,
        "q@binary": adopted["mass_ratio"],
        "sma@binary": adopted["separation_Rsun"],
        "teff@primary": adopted["primary_temperature_K"],
        "ecc@binary": 0.0,
        "t0_supconj@binary": 0.0,
        "syncpar@primary": 1.0,
        "syncpar@secondary": 1.0,
        "gravb_bol@primary": 1.0,
        "gravb_bol@secondary": 0.32,
        "irrad_frac_refl_bol@primary": 1.0,
        "irrad_frac_refl_bol@secondary": 0.5,
    }
    for name, value in settings.items():
        binary.set_value(name, value)
    binary.add_constraint("semidetached", component="secondary")
    binary.add_dataset("lc", times=[0.0], passband="Johnson:V", dataset="lc01")
    binary.set_value_all("atm", "ck2004")
    binary.set_value_all("ld_mode", "interp")
    binary.set_value("irrad_method", "horvat")
    binary.set_value_all("fti_method", "none")  # Evaluate exposure midpoints.
    binary.set_value_all("ntriangles", config["compute"]["ntriangles"])
    return binary


def model_magnitudes(binary, parameters, phases, period):
    """Calculate model magnitudes; a positive phase shift moves eclipse later."""
    inclination, radius_fraction, temperature_ratio, phase_shift = parameters
    separation = binary.get_value("sma@binary@component")
    primary_temperature = binary.get_value("teff@primary@component")
    binary.set_value("incl@binary", inclination)
    binary.set_value("requiv@primary@component", radius_fraction * separation)
    binary.set_value("teff@secondary", temperature_ratio * primary_temperature)
    binary.set_value("t0_supconj@binary", phase_shift * period)

    # Compute sorted unique times, then restore the observations' original order.
    times, original_order = np.unique(np.asarray(phases) * period, return_inverse=True)
    binary.set_value("times@lc01@dataset", times)
    binary.run_compute(model="trial", overwrite=True, progressbar=False)
    flux = np.asarray(binary.get_value("fluxes@lc01@trial"))
    if not np.isfinite(flux).all() or (flux <= 0).any():
        raise ValueError("PHOEBE returned invalid model flux")
    return (-2.5 * np.log10(flux))[original_order]


def profile_offset(observed_mag, errors, model_mag):
    """Align the model's magnitude scale and return residuals divided by sigma."""
    difference = np.asarray(observed_mag) - model_mag
    offset = float(np.average(difference, weights=1 / np.square(errors)))
    return offset, np.asarray((difference - offset) / errors)


def fit_lightcurve(binary, data, config, period):
    """Fit four parameters with bounded Powell optimization."""
    bounds = np.array([config["bounds"][name] for name in PARAMETERS])
    start = np.array([config["initial"][name] for name in PARAMETERS])
    lower, upper = bounds.T
    if np.any(upper <= lower) or np.any((start < lower) | (start > upper)):
        raise ValueError("Check parameter bounds and initial values in model_config.json")
    best_chi2, best_parameters = np.inf, start.copy()

    def objective(scaled_parameters):
        nonlocal best_chi2, best_parameters
        parameters = lower + scaled_parameters * (upper - lower)
        predicted = model_magnitudes(binary, parameters, data["phase"], period)
        _, residuals = profile_offset(data["v_mag"], data["v_mag_err"], predicted)
        chi2 = float(np.sum(residuals**2))
        if chi2 < best_chi2:
            best_chi2, best_parameters = chi2, parameters.copy()
        return chi2

    # Scaling to [0, 1] makes degrees, radius, temperature ratio and phase comparable.
    options = config["optimizer"].copy()
    options["maxfev"] = options.pop("max_evaluations")
    result = minimize(objective, (start - lower) / (upper - lower), method="Powell",
                      bounds=[(0, 1)] * len(start), options=options)
    return best_parameters, {
        "optimizer_success": bool(result.success),
        "message": str(result.message),
        "evaluations": int(result.nfev),
    }


def derive_results(parameters, curve, config, ephemeris, observations):
    """Eclipse depths, a local eclipse epoch, and conditional stellar properties."""
    parameters = parameters.set_index("parameter")
    values = parameters["value"]
    status = ("model-derived" if parameters.loc[list(PARAMETERS), "status"].eq("fitted").all()
              else "trial-derived")
    brightest = float(curve["model_v_mag"].min())
    rows = [("Brightest model V magnitude", brightest, "mag", status,
             "Minimum magnitude of the saved sampled model; reference for both eclipse depths")]
    for name, center in [("Primary", values["phase_shift"]),
                         ("Secondary", values["phase_shift"] + 0.5)]:
        # Select the half-orbit around each eclipse, independent of display wrapping.
        distance = (curve["phase"] - center + 0.5) % 1.0 - 0.5
        faintest = float(curve.loc[distance.abs() < 0.25, "model_v_mag"].max())
        depth = faintest - brightest
        rows.extend([
            (f"{name} eclipse depth", depth, "mag", status,
             "Faintest sampled model magnitude in this half-orbit minus brightest model magnitude"),
            (f"{name} eclipse flux loss", 100 * (1 - 10**(-0.4 * depth)), "%", status,
             "Brightness reduction relative to the brightest sampled model: 100*(1-10**(-0.4*depth))"),
        ])
    primary_temperature = config["adopted"]["primary_temperature_K"]
    rows.append(("Secondary temperature", values["temperature_ratio"] * primary_temperature,
                 "K", status, f"Fitted T2/T1 times adopted T1={primary_temperature:g} K; not an independent temperature measurement"))
    separation = config["adopted"]["separation_Rsun"]
    rows.append(("Primary equivalent radius", values["r1_fraction"] * separation,
                 "solar radii", status, f"Fitted R1/a times adopted a={separation:g} solar radii; volume-equivalent radius, conditional on adopted separation"))

    period, t0 = ephemeris["period_days"], ephemeris["t0_hjd"]
    shift = values["phase_shift"]
    # Choose the fitted primary conjunction nearest the median observation time.
    cycle = int(np.rint((observations["hjd_utc"].median() - t0) / period - shift))
    rows.extend([
        ("Primary eclipse epoch", t0 + (cycle + shift) * period, "HJD (UTC)", status,
         "T0 + (E + phase_shift)*P; fitted primary conjunction nearest median observation HJD; timing uncertainty not estimated"),
        ("Primary eclipse epoch cycle", cycle, "cycles", status,
         "Integer E relative to the adopted T0; identifies the reported model-derived epoch"),
        ("Eclipse timing offset", shift * period * 1440, "min", status,
         "Fitted phase_shift*P*1440 relative to adopted ephemeris; negative means earlier; not a period-change measurement"),
    ])
    return pd.DataFrame(rows, columns=["quantity", "value", "unit", "status", "definition"])


def make_report_table(parameters, config, ephemeris, derived):
    """Round values for presentation while preserving their fitted/adopted status."""
    names = ["Inclination", "Primary radius / separation", "Secondary / primary temperature",
             "Phase shift", "Magnitude offset"]
    decimals = [2, 4, 4, 5, 4]
    parameters = parameters.set_index("parameter").loc[[*PARAMETERS, "magnitude_offset"]]
    rows = [(name, f"{row.value:.{digits}f}", row.unit, row.status)
            for name, digits, row in zip(names, decimals, parameters.itertuples())]
    adopted = config["adopted"]
    rows += [
        ("Mass ratio M2/M1", f"{adopted['mass_ratio']:.4f}", "1", "adopted"),
        ("Primary temperature", f"{adopted['primary_temperature_K']:.0f}", "K", "adopted"),
        ("Orbital separation", f"{adopted['separation_Rsun']:.2f}", "solar radii", "adopted"),
        ("Orbital period", f"{ephemeris['period_days']:.8f}", "days", "adopted"),
        ("Reference epoch", f"{ephemeris['t0_hjd']:.5f}", "HJD (UTC)", "adopted"),
        ("Eccentricity", "0", "1", "fixed"),
        ("Rotation / orbital frequency", "1", "1", "fixed"),
    ]
    for row in derived.itertuples():
        digits = {"mag": 3, "%": 1, "K": 0, "solar radii": 3,
                  "HJD (UTC)": 5, "cycles": 0, "min": 2}[row.unit]
        rows.append((row.quantity, f"{row.value:.{digits}f}", row.unit, row.status))
    return pd.DataFrame(rows, columns=["quantity", "value", "unit", "status"])


def run_modeling(table_path="output/tables/lightcurve_all_frames.csv", output_dir="output/modeling"):
    """Fit the observations with PHOEBE and save parameters, tables and plots."""
    photometry = json.loads(Path("config.json").read_text())
    config = json.loads(Path("model_config.json").read_text())
    period = photometry["ephemeris"]["period_days"]
    data = load_observations(table_path, photometry["ephemeris"])
    binary = build_binary(config, period)
    parameters, info = fit_lightcurve(binary, data, config, period)

    predicted = model_magnitudes(binary, parameters, data["phase"], period)
    offset, residuals = profile_offset(data["v_mag"], data["v_mag_err"], predicted)
    data["model_v_mag"] = predicted + offset
    data["residual_mag"] = data["v_mag"] - data["model_v_mag"]
    data["standardized_residual"] = residuals

    phase_start = photometry["output"]["phase_start"]
    phases = np.linspace(phase_start, phase_start + 1, config["compute"]["curve_points"])
    curve = pd.DataFrame({"phase": phases,
                          "model_v_mag": model_magnitudes(binary, parameters, phases, period) + offset})
    status = "fitted" if info["optimizer_success"] else "trial"
    parameter_table = pd.DataFrame({
        "parameter": [*PARAMETERS, "magnitude_offset"],
        "value": [*parameters, offset],
        "unit": ["deg", "1", "1", "cycles", "mag"],
        "status": [status] * 4 + ["profiled"],
    })
    chi2 = float(np.sum(residuals**2))
    degrees_of_freedom = len(data) - 5  # Four shape parameters and a magnitude offset.
    summary = {
        **info, "mode": "fit",
        "n_observations": len(data), "degrees_of_freedom": degrees_of_freedom,
        "chi2": chi2, "reduced_chi2": chi2 / degrees_of_freedom,
        "rms_residual_mag": float(np.sqrt(np.mean(data["residual_mag"]**2))),
    }
    output_dir = Path(output_dir)
    tables = output_dir / "tables"
    tables.mkdir(parents=True, exist_ok=True)
    data.to_csv(tables / "model_observations.csv", index=False)
    curve.to_csv(tables / "model_curve.csv", index=False)
    parameter_table.to_csv(tables / "parameters.csv", index=False)
    derived = derive_results(parameter_table, curve, config, photometry["ephemeris"], data)
    derived.to_csv(tables / "derived_results.csv", index=False)
    report = make_report_table(parameter_table, config, photometry["ephemeris"], derived)
    report.to_csv(tables / "report_parameters.csv", index=False)
    (output_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    save_modeling_plots(data, curve, output_dir, phase_start)
    print(f"{info['message']}; chi2/dof={summary['reduced_chi2']:.3f}; "
          f"RMS={summary['rms_residual_mag']:.5f} mag")
    print(f"Results saved to {output_dir.resolve()}")
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", default="output/tables/lightcurve_all_frames.csv")
    parser.add_argument("--output", default="output/modeling")
    args = parser.parse_args()
    summary = run_modeling(args.input, args.output)
    if not summary["optimizer_success"]:
        raise SystemExit(2)  # Outputs are saved, but marked as a trial rather than a fit.


if __name__ == "__main__":
    main()
