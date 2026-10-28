"""Light-curve plots and image diagnostics for the photometry pipeline."""

from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from astropy.stats import sigma_clipped_stats


def save_figure(fig, path, dpi=180, pdf=True):
    """Save the PNG and a report-quality PDF, then close the figure."""
    fig.tight_layout()
    fig.savefig(path, dpi=dpi)
    if pdf:
        fig.savefig(Path(path).with_suffix(".pdf"))
    plt.close(fig)


def save_diagnostic(frame, target_xy, comparison_xy, used_mask, row, output_path):
    data = frame.data
    _, median, std = sigma_clipped_stats(data, sigma=3.0)
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.imshow(
        data, origin="lower", cmap="gray",
        vmin=median - 2 * std, vmax=median + 8 * std,
    )
    ax.scatter([target_xy[0]], [target_xy[1]],
               s=180, facecolors="none", edgecolors="C3", linewidths=2,
               label="CZ Aqr")
    if len(comparison_xy):
        ax.scatter(comparison_xy[:, 0], comparison_xy[:, 1],
                   s=90, facecolors="none", edgecolors="C0",
                   label="APASS candidates")
        if np.any(used_mask):
            ax.scatter(comparison_xy[used_mask, 0], comparison_xy[used_mask, 1],
                       s=150, facecolors="none", edgecolors="C2",
                       linewidths=2, label="Used in ZP")
    ax.set_title(
        f"{frame.path.name}\n"
        f"BJD={row.get('bjd_tdb', np.nan):.6f}, "
        f"V={row.get('v_mag', np.nan):.3f}, "
        f"Ncomp={int(row.get('n_comp', 0))}"
    )
    ax.set_xlabel("x [pixel]")
    ax.set_ylabel("y [pixel]")
    ax.legend(loc="best", fontsize=8)
    save_figure(fig, output_path, dpi=150, pdf=False)


def plot_lightcurve(data, x="phase", phase_start=0.75):
    """Plot the supplied measurements; select accepted frames before calling.

    Phase wrapping is for display only: each measurement appears once in
    [phase_start, phase_start + 1), with its original magnitude and error.
    """
    if x == "phase":
        values = phase_start + (data["phase"].to_numpy() - phase_start) % 1.0
        xlabel = "Orbital phase"
        title = "CZ Aqr — phase-folded light curve"
    elif x == "bjd_tdb":
        values = data[x]
        xlabel = "BJD$_{TDB}$"
        title = "CZ Aqr — calibrated light curve"
    else:
        raise ValueError("x must be 'phase' or 'bjd_tdb'")
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.errorbar(
        values, data["v_mag"], yerr=data["v_mag_err"],
        fmt="o", markersize=4, markerfacecolor="none", capsize=3, elinewidth=1
    )
    ax.invert_yaxis()
    if x == "phase":
        ax.set_xlim(phase_start, phase_start + 1.0)
        ax.set_xticks(np.linspace(phase_start, phase_start + 1.0, 5))
    ax.set_xlabel(xlabel)
    ax.set_ylabel("V [mag]")
    ax.set_title(title)
    ax.text(0.02, 0.02, "Error bars: target SNR + reference-pair scatter (quadrature)",
            transform=ax.transAxes, fontsize=8, color="0.35")
    ax.grid(alpha=0.25)
    return fig, ax


def make_summary_plots(frames, output_dir: Path, phase_start=0.75):
    data = frames[frames["quality_flag"].eq("OK") & np.isfinite(frames["v_mag"])]
    data = data.sort_values("bjd_tdb")
    if data.empty:
        return
    plots = output_dir / "plots"
    plots.mkdir(parents=True, exist_ok=True)
    for x, suffix in [("bjd_tdb", "time"), ("phase", "phase")]:
        fig, _ = plot_lightcurve(data, x=x, phase_start=phase_start)
        save_figure(fig, plots / f"lightcurve_{suffix}.png")


def save_zero_point_plot(frames, comparisons, output_dir):
    """Show frame zero points and the residuals of accepted APASS comparisons."""
    frames = frames[frames["quality_flag"].eq("OK")].sort_values("hjd_utc")
    if frames.empty:
        return
    stars = comparisons[comparisons["used_for_zp"]].merge(
        frames[["filename", "zero_point", "instrument"]], on="filename")
    stars["zp_residual"] = stars["Vmag_catalog"] - stars["inst_mag"] - stars["zero_point"]
    time_origin = int(np.floor(frames["hjd_utc"].min()))
    fig, (time_ax, residual_ax) = plt.subplots(2, 1, figsize=(9, 7))
    for index, (instrument, group) in enumerate(frames.groupby("instrument")):
        color = f"C{index % 10}"
        time_ax.scatter(group["hjd_utc"] - time_origin, group["zero_point"],
                        color=color, s=28, label=str(instrument))
        selected = stars[stars["instrument"].eq(instrument)]
        residual_ax.scatter(selected["Vmag_catalog"], selected["zp_residual"],
                            color=color, s=12, alpha=0.5)
    time_ax.set_title("CZ Aqr — APASS zero-point calibration")
    time_ax.set_xlabel(f"HJD$_{{UTC}}$ − {time_origin}")
    time_ax.set_ylabel("Frame zero point [mag]")
    time_ax.legend(title="Instrument", fontsize=8, ncol=4)
    residual_ax.axhline(0, color="0.4", linestyle="--", linewidth=1)
    residual_ax.set_xlabel("APASS V magnitude")
    residual_ax.set_ylabel("Comparison ZP − frame ZP [mag]")
    residual_ax.set_title("Accepted comparison-star measurements", fontsize=10)
    for ax in (time_ax, residual_ax):
        ax.grid(alpha=0.25)
    plots = Path(output_dir) / "plots"
    plots.mkdir(parents=True, exist_ok=True)
    save_figure(fig, plots / "zero_point.png")


def save_reference_plot(pair, output_dir):
    """Illustrate the exact sample standard deviation used for target errors."""
    mean = pair["difference_mag"].mean()
    scatter = pair["difference_mag"].std(ddof=1)
    time_origin = int(np.floor(pair["bjd_tdb"].min()))
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.axhspan(mean - scatter, mean + scatter, color="0.9", label="Mean ± sample standard deviation")
    ax.axhline(mean, color="0.4", linestyle="--", linewidth=1)
    for instrument, group in pair.groupby("instrument"):
        ax.scatter(group["bjd_tdb"] - time_origin, group["difference_mag"],
                   s=28, label=str(instrument))
    ax.set_xlabel(f"BJD$_{{TDB}}$ − {time_origin}")
    ax.set_ylabel("REF1 − REF2 [mag]")
    ax.set_title(f"CZ Aqr — reference-star difference (N={len(pair)}, σ={scatter:.5f} mag)")
    ax.legend(fontsize=8, ncol=3)
    ax.grid(alpha=0.25)
    plots = Path(output_dir) / "plots"
    plots.mkdir(parents=True, exist_ok=True)
    save_figure(fig, plots / "reference_difference.png")


def save_modeling_plots(
    observations, model_curve,
    output_dir: Path, phase_start=0.75,
):
    """Save the fitted light curve and observed-minus-model residuals.

    Inputs contain accepted observations and a model sampled across the
    displayed phase interval. Error bars are measurement errors only.
    """
    plots = output_dir / "plots"
    plots.mkdir(parents=True, exist_ok=True)
    phase = phase_start + (observations["phase"].to_numpy() - phase_start) % 1.0
    fig, (ax, residual_ax) = plt.subplots(
        2, 1, figsize=(9, 7), sharex=True,
        gridspec_kw={"height_ratios": [3, 1]},
    )
    ax.errorbar(
        phase, observations["v_mag"], yerr=observations["v_mag_err"],
        fmt="o", capsize=3, elinewidth=1, label="Observed",
    )
    ax.plot(model_curve["phase"], model_curve["model_v_mag"],
            color="C3", label="PHOEBE model")
    ax.invert_yaxis()
    ax.set_ylabel("V [mag]")
    ax.set_title("CZ Aqr — light-curve model")
    ax.legend()
    residual_ax.errorbar(
        phase, observations["residual_mag"], yerr=observations["v_mag_err"],
        fmt="o", capsize=3, elinewidth=1,
    )
    residual_ax.axhline(0, color="0.4", linestyle="--", linewidth=1)
    residual_ax.set_ylabel("O − C [mag]")
    residual_ax.set_xlabel("Orbital phase")
    residual_ax.set_xlim(phase_start, phase_start + 1.0)
    residual_ax.set_xticks(np.linspace(phase_start, phase_start + 1.0, 5))
    for panel in (ax, residual_ax):
        panel.grid(alpha=0.25)
    save_figure(fig, plots / "lightcurve_model.png")

    fig, ax = plt.subplots(figsize=(9, 5))
    for instrument, group in observations.groupby("instrument", dropna=False):
        ax.errorbar(
            group["hjd_utc"], group["residual_mag"], yerr=group["v_mag_err"],
            fmt="o", capsize=3, elinewidth=1, label=str(instrument),
        )
    ax.axhline(0, color="0.4", linestyle="--", linewidth=1)
    ax.set_xlabel("HJD$_{UTC}$")
    ax.set_ylabel("O − C [mag]")
    ax.set_title("CZ Aqr — model residuals by instrument")
    ax.legend(title="Instrument", fontsize=8, ncol=2)
    ax.grid(alpha=0.25)
    save_figure(fig, plots / "residuals_time.png")
