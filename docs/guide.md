# CZ Aqr project guide

This project measures a calibrated Johnson V light curve of CZ Aqr from reduced LCO images and fits an adopted semi-detached binary model with PHOEBE. The analysis runs through `photometry.py`, `model.py`, and `plots.py`; notebooks and archived versions provide supplementary and historical material.

## Important project files

```text
cz_aqr_photometry_project/
├── README.md                         # commands and quick overview
├── docs/guide.md                     # methods, results, interpretation
├── photometry.py                     # aperture photometry, calibration and errors
├── model.py                          # PHOEBE fit and derived scientific results
├── plots.py                          # PNG/PDF figures and image overlays
├── scripts/plot_photometry_example.py # illustrative calibration and field figures
├── config.json                       # target, ephemeris, apertures, references and quality cuts
├── model_config.json                 # adopted physics, initial values and fit bounds
├── requirements.txt
├── requirements-modeling.txt
├── inputs/                           # reduced FITS.FZ images
├── rejected_inputs/                  # locally excluded FITS images
├── reference/apass_comparison_catalog.ecsv
├── output/
│   ├── tables/
│   │   ├── lightcurve_all_frames.csv
│   │   ├── comparison_star_photometry.csv
│   │   └── reference_pair.csv
│   ├── plots/                        # each chart has PNG and PDF versions
│   │   ├── lightcurve_time.pdf
│   │   ├── lightcurve_phase.pdf
│   │   ├── zero_point.pdf
│   │   ├── zero_point_example.pdf
│   │   ├── field_example.pdf
│   │   └── reference_difference.pdf
│   ├── diagnostics/                  # source-position overlays, PNG only
│   └── modeling/
│       ├── summary.json              # convergence and fit statistics
│       ├── tables/
│       │   ├── parameters.csv         # full-precision fitted parameters
│       │   ├── report_parameters.csv  # rounded fitted, adopted, fixed and derived values
│       │   ├── derived_results.csv    # eclipse depths, epoch, timing offset, conditional T2/R1
│       │   ├── model_observations.csv # observations, predictions and residuals
│       │   └── model_curve.csv        # sampled model used for curves and eclipse depths
│       └── plots/
│           ├── lightcurve_model.pdf
│           └── residuals_time.pdf
├── report/                           # local report sources and PDF
├── tests/                            # numerical measurement checks
└── backup/                           # historical versions
```

## Reproducing the analysis

The commands in [README.md](../README.md#run) run photometry first and then fit the resulting light curve. They assume the reduced LCO FITS/FITS.FZ images are available in `inputs/`, with `SCI`, `ERR`, and `BPM` extensions. The images and Python environments are excluded from Git, so a repository clone alone does not contain everything needed for a full rerun. The fixed APASS catalogue is included at `reference/apass_comparison_catalog.ecsv`.

Photometry uses `.venv` with dependencies from `requirements.txt`; modeling uses `.venv-phoebe` with dependencies from `requirements-modeling.txt`, which pins PHOEBE to 2.5.4. The model also requires locally installed `Johnson:V` passband tables with CK2004 atmosphere and interpolated limb-darkening data. `model.py` disables online passband access by default, so a new environment needs these data installed before fitting. The dependency files do not pin every package version and therefore do not fully specify the environment used for the saved results.

Both commands run from the project directory, where they read `config.json` and `model_config.json`. Photometry writes `output/tables/lightcurve_all_frames.csv`, which is the model's default input. To compare a rerun with the saved results, use a separate photometry output directory with `--output`, then pass its light-curve CSV to `model.py --input` and select a separate modeling directory with `--output`. Each run writes its tables and figures to the selected location, replacing files with the same names.

## Current dataset and fit

The current photometry and PHOEBE results use the 30 images retained after the quality assessment below. The numerical results are saved in `output/tables/` and `output/modeling/`.

- 30 accepted V-band images, each exposed for 25 seconds, from eight instrument identifiers. The fixed catalogue has 15 stars; 13–15 contribute per image.
- The image headers span **2026-10-04 18:01:58 to 2026-10-09 05:22:15 UTC**; the mid-exposure HJD span is 4.4722 days. Filenames can refer to the preceding observing night, so use `date_obs` for calendar dates.
- Reference-pair sample standard deviation: **0.00834 mag**; median final target uncertainty: **0.00865 mag**.
- The fit converged after **385 evaluations** using 30 observations and five fitted quantities (four shape parameters plus the magnitude offset).
- Chi-squared = **22.6116** for **25 degrees of freedom**; reduced chi-squared = **0.9045**; residual RMS = **0.00753 mag**.
- Fitted inclination = **88.79 degrees**, primary radius/separation = **0.3553**, temperature ratio = **0.6725**, phase shift = **−0.00411 cycles**. None is at an allowed bound. Parameter confidence intervals have not been estimated.

### Input selection and excluded images

The quality review on **2026-10-28** established a selection of **30 out of 33 images**, excluding the three exposures below because their comparison-star measurements are inconsistent with the adopted single-zero-point calibration. All three pass the automatic measurement thresholds; their exclusion follows the additional calibration assessment.

| Excluded image | Comparisons used | Comparison zero-point scatter | Main finding |
| --- | ---: | ---: | --- |
| `coj0m416-sq36-20261005-0053-e91.fits.fz` | 6 | 0.050 mag | Reference-pair difference displaced by −0.100 mag from the baseline median; individual comparison stars also disagree. |
| `elp0m414-sq31-20261006-0143-e91.fits.fz` | 4 | 0.113 mag | Strong disagreement even among the bright comparisons; recalibration after omitting one comparison can shift CZ Aqr by 0.049 mag. |
| `elp0m414-sq31-20261006-0121-e91.fits.fz` | 13 | 0.047 mag | Discrepancies across several comparison stars; recalibration after omitting one comparison can shift CZ Aqr by 0.028 mag despite the larger ensemble. |

The comparison zero-point scatter is the sample standard deviation of `V_APASS - m_instrumental` among the stars used in each image's calibration. Its median in the other 30 images is **0.010 mag**. It measures agreement within an image and is distinct from the temporal reference-pair scatter used for the light-curve errors. Checks that calibrated each comparison star without including it in its own zero point confirmed discrepancies across the field. In particular, the problem in `...0121` extends beyond the fixed reference pair, supporting exclusion of the whole observation. The number of comparisons alone does not explain these decisions.

The target and fixed-reference apertures show no obvious saturation, flagged pixels or gross centering problems, and the anomalous reference differences persist under modest aperture changes. The spatial pattern and reduced signals are compatible with nonuniform transparency, although the physical cause is not established. Exclusion of `...0053` is conservative because no other exposure from its instrument is available to separate observing-condition effects from an instrument-dependent response.

The three excluded images are stored in `rejected_inputs/`; the 30 retained images in `inputs/` supply the photometry, reference-pair scatter and model fit reported here. The selection is based on calibration consistency.

## Photometry

`photometry.py` reads reduced LCO images containing `SCI` (image), `ERR` (pixel uncertainty), and `BPM` (bad-pixel mask). It measures CZ Aqr and the 15 stars in `reference/apass_comparison_catalog.ecsv`.

The steps are:

1. Convert the sky coordinates to pixels and refine the centroids.
2. Sum the pixels in a circular aperture; estimate the sky from a sigma-clipped mean in a surrounding annulus and subtract it.
3. Calculate the flux uncertainty from `ERR` and the uncertainty in that sky estimate, then calculate `SNR = flux / flux_error` and `m_instrumental = -2.5 log10(flux / exposure_seconds)`.
4. Reject unusable measurements and calculate a weighted APASS zero point.
5. Add the zero point to CZ Aqr's instrumental magnitude.
6. Calculate the two-reference uncertainty and save the tables and plots.

`config.json` sets aperture and annulus sizes in multiples of the image FWHM. It also sets the minimum SNR, saturation threshold and minimum comparison count. Bad stellar apertures are rejected; bad annulus pixels are excluded from the sky estimate. `quality_flag == "OK"` selects accepted target observations.

The current thresholds are target SNR ≥ 20, comparison SNR ≥ 30, peaks below 90% of the detector linearity limit, and at least three retained comparison stars. Measurements also require valid fluxes and errors and sufficient room for the aperture and sky annulus. `OK` means these configured checks passed; it does not establish that all calibration effects are negligible. Aperture sizes use the header seeing estimate when it lies between 0.7 and 8.0 arcsec and the pixel scale is valid; otherwise the configured FWHM fallback is 4 pixels. If centroid refinement fails, the measurement retains the WCS position, making the diagnostic overlays useful for checking source placement.

For each accepted comparison, the zero point is `V_APASS - m_instrumental`. After sigma clipping, the mean is weighted by `1 / (instrumental_mag_error² + APASS_mag_error²)`. The fixed catalogue is APASS DR9 (`II/336/apass9`); its original selection kept stars with V between 10.5 and 15, catalogue errors at most 0.08 mag, and separation at least 20 arcsec from CZ Aqr. No online catalogue query is needed.

## Uncertainties: the professor's procedure

The same two reference stars are used across the observations:

| Reference | ID (RA and Dec in degrees) | APASS V |
| --- | --- | --- |
| REF1 | `350.7266070_-15.8521530` | 10.820 |
| REF2 | `350.5999100_-15.9874830` | 10.697 |

Their IDs are set in `config.json` under `uncertainty.reference_star_ids`. They are close to CZ Aqr in brightness and are assumed constant.

```text
difference_i = m_REF1,i - m_REF2,i
instrumental_err = std(difference_i, ddof=1)
target_inst_mag_err_i = 1.085736 / SNR_target,i
v_mag_err_i = sqrt(target_inst_mag_err_i² + instrumental_err²)
```

The standard deviation uses accepted frames where both reference measurements pass the quality checks; at least three paired observations are required. It is not divided by sqrt(N), and pair differences are not sigma-clipped. The common zero point cancels in their difference, so instrumental magnitudes are sufficient. Zero-point clipping does not remove a star from this calculation.

The full reference scatter supplies one instrumental term for the run. The approximately 15 APASS stars still calibrate magnitudes, but no additional ensemble-error term is added. Both plots and PHOEBE use `v_mag_err`.

This is a simple empirical approximation: the scatter includes reference-star noise and differential observing effects. The input-quality assessment above identifies three exposures whose calibration discrepancies make them unsuitable for the adopted reduction and shared uncertainty prescription.

The adopted errors do not constitute a complete systematic uncertainty budget. A constant APASS calibration offset cancels in the reference-star difference, so that scatter does not measure absolute calibration uncertainty. Shared catalogue errors, colour-dependent extinction and instrument response, residual detector effects, contamination, and aperture or centroid choices are not separately modeled. Some differential effects may contribute to the measured reference scatter, but their contributions are not isolated and the fit treats the final errors as independent.

`output/plots/reference_difference.png` plots these same differences against BJD (TDB), colored by instrument. The line marks their mean, **0.11432 mag**; the shaded band is mean ± sample standard deviation, **0.00834 mag**, using all 30 retained pairs. This band is not the uncertainty of the mean. Reference IDs are recorded in the pair table.

## Tables, timing and phase

`output/tables/` contains:

- `lightcurve_all_frames.csv`: times, phase, calibrated V magnitude and error, target flux and SNR, zero point, comparison count and quality flag.
- `comparison_star_photometry.csv`: each star's instrumental magnitude, error, SNR, peak, aperture flag and whether it contributed to the zero point.
- `reference_pair.csv`: the two reference magnitudes and their differences, with identifiers, times and instrument.

The `_counts` suffix on flux, peak and detector-limit columns refers to the processed image's count scale. The calculations use matching SCI/ERR image units without an extra gain conversion. Times refer to mid-exposure. Orbital phase is `((hjd_utc - t0_hjd) / period_days) % 1` using the configured HJD ephemeris. `bjd_tdb` is saved separately. The plots wrap phases into [0.75, 1.75) by default, so each observation appears once. The period is adopted, not fitted.

The fixed plotting ephemeris uses `T0 = 2460232.18988` and `P = 0.86275208 days`. The epoch matches the final TESS-derived minimum in [Zeng et al. (2024), Table 6](https://arxiv.org/html/2411.06803v1#S0.T6), listed as `2460232.18988(09)` under HJD. The table does not specify UTC/TDB; HJD_UTC remains the project's working convention, not a verified time-scale statement from that source. [AAVSO VSX, CZ Aqr](https://vsx.aavso.org/index.php?view=detail.top&oid=971) (checked 2026-10-27) and Zeng et al., Eq. 1, give `P = 0.86275209 days`. The retained working period differs by only `-0.00000001 days`, accumulating to about 1.1 seconds over 1262 cycles. This approximation preserves the existing phases and fit. We combine a recent tabulated minimum with that near-identical period; this is not the historical reference epoch used in the paper's Eq. 1. The model-derived epoch below is an output and is not fed back into phase plots.

`output/plots/zero_point.png` shows the frame zero point over time, colored by instrument. Its lower panel shows `V_APASS - m_instrumental - frame_zero_point` against APASS magnitude for each comparison measurement used in calibration. These residuals illustrate agreement within the ensemble; they are not the light-curve error bars. Frame zero points can change with instrument throughput and observing conditions, so their variation is not target variability.

Each plotted comparison contributes to its own frame zero point, so these residuals are not an independent test of that star's temporal stability. Such a test would calibrate the star with an ensemble that excludes it.

All light-curve, calibration, reference-difference and model figures also have PDF copies for the report. The image overlays in `output/diagnostics/` remain PNGs because their backgrounds are raster images.

### Illustrative photometry figures

`output/plots/zero_point_example.pdf` and `field_example.pdf` illustrate the same retained exposure, `cpt0m439-sq39-20261004-0106-e91.fits.fz` (2026-10-04, 19:30:16 UTC), with 15 comparison stars. Both also have PNG copies. The calibration example shows catalogue versus instrumental magnitudes with the adopted slope-one relation and `ZP = 20.7336 mag`, followed by the individual zero points. Its shaded band is the weighted mean zero point ± the comparison sample standard deviation, **0.0107 mag**; it represents calibration scatter, not the final target uncertainty.

The field image identifies CZ Aqr, the APASS comparisons and the two fixed reference stars. The large field markers identify sources; only the enlarged inset shows the actual stellar aperture and sky annulus, with radii of 1.5 FWHM and 2.5–4 FWHM, respectively. These two figures can be generated from the saved measurements and original image with `.venv/bin/python scripts/plot_photometry_example.py`, without rerunning photometry or fitting.

## PHOEBE model

`model.py` uses a circular, synchronous, semi-detached binary whose secondary fills its Roche lobe. The adopted mass ratio (0.4866), primary temperature (8200 K) and separation (5.09 solar radii) are recorded in `model_config.json`, following [Zeng et al. (2024), Tables 2–3](https://arxiv.org/html/2411.06803v1#S3). These are inputs, not measurements made by this project.

The model uses Johnson V, CK2004 atmospheres, interpolated limb darkening, Horvat irradiation, and 3000 surface triangles per star. It evaluates exposure midpoints. There are no spots or third light.

Four parameters are fitted: inclination, primary radius divided by separation, secondary/primary temperature ratio, and phase shift. Initial values and bounds are in `model_config.json`. A positive phase shift moves primary eclipse later. A common magnitude offset is solved analytically at each trial, and Powell optimization minimizes `sum(((observed - model) / v_mag_err)²)`.

The model output contains:

- `tables/parameters.csv`: the four shape parameters and magnitude offset.
- `tables/report_parameters.csv`: rounded fitted parameters, adopted inputs, fixed assumptions and model-derived results, with their status labelled.
- `tables/derived_results.csv`: full-precision derived values with units, status and the definition of each quantity.
- `tables/model_observations.csv`: observations, predictions and residuals.
- `tables/model_curve.csv`: the model sampled over one displayed orbital cycle.
- `plots/lightcurve_model.png` and `plots/residuals_time.png`.
- `summary.json`: optimizer status, evaluation count, chi-squared and RMS.

Every run fits the observations. Set `optimizer.max_evaluations` in `model_config.json` to limit fit work. Only the final fit summary is printed; if optimization stops early, parameters are labelled `trial` and the command exits with code 2 after saving the results. A converged run is labelled `fitted`. The offset is labelled `profiled` in either case. Rounding in the report table is only for presentation; it does not imply parameter confidence intervals. The full-precision fit table remains unchanged.

These are exploratory fits with approximate independent errors. Optimizer convergence indicates that its stopping criteria were met; it does not establish a unique physical solution. Parameter confidence intervals and instrument-dependent offsets are not fitted. Masses, absolute radii and distance are not independently measured by this project. The secondary temperature below is conditional on the adopted primary temperature. The absolute primary radius is likewise conditional on the adopted separation.

## Derived scientific results

`model.py` calculates these quantities from the fitted parameters and sampled model curve and writes them to `output/modeling/tables/derived_results.csv` after each fit.

### Eclipse depths and brightness reductions

Both eclipse depths use the brightest point of the **sampled fitted curve** as a common reference:

```text
V_bright = min(model_V)
depth = V_faintest_in_eclipse - V_bright
flux_loss_percent = 100 * (1 - 10**(-0.4 * depth))
```

The primary half-orbit is centered on the fitted phase shift; the secondary half-orbit is centered half a cycle later. Phase differences are wrapped so the selection does not depend on how the graph displays phase. Within each half-orbit the largest magnitude is the faintest sampled model point. Both depths use the same brightest-model reference, rather than a fitted local continuum at eclipse ingress. They therefore describe the full change from maximum model brightness, including the variation outside eclipse.

| Derived quantity | Current value | Interpretation |
| --- | ---: | --- |
| Brightest model V magnitude | 10.660 mag | Reference level for both depths |
| Primary eclipse depth | 1.432 mag | Model-derived depth |
| Primary eclipse flux loss | 73.2% | Relative to brightest model flux |
| Secondary eclipse depth | 0.258 mag | Model-derived depth |
| Secondary eclipse flux loss | 21.1% | Relative to brightest model flux |

The curve has 301 samples across one cycle, with spacing about 0.00333 cycles. The extrema are sampled values, not exact continuous extrema. In particular, the bottom of primary eclipse is not directly sampled by the observations; its depth comes from the fitted model. The flux-loss percentages measure the reduction in total system brightness, which also depends on the stars' surface brightnesses and cannot be identified with the fraction of stellar surface area covered.

### Conditional secondary temperature

```text
T2 = fitted_temperature_ratio * adopted_T1
   = 0.6725081667 * 8200 K
   = 5514.57 K  (about 5515 K for the report)
```

This converts the fitted ratio into a temperature using an external input. It is not an independent temperature measurement or a spectroscopic result. No temperature uncertainty has been propagated. The adopted 8200 K value and its literature basis are identified in the PHOEBE section above.

### Conditional primary radius

The primary's volume-equivalent radius is `R1 = fitted_R1/a * adopted_a`: `0.3553371180 * 5.09 = 1.80867 solar radii`, or **1.809 solar radii** for the report. This is the radius used by PHOEBE, not a polar or side radius of the distorted star. It depends on the adopted separation and is not an independent absolute-size measurement. No radius uncertainty has been propagated.

### Model-derived primary-eclipse epoch

The fitted phase shift determines a primary conjunction time. The calculation selects the conjunction nearest the median observation HJD, placing the reported epoch near this dataset:

```text
E = round((median(observation_HJD) - adopted_T0) / P - fitted_phase_shift)
T_primary = adopted_T0 + (E + fitted_phase_shift) * P
timing_offset_minutes = fitted_phase_shift * P * 1440
```

For the saved fit, **E = 1262**, **T_primary = 2461320.97946 HJD_UTC**, and the timing offset is **-5.11 minutes** (earlier than the adopted linear ephemeris). This is a conjunction inferred from the combined phase-folded light curve; it need not coincide with a directly observed eclipse. In our circular model, it defines the primary-eclipse center. The epoch is calculated directly from the fitted shift, independently of the 301-point plotting grid.

The period remains fixed, and the adopted ephemeris still controls every phase plot. There is no additional fitted parameter or change to fit degrees of freedom. No timing confidence interval has been calculated: sparse eclipse sampling and adopted geometry limit accuracy. The reported decimals are computational precision. The offset is not evidence of period change, and its physical interpretation also requires a consistent external time convention.

All these quantities are marked `model-derived` after a converged fit, or `trial-derived` if the optimizer stops early. Rounding in the report table is for presentation and does not supply missing uncertainty estimates.

## Data products and academic report

The report LaTeX sources and compiled PDF stay in `report/`, which is excluded from Git.

- **Methods and inputs:** `photometry.py`, `model.py`, `plots.py`, both JSON configuration files, and `reference/apass_comparison_catalog.ecsv`.
- **Measurements:** the three CSV files in `output/tables/`.
- **Fit and derived results:** `output/modeling/summary.json` and the CSV files in `output/modeling/tables/`. The full-precision tables support further calculations; `report_parameters.csv` provides rounded values for presentation.
- **Figures:** the PDFs in `output/plots/` and `output/modeling/plots/`. Local PNG overlays in `output/diagnostics/` can illustrate aperture placement.

These files connect the methods described in this guide to the measurements, fitted parameters and figures used in the report. The saved tables support inspection and further analysis without repeating the image measurements; reproducing the photometry itself requires the original reduced FITS images.
