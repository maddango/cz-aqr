# CZ Aqr light curve

Three Python modules make up the project:

- `photometry.py`: aperture photometry, APASS calibration and uncertainties.
- `model.py`: PHOEBE light-curve fitting.
- `plots.py`: image overlays, light curves and model residual plots.

## Run

Run these commands from the project directory:

```bash
# Measure the reduced FITS images.
.venv/bin/python photometry.py inputs --output output

# Fit the light curve.
.venv-phoebe/bin/python model.py
```

Use `config.json` for photometry, reference stars and ephemeris; `model_config.json` contains the adopted physical parameters and fit settings. Both scripts read these fixed filenames from the project directory. Command-line options select input and output locations. Dependencies are listed in `requirements.txt` and `requirements-modeling.txt`. These commands assume the reduced images and both Python environments are already available; see the guide's [reproduction prerequisites](docs/guide.md#reproducing-the-analysis) for details.

## Results

Photometry saves three CSV tables in `output/tables/`, light curves, a zero-point calibration figure and the reference-star difference plot in `output/plots/`. These figures are saved as PNG and PDF. Optional image overlays remain PNGs in `output/diagnostics/`. The PHOEBE fit saves parameters, model predictions, residuals, plots and a short `summary.json` in `output/modeling/`. Check the optimizer status before using the fit. `tables/report_parameters.csv` presents rounded fitted, adopted, fixed and model-derived quantities. `tables/parameters.csv` retains the full-precision fit; `tables/derived_results.csv` records model eclipse depths, brightness reductions, the fitted eclipse epoch and timing offset, and conditional secondary temperature and primary radius. Phase plots retain the adopted ephemeris; sources and the small working-period approximation are documented in `config.json` and the guide.

See [docs/guide.md](docs/guide.md) for the project tree, reproduction prerequisites, methods, scientific results and limitations. The report is maintained in `report/`, which is excluded from version control. The notebooks are optional and are not part of the pipeline. Earlier versions are preserved in `backup/`.

Generate the illustrative calibration plot and field image with aperture inset from the saved photometry and one retained FITS image:

```bash
.venv/bin/python scripts/plot_photometry_example.py
```

This writes `zero_point_example` and `field_example` as PNG and PDF files in `output/plots/`.

Run the measurement checks with:

```bash
.venv/bin/python -m unittest discover -s tests -v
```
