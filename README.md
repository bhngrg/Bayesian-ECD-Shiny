# Bayesian-ECD Shiny App

This repository contains the Bayesian-ECD Shiny application, a local graphical interface for applying the Bayesian-ECD historical-borrowing workflow to a current randomized clinical trial (RCT). The application provides posterior survival and treatment-effect estimation, restricted mean survival time (RMST) summaries, prespecified subpopulation analyses, prediction in new patient populations, and an optional concurrent-control compatibility diagnostic.

The application runs locally in R and uses stored posterior information from the historical Bayesian-ECD model. Current RCT data therefore do not need to be uploaded to a separately hosted Bayesian-ECD service.

## User Guide

For complete installation instructions, platform-specific setup, data requirements, analysis workflows, interpretation guidance, output descriptions, reproducibility recommendations, and troubleshooting, see:

- [Bayesian-ECD Shiny User Guide](docs/Bayesian-ECD-Shiny-User-Guide.md)

A rendered PDF version of the user guide is provided separately with the journal supplementary materials.

The user guide is the primary reference for using the application. The example below is intended only as a quick orientation to the workflow.

## Quick Example

### 1. Launch the application

After completing the installation steps in the user guide, open:

```text
app/app.R
```

in RStudio and select:

```text
Run App > Run in Window
```

Alternatively, from the repository root, run:

```r
shiny::runApp("app")
```

The **Uploaded Data** tab is the starting point for the primary Bayesian-ECD workflow.

### 2. Upload the example current-RCT dataset

Use:

```text
example_data/example_model_data.csv
```

Upload the file in the **Uploaded Data** tab and enter the requested column specifications.

The example dataset contains the information required by the Bayesian-ECD analysis, including patient identifiers, survival outcomes, censoring indicators, treatment and cohort information, and the baseline covariates used by the model.

Detailed definitions, coding requirements, and input-format guidance are provided in the user guide.

### 3. Optionally run the control compatibility diagnostic

If the uploaded dataset contains a concurrent control arm labeled:

```text
Control
```

the **Control Compatibility** tab can be used to compare the observed concurrent-control survival experience with the historical Bayesian-ECD posterior predictive control distribution.

This diagnostic is optional and is separate from the primary Stage 2 Bayesian-ECD model extension. A concurrent control arm is not required to run the primary Bayesian-ECD analysis.

### 4. Review the primary Bayesian-ECD outputs

After fitting the model, the main analysis results are available through:

- **Plot Output** for posterior survival and hazard-ratio summaries;
- **RMST Output** for restricted mean survival time comparisons;
- **Subpopulation Analysis** for survival or hazard-ratio summaries within a selected subpopulation; and
- **Subpopulation RMST Output** for RMST comparisons within a selected subpopulation.

The user guide describes the available treatment selections, plot types, interpretation, and download options in detail.

### 5. Upload the example prediction dataset

To evaluate treatment-specific survival in a new covariate distribution, use:

```text
example_data/example_prediction_data.csv
```

The prediction dataset contains patient identifiers and the baseline covariates required by the stored Bayesian-ECD model. It does not require observed survival outcomes, censoring indicators, treatment assignments, or cohort labels.

Prediction results are available through:

- **Prediction Output**; and
- **Prediction RMST Output**.

Optional subpopulation filters can also be applied to the prediction population.

### 6. Optional posterior-probability analyses

Additional posterior-probability utilities are available outside the Shiny interface for analyses such as:

- posterior probabilities for prespecified hazard-ratio thresholds;
- posterior probabilities for RMST-difference thresholds; and
- posterior probabilities for RMST-ratio thresholds.

The utility functions are located in:

```text
analysis_utils/posterior-probability-utils.R
```

A worked example is provided in:

```text
examples/run_posterior_probabilities.R
```

From the repository root, the example can be run with:

```r
source("examples/run_posterior_probabilities.R")
```

Generated posterior-probability results are written to the Git-ignored:

```text
outputs/
```

directory.

## Repository Components

The main repository directories are:

- `app/`: Shiny application files.
- `app_core/`: Core modeling, plotting, RMST, subpopulation, and C++ routines used by the application.
- `app/data/storedMCMCiter/`: Stored Stage 1 Bayesian-ECD posterior components used by the application.
- `analysis_utils/`: Optional analysis utilities that are not part of the Shiny interface.
- `examples/`: Worked R examples for optional analyses.
- `example_data/`: Example model-fitting and prediction datasets.
- `docs/`: Bayesian-ECD Shiny User Guide documentation.

## Reproducibility

The stored Stage 1 posterior components under:

```text
app/data/storedMCMCiter/
```

are part of the computational specification of the application. Reproducible analyses should identify the repository version used and retain the corresponding input datasets, analysis settings, and outputs.

See the [Bayesian-ECD Shiny User Guide](docs/Bayesian-ECD-Shiny-User-Guide.md) for the complete reproducibility recommendations.
