# Bayesian-ECD Shiny App

This repository contains a Shiny application for Bayesian-ECD survival modeling, subgroup analysis, restricted mean survival time (RMST) summaries, and prediction on new covariate-only datasets.

## Repository structure

- `app/`: Shiny application files.
- `app_core/`: Core modeling, plotting, RMST, subgroup, and C++ helper routines used by the app.
- `app/data/storedMCMCiter/`: Stored baseline MCMC components loaded by the app.
- `analysis_utils/`: Additional user-facing analysis utilities that are not part of the Shiny interface.
- `examples/`: Example R scripts showing how to run optional analyses outside the Shiny app.
- `example_data/`: Example model-fitting and prediction datasets.
- `notes/`: Development notes and dependency map.

## First-time R/RStudio setup

If you are new to R, install R first, then install RStudio Desktop.

1. Download and install R from CRAN:  
   https://cran.r-project.org/

2. Download and install RStudio Desktop from Posit:  
   https://posit.co/download/rstudio-desktop/

After installing both, open RStudio and run the following in the RStudio Console:

```r
R.version.string
```

This should print your installed R version.

## Where do I run commands?

This README uses two types of commands.

### R commands

Commands marked as `r` should be run in the RStudio Console.

Example:

```r
install.packages("shiny")
```

### Terminal or bash commands

Commands marked as `bash` should be run in your system terminal.

- macOS: use the Terminal app.
- Linux: use Terminal.
- Windows: use Git Bash, PowerShell, or the terminal inside RStudio.

Example:

```bash
cd path/to/project-folder
```

If you are new to command-line tools, downloading the repository as a ZIP file is simpler than cloning with Git.

## Download this repository

You can get the repository onto your computer using either the ZIP download method or Git.

### Option A: Download ZIP from GitHub

This is the simplest option for first-time users.

1. Go to the GitHub repository page.
2. Click the green **Code** button.
3. Click **Download ZIP**.
4. Unzip the downloaded file.
5. Open RStudio.
6. In RStudio, go to **File > Open Project** if an `.Rproj` file is available, or use **Session > Set Working Directory > Choose Directory** and select the unzipped repository folder.

### Option B: Clone with Git

If you use Git, run the following in Terminal or Git Bash:

```bash
git clone https://github.com/USERNAME/REPOSITORY.git
cd REPOSITORY
```

Replace `USERNAME` and `REPOSITORY` with the actual GitHub username and repository name.

## Install dependencies

This app uses a mixture of regular R packages and C++-backed R packages. Most dependencies install directly from CRAN, but `RcppGSL` may require an additional system library called the GNU Scientific Library (GSL).

### Step 1: Install system requirements for compiled R packages

Some R packages need compilation. Make sure your system has the usual R build tools installed.

#### macOS

Install Xcode command line tools by running this in Terminal:

```bash
xcode-select --install
```

If you use Homebrew, install GSL:

```bash
brew install gsl
```

#### Ubuntu/Debian Linux

Install system build tools and GSL by running this in Terminal:

```bash
sudo apt-get update
sudo apt-get install build-essential libcurl4-openssl-dev libssl-dev libxml2-dev libgsl-dev
```

#### Windows

Install Rtools for your version of R:

https://cran.r-project.org/bin/windows/Rtools/

After installing Rtools, restart RStudio before installing packages.

### Step 2: Install regular R package dependencies

Run this in the RStudio Console:

```r
install.packages(c(
  "shiny",
  "bslib",
  "tibble",
  "dplyr",
  "readr",
  "ggplot2",
  "reshape2",
  "plyr",
  "foreach",
  "doParallel",
  "matrixStats",
  "DT",
  "patchwork",
  "gridExtra",
  "gt",
  "writexl",
  "zip",
  "Rcpp",
  "RcppArmadillo",
  "RcppDist",
  "purrr",
  "rstudioapi"
))
```

### Step 3: Check whether `RcppGSL` can be installed

Before installing all remaining dependencies at once, test `RcppGSL` separately.

Run this in the RStudio Console:

```r
install.packages("RcppGSL")
library(RcppGSL)
```

If this works, continue to the next section.

If this fails, the issue is usually that GSL or the system compiler tools are missing. Check the error message and confirm that the system setup steps above were completed for your operating system.

### Step 4: Confirm package loading

After installation, run this in the RStudio Console:

```r
required_packages <- c(
  "shiny",
  "bslib",
  "tibble",
  "dplyr",
  "readr",
  "ggplot2",
  "reshape2",
  "plyr",
  "foreach",
  "doParallel",
  "matrixStats",
  "DT",
  "patchwork",
  "gridExtra",
  "gt",
  "writexl",
  "zip",
  "Rcpp",
  "RcppArmadillo",
  "RcppDist",
  "RcppGSL",
  "purrr",
  "rstudioapi"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) == 0) {
  message("All required R packages are installed.")
} else {
  message("Missing packages: ", paste(missing_packages, collapse = ", "))
}
```

## Run the app

Open RStudio and make sure your working directory is the repository root folder.

You can check your current working directory by running:

```r
getwd()
```

If needed, set the working directory manually. For example:

```r
setwd("path/to/your/repository-folder")
```

Then run:

```r
shiny::runApp("app")
```

Alternatively, open `app/app.R` in RStudio and click **Run App**.

## Example workflow

### 1. Upload model-fitting data

Use the example model dataset in:

```text
example_data/
```

In the **Uploaded Data** tab, enter the requested column names.

The model-fitting dataset should contain:

- patient ID
- overall survival time
- censoring indicator
- treatment
- cohort/source indicator
- covariates used for modeling, such as age, sex, KPS, and extent of resection

### 2. View model outputs

Use the following tabs:

- **Plot Output**
- **RMST Output**
- **Subgroup Analysis**
- **Subgroup RMST Output**

### 3. Upload prediction data

Use the example prediction dataset in:

```text
example_data/
```

The prediction dataset should contain:

- prediction patient ID
- the same covariates used in model fitting

It does not need outcome, censoring, treatment, or cohort/source columns.

### 4. View prediction outputs

Use the following tabs:

- **Prediction Output**
- **Prediction RMST Output**

The prediction workflow reruns `cappmx_extend_approx_fit()` with `input_df_pred = pred_data()` and then uses `subgroup_data(..., use_pred = TRUE)`.

## Optional posterior probability utilities

The Shiny app focuses on survival curves, hazard-ratio curves, RMST summaries, subgroup analyses, and prediction outputs.

Additional posterior probability utilities are provided in:

```text
analysis_utils/posterior-probability-utils.R
```

An example script is provided in:

```text
examples/run_posterior_probabilities.R
```

This script shows how to compute quantities such as:

- `P(HR(t) < threshold | Data)`
- `P(RMST difference > threshold | Data)`
- `P(RMST ratio > threshold | Data)`

The script can be opened in RStudio and run directly. It automatically locates the repository root when sourced from RStudio.

Generated posterior probability outputs are saved to:

```text
outputs/
```

The `outputs/` folder is ignored by Git because it contains user-generated results.

## Cross-platform notes

The app is organized to run on Windows, macOS, and Linux.

Parallel computations use `foreach`, `doParallel`, and `parallel::makeCluster()`, which use PSOCK-style workers and are more portable than Unix-only fork-based parallelization.

R-level density and survival helper functions use `matrixStats::logSumExp()` rather than the Rcpp `log_sum_exp()` function inside parallel workers to avoid external-pointer issues.

## Development notes

See:

- `notes/work_log.txt`
- `notes/app_core_dependency_map.txt`
