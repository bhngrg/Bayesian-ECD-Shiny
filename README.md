# Bayesian-ECD Shiny App

This repository contains the Bayesian-ECD Shiny application, which provides a local graphical interface for applying the Bayesian-ECD historical-borrowing workflow to a current randomized clinical trial (RCT). The application provides posterior survival and treatment-effect estimation, restricted mean survival time (RMST) summaries, prespecified subpopulation analyses, and prediction in new patient populations.

The application runs locally in R and uses stored posterior information from the historical Bayesian-ECD model, so current RCT data do not need to be uploaded to a separately hosted Bayesian-ECD service.

## User guide

Detailed installation instructions, data requirements, analysis workflows, output descriptions, and troubleshooting guidance are provided in:

```text
docs/Bayesian-ECD-Shiny-User-Guide.Rmd
```

A rendered PDF version of the user guide is provided separately with the journal supplementary materials.

## Repository structure

- `app/`: Shiny application files.
- `app_core/`: Core modeling, plotting, RMST, subpopulation, and C++ helper routines used by the app.
- `app/data/storedMCMCiter/`: Stored baseline MCMC components loaded by the app.
- `analysis_utils/`: Additional user-facing analysis utilities that are not part of the Shiny interface.
- `examples/`: Example R scripts showing how to run optional analyses outside the Shiny app.
- `example_data/`: Example model-fitting and prediction datasets.
- `docs/`: Source for the Bayesian-ECD Shiny supplementary user manual.

## First-time R/RStudio setup

The application requires R version 4.2.2 or later. If you are new to R, install R first, then install RStudio Desktop.

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
git clone https://github.com/bhngrg/Bayesian-ECD-Shiny.git
cd Bayesian-ECD-Shiny
```

## Install dependencies

This app uses a mixture of regular R packages and C++-backed R packages. All listed R package dependencies are available from CRAN. The C++ code is compiled through `Rcpp`, `RcppArmadillo`, and `RcppDist`.

### Step 1: Install system requirements for compiled R packages

Most users will receive precompiled CRAN binaries through `install.packages()`, but the system requirements below are useful when R attempts to compile packages from source. Make sure your system has the usual R build tools installed.

#### macOS

Install Xcode command line tools by running this in Terminal:

```bash
xcode-select --install
```

On Apple Silicon Macs, such as M1, M2, M3, or M4 machines, some R packages may also require GNU Fortran if R tries to compile packages from source. This can show up as installation errors mentioning gfortran, Fortran, or ld: library not found for -lgfortran.

If that happens, install the GNU Fortran tools recommended for R from:

```text
https://mac.r-project.org/tools/
```

After installing GNU Fortran, restart RStudio and rerun the R package installation commands below.

You can check whether gfortran is visible by running this in Terminal:

```bash
which gfortran
gfortran --version
```

If these commands do not find gfortran, R may not be able to compile packages that require Fortran code.

Note: The GNU Fortran issue is a general macOS/R package installation issue that may occur when CRAN packages need to be compiled locally.

#### Ubuntu/Debian Linux

Install system build tools by running this in Terminal:

```bash
sudo apt-get update
sudo apt-get install build-essential libcurl4-openssl-dev libssl-dev libxml2-dev
```

#### Windows

Install Rtools for your version of R:

https://cran.r-project.org/bin/windows/Rtools/

After installing Rtools, restart RStudio before installing packages.

### Step 2: Install R package dependencies

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
  "rstudioapi",
  "mcclust",
  "cluster",
  "e1071",
  "randomForest",
  "survival"
))
```

### Step 3: Confirm package loading

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
  "purrr",
  "rstudioapi",
  "mcclust",
  "cluster",
  "e1071",
  "randomForest",
  "survival"
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

The recommended way to run the app is from **RStudio** using **Run in Window**.

First, open the repository folder in RStudio. Make sure your working directory is the repository root folder.

You can check your current working directory by running:

```r
getwd()
```

If needed, set the working directory manually. For example:

```r
setwd("path/to/your/repository-folder")
```

Then open:

```text
app/app.R
```

In RStudio, click the drop-down arrow next to **Run App** and select:

```text
Run in Window
```

This is the preferred launch method for local testing. External browser mode can sometimes open multiple browser tabs if older Shiny sessions are still active.

When finished, stop the app using the red **Stop** button in RStudio before launching it again.

Alternatively, the app can be launched from the R console with:

```r
shiny::runApp("app")
```

However, for most users, **Run App > Run in Window** is recommended.

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

### 2. Optional control compatibility diagnostic

After uploading the model-fitting dataset, the app includes an optional **Control Compatibility** tab for datasets containing a concurrent control arm labeled exactly:

```text
Control
```

The diagnostic compares the observed concurrent-control survival experience with the historical Bayesian-ECD posterior predictive control distribution. It does not perform the Stage 2 Bayesian-ECD model extension used for the primary treatment analysis.

A concurrent control arm is not required to use the primary Bayesian-ECD analysis. See the user guide for details on running and interpreting the compatibility diagnostic.

### 3. View model outputs

Use the following tabs:

- **Plot Output**
- **RMST Output**
- **Subpopulation Analysis**
- **Subpopulation RMST Output**

### 4. Upload prediction data

Use the example prediction dataset in:

```text
example_data/
```

The prediction dataset should contain:

- prediction patient ID
- the same covariates used in model fitting

It does not need outcome, censoring, treatment, or cohort/source columns.

### 5. View prediction outputs

Use the following tabs:

- **Prediction Output**
- **Prediction RMST Output**


## Optional posterior probability utilities

The Shiny app focuses on survival curves, hazard-ratio curves, RMST summaries, subpopulation analyses, and prediction outputs.

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

The application is organized to run on Windows, macOS, and Linux. Platform-specific build requirements for R packages compiled from source are described above and in the user guide.
