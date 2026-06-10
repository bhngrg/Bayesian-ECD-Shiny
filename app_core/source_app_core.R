# Source all core functions required by the Bayesian-ECD Shiny app.
#
# This file is designed to work whether it is sourced from:
#   1. app/app.R
#   2. the repository root
#   3. app_core/source_app_core.R directly

required_packages <- c(
  "Rcpp",
  "RcppArmadillo",
  "RcppDist",
  "RcppGSL",
  "matrixStats",
  "mcclust"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "The following required package(s) are missing or not visible to this R session:\n",
    paste(missing_packages, collapse = ", "),
    "\n\nPlease install them with:\n",
    "install.packages(c(",
    paste(sprintf('"%s"', missing_packages), collapse = ", "),
    "))\n\n",
    "If RcppGSL fails to install, check that the GNU Scientific Library (GSL) ",
    "and system build tools are installed for your operating system.",
    call. = FALSE
  )
}

get_current_source_file <- function() {
  frames <- sys.frames()
  
  source_files <- vapply(
    frames,
    function(frame) {
      if (!is.null(frame$ofile)) {
        return(frame$ofile)
      }
      NA_character_
    },
    character(1)
  )
  
  source_files <- source_files[!is.na(source_files)]
  
  if (length(source_files) == 0) {
    return(NA_character_)
  }
  
  normalizePath(source_files[length(source_files)], winslash = "/", mustWork = TRUE)
}

find_repo_root <- function() {
  this_file <- get_current_source_file()
  
  candidate_dirs <- c(
    getwd(),
    file.path(getwd(), ".."),
    if (!is.na(this_file)) dirname(dirname(this_file)) else NA_character_,
    if (!is.na(this_file)) dirname(this_file) else NA_character_
  )
  
  candidate_dirs <- unique(normalizePath(
    candidate_dirs[!is.na(candidate_dirs)],
    winslash = "/",
    mustWork = FALSE
  ))
  
  for (candidate in candidate_dirs) {
    if (
      file.exists(file.path(candidate, "app_core", "Robust_stage_1and2.cpp")) &&
      file.exists(file.path(candidate, "app_core", "survival-utils.R")) &&
      file.exists(file.path(candidate, "app", "app.R"))
    ) {
      return(candidate)
    }
  }
  
  stop(
    "Could not find the repository root. Please run the app from the downloaded ",
    "Bayesian-ECD-Shiny repository folder.",
    call. = FALSE
  )
}

repo_root <- find_repo_root()
app_core_dir <- file.path(repo_root, "app_core")

message("Using repository root: ", repo_root)

Rcpp::sourceCpp(file.path(app_core_dir, "Robust_stage_1and2.cpp"))
source(file.path(app_core_dir, "survival-utils.R"))
source(file.path(app_core_dir, "cappmx-extend-approx-fit.R"))
source(file.path(app_core_dir, "plotting.R"))
source(file.path(app_core_dir, "rmst-utils.R"))
source(file.path(app_core_dir, "subgroup-analysis.R"))