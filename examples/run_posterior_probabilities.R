# ------------------------------------------------------------
# 0. Find repository root
# ------------------------------------------------------------
#
# This script is designed to work when opened and sourced from
# RStudio, even if getwd() is not the repository root.

get_script_path <- function() {
  # Works when the script is opened in RStudio and sourced.
  if (requireNamespace("rstudioapi", quietly = TRUE)) {
    ctx <- tryCatch(
      rstudioapi::getActiveDocumentContext(),
      error = function(e) NULL
    )
    
    if (!is.null(ctx) && nzchar(ctx$path)) {
      return(normalizePath(ctx$path, winslash = "/", mustWork = TRUE))
    }
  }
  
  # Fallback for command-line execution with Rscript.
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  
  if (length(file_arg) > 0) {
    return(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = TRUE))
  }
  
  stop(
    "Could not determine script location. ",
    "Please open examples/run_posterior_probabilities.R in RStudio and click Source, ",
    "or run it with Rscript from the command line."
  )
}

find_repo_root <- function(script_path) {
  script_dir <- dirname(script_path)
  
  candidate_dirs <- c(
    script_dir,
    normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = FALSE),
    normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  )
  
  for (candidate in unique(candidate_dirs)) {
    if (
      dir.exists(file.path(candidate, "app_core")) &&
      dir.exists(file.path(candidate, "analysis_utils")) &&
      dir.exists(file.path(candidate, "example_data")) &&
      file.exists(file.path(candidate, "README.md"))
    ) {
      return(candidate)
    }
  }
  
  stop(
    "Could not find repository root. ",
    "Please make sure this script is inside the examples/ folder of the downloaded repository."
  )
}

script_path <- get_script_path()
repo_root <- find_repo_root(script_path)

message("Script path: ", script_path)
message("Using repository root: ", repo_root)

# ============================================================
# Example: Posterior probability calculations
# ============================================================
#
# This script shows how to compute posterior probabilities from
# a Bayesian-ECD result object outside the Shiny app.
#
# Examples:
#   P(HR(t) < threshold | Data)
#   P(RMST difference > threshold | Data)
#   P(RMST ratio > threshold | Data)
#
# This script can be run from RStudio after opening:
#   examples/run_posterior_probabilities.R
# ============================================================


# ------------------------------------------------------------
# 1. Source required code
# ------------------------------------------------------------

# source_app_core.R and helpers.R use app-relative paths, so source
# them from inside the app/ directory and then restore the old directory.
old_wd <- getwd()
on.exit(setwd(old_wd), add = TRUE)

setwd(file.path(repo_root, "app"))
source("../app_core/source_app_core.R")
source("helpers.R")
setwd(old_wd)

# Source posterior probability utilities.
source(file.path(repo_root, "analysis_utils", "posterior-probability-utils.R"))


# ------------------------------------------------------------
# 2. Load example data
# ------------------------------------------------------------

model_data_path <- file.path(repo_root, "example_data", "example_model_data.csv")
prediction_data_path <- file.path(repo_root, "example_data", "example_prediction_data.csv")

if (!file.exists(model_data_path)) {
  stop("Could not find model data file: ", model_data_path)
}

if (!file.exists(prediction_data_path)) {
  stop("Could not find prediction data file: ", prediction_data_path)
}

model_data <- readr::read_csv(model_data_path, show_col_types = FALSE)
prediction_data <- readr::read_csv(prediction_data_path, show_col_types = FALSE)

message("Loaded model data with columns:")
print(names(model_data))

message("Loaded prediction data with columns:")
print(names(prediction_data))


# ------------------------------------------------------------
# 3. Define model-input specifications
# ------------------------------------------------------------
#
# IMPORTANT:
# Edit these column names if your example CSV files use different names.
#
# These should match the fields users enter in the Shiny app.

input_specs <- list(
  response = "os",
  censor_ind = "os_status",
  trt_type = "trt",
  dat_type = "cohort",
  cont_vars = "age",
  cat_vars = c("sex", "kps", "eor")
)

required_model_cols <- as.character(c(
  input_specs$response,
  input_specs$censor_ind,
  input_specs$trt_type,
  input_specs$dat_type,
  input_specs$cont_vars,
  input_specs$cat_vars
))

required_prediction_cols <- as.character(c(
  input_specs$cont_vars,
  input_specs$cat_vars
))

# Keep ID checks separate because the model function itself does not need ID.
id_col <- "ID"

if (!id_col %in% names(model_data)) {
  stop("The model dataset is missing the ID column: ", id_col)
}

if (!id_col %in% names(prediction_data)) {
  stop("The prediction dataset is missing the ID column: ", id_col)
}

missing_model_cols <- setdiff(required_model_cols, names(model_data))
missing_prediction_cols <- setdiff(required_prediction_cols, names(prediction_data))

if (length(missing_model_cols) > 0) {
  stop(
    "The model dataset is missing these required columns: ",
    paste(missing_model_cols, collapse = ", "),
    "\nEdit input_specs in this script to match your CSV file."
  )
}

if (length(missing_prediction_cols) > 0) {
  stop(
    "The prediction dataset is missing these required columns: ",
    paste(missing_prediction_cols, collapse = ", "),
    "\nEdit input_specs in this script to match your CSV file."
  )
}

# ------------------------------------------------------------
# 3b. Set categorical variable levels consistently
# ------------------------------------------------------------
#
# These levels should match the coding used by the stored baseline
# MCMC object and the Shiny app examples.
#
# The order matters because subgroup filters use numeric level codes:
#   first level  -> "0"
#   second level -> "1"
#   third level  -> "2"

cat_levels <- result.CAPPMx.base$Cat_Levels

set_example_factor_levels <- function(df, cat_levels, dataset_name = "dataset") {
  for (var_name in names(cat_levels)) {
    if (!var_name %in% names(df)) {
      stop(dataset_name, " is missing categorical variable: ", var_name)
    }
    
    observed_values <- unique(as.character(df[[var_name]]))
    allowed_values <- cat_levels[[var_name]]
    bad_values <- setdiff(observed_values, allowed_values)
    
    if (length(bad_values) > 0) {
      stop(
        dataset_name, " has unexpected value(s) in ", var_name, ": ",
        paste(bad_values, collapse = ", "),
        "\nAllowed values are: ",
        paste(allowed_values, collapse = ", ")
      )
    }
    
    df[[var_name]] <- factor(as.character(df[[var_name]]), levels = allowed_values)
  }
  
  df
}

model_data <- set_example_factor_levels(
  df = model_data,
  cat_levels = cat_levels,
  dataset_name = "model_data"
)

prediction_data <- set_example_factor_levels(
  df = prediction_data,
  cat_levels = cat_levels,
  dataset_name = "prediction_data"
)

str(prediction_data[, c("sex", "kps", "eor")])


# ------------------------------------------------------------
# 4. Check stored baseline MCMC object
# ------------------------------------------------------------
#
# helpers.R creates the object:
#   result.CAPPMx.base

if (!exists("result.CAPPMx.base")) {
  stop("result.CAPPMx.base was not created. Check app/helpers.R and app/data/result.CAPPMx.base/.")
}

message("Stored MCMC object loaded successfully.")

message("Available treatment levels in result.CAPPMx.base:")
print(result.CAPPMx.base$Treatment_Levels)

control_treatment <- result.CAPPMx.base$Treatment_Levels

# ------------------------------------------------------------
# 5. Run Bayesian-ECD model with prediction data
# ------------------------------------------------------------
#
# This mirrors the Shiny app prediction workflow.
#
# The key argument is:
#   input_df_pred = prediction_data
#
# That creates Predicted_Allocation_variables in pred_result.

message("Running Bayesian-ECD prediction workflow. This may take a moment.")

pred_result <- cappmx_extend_approx_fit(
  result_CAPPMx = result.CAPPMx.base,
  input_df = model_data,
  input_specs = input_specs,
  ref_trt = control_treatment,
  input_df_pred = prediction_data,
  del_range_response_1 = c(0.005, 0.02) * 8,
  del_range_response_2 = c(0.005, 0.02) * 9,
  del_range_alp1 = c(0.1, 0.3) * 2.8
)
 
if (!("Predicted_Allocation_variables" %in% names(pred_result))) {
  stop("pred_result does not contain Predicted_Allocation_variables.")
}

message("Prediction result created successfully.")


# ------------------------------------------------------------
# 6. Posterior probability for hazard ratio
# ------------------------------------------------------------
#
# Example:
#   P(HR at 2 years < 0.8 | Data)
#
# use_pred = TRUE means the function uses Predicted_Allocation_variables.

print(pred_result$Treatment_Levels)

comparison_treatment <- "Drug A"

hr_prob <- posterior_prob_hr_lt(
  result = pred_result,
  input_df = prediction_data,
  cntrl = control_treatment,
  trt = comparison_treatment,
  timepoint_days = 730.5,
  threshold = c(0.6, 0.8, 1.0),
  burnin = 200,
  subgroup_list = NULL,
  use_pred = TRUE
)

message("Posterior probability results for HR:")
print(hr_prob$posterior_probability)

message("Posterior HR summary:")
print(hr_prob$hr_summary)


# ------------------------------------------------------------
# 7. Posterior probability for RMST
# ------------------------------------------------------------
#
# Examples:
#   P(RMST difference > 60 days | Data)
#   P(RMST ratio > 1.2 | Data)

rmst_prob <- posterior_prob_rmst(
  result = pred_result,
  input_df = prediction_data,
  cntrl = control_treatment,
  trt = comparison_treatment,
  time_horizons = c(365.25, 730.5),
  ratio_threshold = c(1.0, 1.2, 1.4),
  diff_threshold = c(0, 60, 120),
  burnin = 200,
  subgroup_list = NULL,
  use_pred = TRUE,
  n_time_grid = 200
)

message("Posterior probability results for RMST difference at 730.5 days:")
print(rmst_prob$results_by_horizon[["730.5_days"]]$posterior_probability_diff)

message("Posterior probability results for RMST ratio at 730.5 days:")
print(rmst_prob$results_by_horizon[["730.5_days"]]$posterior_probability_ratio)

message("Posterior RMST summary table at 730.5 days:")
print(rmst_prob$results_by_horizon[["730.5_days"]]$summary_table)


# ------------------------------------------------------------
# 8. Optional subgroup example
# ------------------------------------------------------------
#
# This section is optional.
#
# subgroup_list uses the same convention as the utility functions:
#
#   categorical variables:
#     use character-coded level indices such as "0", "01", or "012"
#
#   continuous variables:
#     use a numeric interval c(lower, upper)
#
# Example below:
#   sex = "0"
#   eor = "01"
#   kps = NA means no KPS filtering
#   age = c(50, 70)
#
# If the example dataset does not use these coding conventions,
# edit or skip this section.

subgroup_list <- list(
  sex = "0",
  eor = "01",
  kps = NA,
  age = c(50, 70)
)

message("Running optional subpopulation posterior probability example.")

hr_prob_subgroup <- posterior_prob_hr_lt(
  result = pred_result,
  input_df = prediction_data,
  cntrl = control_treatment,
  trt = comparison_treatment,
  timepoint_days = 730.5,
  threshold = c(0.6, 0.8, 1.0),
  burnin = 200,
  subgroup_list = subgroup_list,
  use_pred = TRUE
)

message("Posterior probability results for subpopulation HR:")
print(hr_prob_subgroup$posterior_probability)

message("Posterior subpopulation HR summary:")
print(hr_prob_subgroup$hr_summary)


rmst_prob_subgroup <- posterior_prob_rmst(
  result = pred_result,
  input_df = prediction_data,
  cntrl = control_treatment,
  trt = comparison_treatment,
  time_horizons = c(365.25, 730.5),
  ratio_threshold = c(1.0, 1.2, 1.4),
  diff_threshold = c(0, 60, 120),
  burnin = 200,
  subgroup_list = subgroup_list,
  use_pred = TRUE,
  n_time_grid = 200
)

message("Posterior probability results for subpopulation RMST difference at 730.5 days:")
print(rmst_prob_subgroup$results_by_horizon[["730.5_days"]]$posterior_probability_diff)

message("Posterior probability results for subpopulation RMST ratio at 730.5 days:")
print(rmst_prob_subgroup$results_by_horizon[["730.5_days"]]$posterior_probability_ratio)

message("Posterior subpopulation RMST summary table at 730.5 days:")
print(rmst_prob_subgroup$results_by_horizon[["730.5_days"]]$summary_table)


# ------------------------------------------------------------
# 9. Save posterior probability outputs
# ------------------------------------------------------------

output_dir <- file.path(repo_root, "outputs")
dir.create(output_dir, showWarnings = FALSE)

saveRDS(
  hr_prob,
  file.path(output_dir, "posterior_probability_hr_prediction.rds")
)

saveRDS(
  rmst_prob,
  file.path(output_dir, "posterior_probability_rmst_prediction.rds")
)

write.csv(
  hr_prob$posterior_probability,
  file.path(output_dir, "posterior_probability_hr_prediction.csv"),
  row.names = TRUE
)

write.csv(
  rmst_prob$results_by_horizon[["730.5_days"]]$posterior_probability_diff,
  file.path(output_dir, "posterior_probability_rmst_difference_prediction.csv"),
  row.names = TRUE
)

write.csv(
  rmst_prob$results_by_horizon[["730.5_days"]]$posterior_probability_ratio,
  file.path(output_dir, "posterior_probability_rmst_ratio_prediction.csv"),
  row.names = TRUE
)

message("Saved posterior probability outputs to: ", output_dir)


# ------------------------------------------------------------
# 11. Done
# ------------------------------------------------------------

message("Posterior probability example completed successfully.")



