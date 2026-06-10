# Source all app-core code in the correct order.
# This script assumes it is called from cappmx_package_prep/app/app.R.

Rcpp::sourceCpp("../app_core/Robust_stage_1and2.cpp")

source("../app_core/survival-utils.R")
source("../app_core/cappmx-extend-approx-fit.R")
source("../app_core/plotting.R")
source("../app_core/rmst-utils.R")
source("../app_core/subgroup-analysis.R")

core_functions <- c(
  "log_sum_exp",
  "surv_fn_lognorm",
  "post_t_dens",
  "stage2_cov_loglik_precompute",
  "stage2_init_alloc_cov_only",
  "background_MCMC_storage",
  "common_atoms_cat_lognormal_shared_approx",
  "extract_MCMC_storage",
  "cappmx_extend_approx_fit",
  "density.fn.lognorm",
  "survival.fn.lognorm",
  "plots_lognormal",
  "rmst_lognormal",
  "RMST_result",
  "subgroup_data"
)

missing_functions <- core_functions[!vapply(core_functions, exists, logical(1))]

if (length(missing_functions) > 0) {
  stop(
    "Missing app-core functions: ",
    paste(missing_functions, collapse = ", ")
  )
}

message("All app-core functions sourced successfully.")