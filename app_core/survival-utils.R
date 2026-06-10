# Survival and density helper functions for lognormal mixture output.
# These functions depend on the Rcpp function log_sum_exp().

# Survival and density helper functions.

density.fn.lognorm <- function(time, mu, sig, weights) {
  log_terms <- log(weights) + stats::dnorm(
    x = time,
    mean = mu,
    sd = sig,
    log = TRUE
  )
  
  matrixStats::logSumExp(log_terms)
}


survival.fn.lognorm <- function(time, mu, sig, weights) {
  log_terms <- log(weights) + stats::pnorm(
    q = time,
    mean = mu,
    sd = sig,
    lower.tail = FALSE,
    log.p = TRUE
  )
  
  matrixStats::logSumExp(log_terms)
}