# Cross-platform parallel-computing helpers for the Bayesian-ECD Shiny app.
#
# Application-level parallelism is provided explicitly through PSOCK workers.
# Native numerical libraries are restricted to one thread per R process so
# that each worker does not independently create additional BLAS/OpenMP
# threads.

bayesian_ecd_single_thread_env <- c(
  OMP_NUM_THREADS = "1",
  OMP_THREAD_LIMIT = "1",
  OMP_MAX_ACTIVE_LEVELS = "1",
  OPENBLAS_NUM_THREADS = "1",
  GOTO_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  MKL_DYNAMIC = "FALSE",
  VECLIB_MAXIMUM_THREADS = "1",
  BLIS_NUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1"
)


set_bayesian_ecd_single_thread <- function() {
  do.call(
    Sys.setenv,
    as.list(bayesian_ecd_single_thread_env)
  )

  if (
    requireNamespace("RcppArmadillo", quietly = TRUE) &&
    exists(
      "armadillo_set_number_of_omp_threads",
      envir = asNamespace("RcppArmadillo"),
      inherits = FALSE
    )
  ) {
    try(
      RcppArmadillo::armadillo_set_number_of_omp_threads(1L),
      silent = TRUE
    )
  }

  invisible(NULL)
}


get_bayesian_ecd_worker_count <- function(n_tasks = Inf) {
  logical_cores <- suppressWarnings(
    parallel::detectCores(logical = TRUE)
  )

  physical_cores <- suppressWarnings(
    parallel::detectCores(logical = FALSE)
  )

  if (
    length(logical_cores) != 1L ||
    is.na(logical_cores) ||
    !is.finite(logical_cores) ||
    logical_cores < 1
  ) {
    logical_cores <- 1L
  }

  logical_cores <- as.integer(logical_cores)

  physical_is_usable <-
    length(physical_cores) == 1L &&
    !is.na(physical_cores) &&
    is.finite(physical_cores) &&
    physical_cores >= 1 &&
    physical_cores < logical_cores

  if (physical_is_usable) {
    available_cores <- as.integer(physical_cores)
  } else if (logical_cores > 1L) {
    # Some platforms do not reliably distinguish physical from logical
    # processors. Use a conservative SMT/Hyper-Threading fallback instead
    # of creating one worker for nearly every reported logical processor.
    available_cores <- max(
      as.integer(floor(logical_cores / 2)),
      1L
    )
  } else {
    available_cores <- 1L
  }

  # Leave one estimated execution core available for the Shiny/R process.
  nworkers <- max(available_cores - 1L, 1L)

  if (!is.infinite(n_tasks)) {
    if (
      length(n_tasks) != 1L ||
      is.na(n_tasks) ||
      !is.finite(n_tasks) ||
      n_tasks < 1
    ) {
      stop(
        "n_tasks must be a positive finite number or Inf.",
        call. = FALSE
      )
    }

    nworkers <- min(
      nworkers,
      as.integer(ceiling(n_tasks))
    )
  }

  as.integer(max(nworkers, 1L))
}


make_bayesian_ecd_cluster <- function(n_tasks = Inf) {
  nworkers <- get_bayesian_ecd_worker_count(
    n_tasks = n_tasks
  )

  cl <- parallel::makeCluster(nworkers)

  parallel::clusterCall(
    cl,
    function(thread_env) {
      do.call(
        Sys.setenv,
        as.list(thread_env)
      )

      if (
        requireNamespace("RcppArmadillo", quietly = TRUE) &&
        exists(
          "armadillo_set_number_of_omp_threads",
          envir = asNamespace("RcppArmadillo"),
          inherits = FALSE
        )
      ) {
        try(
          RcppArmadillo::armadillo_set_number_of_omp_threads(1L),
          silent = TRUE
        )
      }

      invisible(NULL)
    },
    thread_env = bayesian_ecd_single_thread_env
  )

  cl
}
