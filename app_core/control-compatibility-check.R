# Posterior predictive concurrent-control compatibility check.
#
# This optional diagnostic compares the uploaded concurrent-control arm
# with Historical-Control posterior predictions standardized to the
# concurrent-Control covariate population.
#
# Primary compatibility statistic:
#   1. Fit Stage 2 with concurrent-Control covariates supplied for prediction.
#   2. Retain 1000 posterior draws after 200 stored-draw burn-in.
#   3. Generate one uncensored Historical-Control potential outcome per
#      concurrent-Control covariate profile and retained posterior draw.
#   4. Compute a log-rank chi-square statistic for every posterior draw.
#   5. Use the MEAN posterior log-rank statistic as the observed statistic.
#   6. Compare it with a paired nonparametric bootstrap null distribution
#      generated from the uploaded concurrent Controls.
#   7. Report p = r / B, where r is the number of bootstrap statistics at
#      least as large as the mean posterior statistic.
#
# RNG streams are deterministic so unchanged data/settings reproduce the
# same compatibility result even when bootstrap computation is parallel.

CONTROL_COMPAT_STAGE2_N_STORED <- 1200L
CONTROL_COMPAT_STAGE2_BURNIN <- 200L
CONTROL_COMPAT_N_POSTERIOR_DRAWS <-
  CONTROL_COMPAT_STAGE2_N_STORED - CONTROL_COMPAT_STAGE2_BURNIN

CONTROL_COMPAT_BOOTSTRAP_REPLICATES <- 500L

CONTROL_COMPAT_STAGE2_SEED <- 202609031L
CONTROL_COMPAT_POTENTIAL_OUTCOME_SEED <- 702609031L
CONTROL_COMPAT_BOOTSTRAP_TRT1_SEED_BASE <- 502609030L
CONTROL_COMPAT_BOOTSTRAP_TRT2_SEED_BASE <- 552609030L


compute_uploaded_control_km_curve <- function(
    control_df,
    time_col,
    status_col,
    conf_int = 0.95
) {
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Package 'survival' is required for the control compatibility check.")
  }

  if (!time_col %in% names(control_df)) {
    stop("time_col was not found in the uploaded control data: ", time_col)
  }

  if (!status_col %in% names(control_df)) {
    stop("status_col was not found in the uploaded control data: ", status_col)
  }

  time_values <- suppressWarnings(as.numeric(control_df[[time_col]]))
  status_values <- as.logical(control_df[[status_col]])

  keep <- is.finite(time_values) &
    time_values > 0 &
    !is.na(status_values)

  if (sum(keep) < 5L) {
    stop(
      "Fewer than 5 uploaded concurrent-control patients have usable survival ",
      "time and censoring information."
    )
  }

  tmp_df <- data.frame(
    time = time_values[keep],
    event = status_values[keep]
  )

  fit <- survival::survfit(
    survival::Surv(time = time, event = event) ~ 1,
    data = tmp_df,
    conf.int = conf_int
  )

  km_summary <- summary(fit)

  km_curve <- data.frame(
    time = c(0, km_summary$time),
    survival = c(1, km_summary$surv),
    lower = c(1, km_summary$lower),
    upper = c(1, km_summary$upper),
    curve = "Uploaded concurrent control KM"
  )

  fit_table <- summary(fit)$table

  observed_median <- suppressWarnings({
    if (is.null(dim(fit_table))) {
      as.numeric(unname(fit_table["median"]))
    } else if ("median" %in% colnames(fit_table)) {
      as.numeric(fit_table[1, "median"])
    } else {
      NA_real_
    }
  })

  if (length(observed_median) == 0L) {
    observed_median <- NA_real_
  }

  list(
    fit = fit,
    curve = km_curve,
    observed_median = observed_median,
    n_evaluable = nrow(tmp_df),
    evaluable_data = tmp_df,
    keep = keep
  )
}


make_control_compatibility_seed <- function(
    base_seed,
    index
) {
  seed <- as.double(base_seed) + as.double(index)

  if (
    !is.finite(seed) ||
    seed <= 0 ||
    seed >= .Machine$integer.max
  ) {
    stop("Invalid deterministic compatibility seed.")
  }

  as.integer(seed)
}


fit_control_compatibility_stage2 <- function(
    result_CAPPMx,
    input_df,
    input_specs,
    prediction_data,
    control_label = "Control",
    seed = CONTROL_COMPAT_STAGE2_SEED
) {
  set.seed(seed)

  invisible(utils::capture.output(
    result <- cappmx_extend_approx_fit(
      result_CAPPMx = result_CAPPMx,
      input_df = as.data.frame(input_df),
      input_specs = input_specs,
      ref_trt = control_label,
      input_df_pred = as.data.frame(prediction_data),
      del_range_response_1 = c(0.005, 0.02) * 8,
      del_range_response_2 = c(0.005, 0.02) * 9,
      del_range_alp1 = c(0.1, 0.3) * 2.3,
      freeze_control = TRUE,
      eb_beta_unique = TRUE,
      a0 = 10.1,
      cc = 2,
      beta_temper_tau = 10,
      beta_cap_q = c(0.10, 0.90),
      beta_var_floor_mult = 1.0,
      min_failures_for_EB = 3,
      stage2_n_stored = CONTROL_COMPAT_STAGE2_N_STORED
    ),
    type = "output"
  ))

  required_fields <- c(
    "Treatment_Levels",
    "Treatment_Indices",
    "Predicted_Allocation_variables",
    "Lognormal_Mu_Cube",
    "Lognormal_Sig_Cube"
  )

  missing_fields <- setdiff(required_fields, names(result))

  if (length(missing_fields) > 0L) {
    stop(
      "Compatibility Stage-2 result is missing field(s): ",
      paste(missing_fields, collapse = ", ")
    )
  }

  if (
    nrow(result$Predicted_Allocation_variables) !=
      CONTROL_COMPAT_STAGE2_N_STORED ||
    ncol(result$Predicted_Allocation_variables) != nrow(prediction_data)
  ) {
    stop("Unexpected Predicted_Allocation_variables dimensions.")
  }

  result
}


extract_control_compatibility_posterior_state <- function(
    result,
    prediction_data,
    control_label = "Control",
    burnin = CONTROL_COMPAT_STAGE2_BURNIN
) {
  alloc <- result$Predicted_Allocation_variables

  keep <- seq.int(
    burnin + 1L,
    nrow(alloc)
  )

  if (length(keep) != CONTROL_COMPAT_N_POSTERIOR_DRAWS) {
    stop("Unexpected retained posterior-draw count.")
  }

  if (ncol(alloc) != nrow(prediction_data)) {
    stop("Prediction population does not match predicted allocations.")
  }

  allocation_0based <- matrix(
    as.integer(alloc[keep, , drop = FALSE]),
    nrow = length(keep),
    ncol = ncol(alloc)
  )

  if (
    anyNA(allocation_0based) ||
    any(allocation_0based < 0L)
  ) {
    stop("Invalid retained 0-based predicted allocations.")
  }

  treatment_ind <- unique(
    result$Treatment_Indices[
      result$Treatment_Levels == control_label
    ]
  )

  if (length(treatment_ind) != 1L) {
    stop(
      "Could not identify one treatment index for ",
      control_label,
      "."
    )
  }

  n_draws <- length(keep)
  n_pred <- nrow(prediction_data)

  mu <- matrix(
    NA_real_,
    nrow = n_draws,
    ncol = n_pred
  )

  sig2 <- matrix(
    NA_real_,
    nrow = n_draws,
    ncol = n_pred
  )

  for (draw_position in seq_len(n_draws)) {
    m <- keep[[draw_position]]

    cluster_1b <-
      allocation_0based[draw_position, ] + 1L

    mu_draw <- as.numeric(
      result$Lognormal_Mu_Cube[
        m,
        ,
        treatment_ind
      ]
    )

    sig2_draw <- as.numeric(
      result$Lognormal_Sig_Cube[
        m,
        ,
        treatment_ind
      ]
    )

    if (
      any(cluster_1b < 1L) ||
      any(cluster_1b > length(mu_draw))
    ) {
      stop(
        "Invalid predicted allocation while extracting Historical-Control ",
        "posterior state."
      )
    }

    mu[draw_position, ] <- mu_draw[cluster_1b]
    sig2[draw_position, ] <- sig2_draw[cluster_1b]
  }

  if (
    any(!is.finite(mu)) ||
    any(!is.finite(sig2)) ||
    any(sig2 <= 0)
  ) {
    stop("Invalid retained Historical-Control posterior parameters.")
  }

  list(
    retained_stage2_rows = keep,
    allocation_0based = allocation_0based,
    treatment_label = control_label,
    treatment_index = as.integer(treatment_ind),
    mu = mu,
    sig2 = sig2
  )
}


generate_control_compatibility_potential_outcomes <- function(
    posterior_state,
    seed = CONTROL_COMPAT_POTENTIAL_OUTCOME_SEED
) {
  mu <- posterior_state$mu
  sig2 <- posterior_state$sig2

  if (!identical(dim(mu), dim(sig2))) {
    stop("Historical-Control mu and sig2 dimensions do not agree.")
  }

  n_draws <- nrow(mu)
  n_pred <- ncol(mu)

  set.seed(seed)

  z <- matrix(
    stats::rnorm(n_draws * n_pred),
    nrow = n_draws,
    ncol = n_pred
  )

  potential_event_time <- exp(
    mu + sqrt(sig2) * z
  )

  if (
    any(!is.finite(potential_event_time)) ||
    any(potential_event_time <= 0)
  ) {
    stop("Generated Historical-Control potential event times are invalid.")
  }

  list(
    retained_stage2_rows = posterior_state$retained_stage2_rows,
    allocation = posterior_state$allocation_0based,
    mu = mu,
    sig2 = sig2,
    z = z,
    potential_event_time = potential_event_time,
    potential_event_status = matrix(
      1L,
      nrow = n_draws,
      ncol = n_pred
    )
  )
}


compute_logrank_statistic <- function(
    time_1,
    status_1,
    time_2,
    status_2
) {
  time_1 <- as.numeric(time_1)
  status_1 <- as.integer(status_1)
  time_2 <- as.numeric(time_2)
  status_2 <- as.integer(status_2)

  if (
    length(time_1) != length(status_1) ||
    length(time_2) != length(status_2)
  ) {
    stop("Log-rank time/status vectors have incompatible lengths.")
  }

  if (
    length(time_1) < 1L ||
    length(time_2) < 1L
  ) {
    stop("Log-rank comparison requires non-empty groups.")
  }

  if (
    any(!is.finite(c(time_1, time_2))) ||
    any(c(time_1, time_2) <= 0)
  ) {
    stop("Log-rank comparison contains invalid survival times.")
  }

  if (
    any(!status_1 %in% c(0L, 1L)) ||
    any(!status_2 %in% c(0L, 1L))
  ) {
    stop("Log-rank status must be coded 0/1.")
  }

  logrank_data <- data.frame(
    time = c(time_1, time_2),
    status = c(status_1, status_2),
    group = factor(
      c(
        rep("group_1", length(time_1)),
        rep("group_2", length(time_2))
      ),
      levels = c("group_1", "group_2")
    )
  )

  fit <- survival::survdiff(
    survival::Surv(time, status) ~ group,
    data = logrank_data,
    rho = 0
  )

  statistic <- as.numeric(fit$chisq)

  if (
    length(statistic) != 1L ||
    !is.finite(statistic) ||
    statistic < 0
  ) {
    stop("Invalid log-rank chi-square statistic.")
  }

  statistic
}


compute_control_compatibility_logrank_draws <- function(
    observed_time,
    observed_status,
    historical_potential_outcomes
) {
  potential_time <-
    historical_potential_outcomes$potential_event_time

  potential_status <-
    historical_potential_outcomes$potential_event_status

  if (!identical(dim(potential_time), dim(potential_status))) {
    stop("Historical-Control potential-outcome dimensions do not agree.")
  }

  if (ncol(potential_time) != length(observed_time)) {
    stop(
      "Observed concurrent Controls do not match the Historical-Control ",
      "prediction population."
    )
  }

  statistics <- vapply(
    seq_len(nrow(potential_time)),
    function(draw_index) {
      compute_logrank_statistic(
        time_1 = observed_time,
        status_1 = observed_status,
        time_2 = potential_time[draw_index, ],
        status_2 = potential_status[draw_index, ]
      )
    },
    numeric(1)
  )

  list(
    draws = statistics,
    median = stats::median(statistics),
    mean = mean(statistics)
  )
}


compute_covariate_specific_historical_control_curve <- function(
    posterior_state,
    time_grid,
    conf_int = 0.95
) {
  mu <- posterior_state$mu
  sig2 <- posterior_state$sig2

  if (!identical(dim(mu), dim(sig2))) {
    stop("Historical-Control posterior parameter dimensions do not agree.")
  }

  time_grid <- as.numeric(time_grid)

  if (
    any(!is.finite(time_grid)) ||
    any(time_grid <= 0)
  ) {
    stop("Historical-Control curve requires positive finite time points.")
  }

  survival_draw_matrix <- vapply(
    time_grid,
    function(time_value) {
      rowMeans(
        stats::pnorm(
          (log(time_value) - mu) / sqrt(sig2),
          lower.tail = FALSE
        )
      )
    },
    numeric(nrow(mu))
  )

  # vapply returns posterior draws x time; transpose to time x draws.
  survival_draw_matrix <- t(survival_draw_matrix)

  low_quant <- (1 - conf_int) / 2
  high_quant <- 1 - low_quant

  data.frame(
    time = time_grid,
    survival = apply(
      survival_draw_matrix,
      1,
      stats::median,
      na.rm = TRUE
    ),
    lower = apply(
      survival_draw_matrix,
      1,
      stats::quantile,
      probs = low_quant,
      names = FALSE,
      na.rm = TRUE
    ),
    upper = apply(
      survival_draw_matrix,
      1,
      stats::quantile,
      probs = high_quant,
      names = FALSE,
      na.rm = TRUE
    ),
    curve = "Historical posterior predictive control"
  )
}


run_one_control_compatibility_bootstrap <- function(
    bootstrap_index,
    observed_time,
    observed_status,
    trt1_seed_base = CONTROL_COMPAT_BOOTSTRAP_TRT1_SEED_BASE,
    trt2_seed_base = CONTROL_COMPAT_BOOTSTRAP_TRT2_SEED_BASE
) {
  n_control <- length(observed_time)

  trt1_seed <- make_control_compatibility_seed(
    trt1_seed_base,
    bootstrap_index
  )

  trt2_seed <- make_control_compatibility_seed(
    trt2_seed_base,
    bootstrap_index
  )

  set.seed(trt1_seed)

  trt1_positions <- sample.int(
    n = n_control,
    size = n_control,
    replace = TRUE
  )

  set.seed(trt2_seed)

  trt2_positions <- sample.int(
    n = n_control,
    size = n_control,
    replace = TRUE
  )

  statistic <- compute_logrank_statistic(
    time_1 = observed_time[trt1_positions],
    status_1 = observed_status[trt1_positions],
    time_2 = observed_time[trt2_positions],
    status_2 = observed_status[trt2_positions]
  )

  list(
    bootstrap_index = bootstrap_index,
    logrank_statistic = statistic,
    trt1_seed = trt1_seed,
    trt2_seed = trt2_seed,
    trt1_positions = trt1_positions,
    trt2_positions = trt2_positions
  )
}


run_control_compatibility_bootstrap <- function(
    observed_time,
    observed_status,
    n_bootstrap = CONTROL_COMPAT_BOOTSTRAP_REPLICATES
) {
  n_bootstrap <- as.integer(n_bootstrap)

  if (
    length(n_bootstrap) != 1L ||
    is.na(n_bootstrap) ||
    n_bootstrap < 1L
  ) {
    stop("n_bootstrap must be a positive integer.")
  }

  ncores <- max(parallel::detectCores() - 1L, 1L)
  ncores <- min(ncores, n_bootstrap)

  cl <- parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl)

  on.exit({
    try(parallel::stopCluster(cl), silent = TRUE)
    try(foreach::registerDoSEQ(), silent = TRUE)
  }, add = TRUE)

  bootstrap_objects <- foreach::`%dopar%`(
    foreach::foreach(
      bootstrap_index = seq_len(n_bootstrap),
      .packages = "survival",
      .export = c(
        "make_control_compatibility_seed",
        "compute_logrank_statistic",
        "run_one_control_compatibility_bootstrap",
        "CONTROL_COMPAT_BOOTSTRAP_TRT1_SEED_BASE",
        "CONTROL_COMPAT_BOOTSTRAP_TRT2_SEED_BASE"
      )
    ),
    {
      run_one_control_compatibility_bootstrap(
        bootstrap_index = bootstrap_index,
        observed_time = observed_time,
        observed_status = observed_status
      )
    }
  )

  statistics <- vapply(
    bootstrap_objects,
    function(x) x$logrank_statistic,
    numeric(1)
  )

  if (any(!is.finite(statistics))) {
    stop("Non-finite bootstrap log-rank statistic.")
  }

  list(
    objects = bootstrap_objects,
    statistics = statistics,
    q95 = as.numeric(
      stats::quantile(
        statistics,
        probs = 0.95,
        names = FALSE,
        type = 7
      )
    ),
    n_bootstrap = n_bootstrap,
    ncores = ncores
  )
}


run_control_compatibility_check <- function(
    result_CAPPMx,
    uploaded_data,
    input_specs,
    control_label = "Control",
    time_grid = seq(1, 1251, length.out = 201),
    conf_int = 0.95,
    n_bootstrap = CONTROL_COMPAT_BOOTSTRAP_REPLICATES,
    stage2_seed = CONTROL_COMPAT_STAGE2_SEED,
    potential_outcome_seed = CONTROL_COMPAT_POTENTIAL_OUTCOME_SEED
) {
  required_cols <- c(
    input_specs$response,
    input_specs$censor_ind,
    input_specs$trt_type
  )

  missing_cols <- setdiff(required_cols, names(uploaded_data))

  if (length(missing_cols) > 0L) {
    stop(
      "Uploaded data are missing required column(s): ",
      paste(missing_cols, collapse = ", ")
    )
  }

  control_rows <- which(
    as.character(uploaded_data[[input_specs$trt_type]]) ==
      control_label
  )

  if (length(control_rows) == 0L) {
    stop(
      "No concurrent-control patients were found for control label '",
      control_label,
      "'."
    )
  }

  control_df <- uploaded_data[
    control_rows,
    ,
    drop = FALSE
  ]

  km_result <- compute_uploaded_control_km_curve(
    control_df = control_df,
    time_col = input_specs$response,
    status_col = input_specs$censor_ind,
    conf_int = conf_int
  )

  # Prediction and log-rank calculations must use the same evaluable
  # concurrent-Control patients represented in the KM curve.
  control_df <- control_df[
    km_result$keep,
    ,
    drop = FALSE
  ]

  prediction_variables <- unique(c(
    input_specs$cat_vars,
    input_specs$cont_vars
  ))

  prediction_variables <- prediction_variables[
    !is.na(prediction_variables) &
      nzchar(prediction_variables)
  ]

  missing_prediction_variables <- setdiff(
    prediction_variables,
    names(control_df)
  )

  if (length(missing_prediction_variables) > 0L) {
    stop(
      "Concurrent-control prediction data are missing covariate(s): ",
      paste(missing_prediction_variables, collapse = ", ")
    )
  }

  prediction_data <- as.data.frame(
    control_df[
      ,
      prediction_variables,
      drop = FALSE
    ]
  )

  compatibility_stage2_result <-
    fit_control_compatibility_stage2(
      result_CAPPMx = result_CAPPMx,
      input_df = uploaded_data,
      input_specs = input_specs,
      prediction_data = prediction_data,
      control_label = control_label,
      seed = stage2_seed
    )

  posterior_state <-
    extract_control_compatibility_posterior_state(
      result = compatibility_stage2_result,
      prediction_data = prediction_data,
      control_label = control_label,
      burnin = CONTROL_COMPAT_STAGE2_BURNIN
    )

  historical_potential_outcomes <-
    generate_control_compatibility_potential_outcomes(
      posterior_state = posterior_state,
      seed = potential_outcome_seed
    )

  observed_time <- as.numeric(
    control_df[[input_specs$response]]
  )

  observed_status <- as.integer(
    as.logical(control_df[[input_specs$censor_ind]])
  )

  posterior_logrank <-
    compute_control_compatibility_logrank_draws(
      observed_time = observed_time,
      observed_status = observed_status,
      historical_potential_outcomes =
        historical_potential_outcomes
    )

  bootstrap_result <-
    run_control_compatibility_bootstrap(
      observed_time = observed_time,
      observed_status = observed_status,
      n_bootstrap = n_bootstrap
    )

  exceed_mean <- sum(
    bootstrap_result$statistics >=
      posterior_logrank$mean
  )

  p_boot_mean <-
    exceed_mean /
    bootstrap_result$n_bootstrap

  compatible <- p_boot_mean >= 0.05

  posterior_curve <-
    compute_covariate_specific_historical_control_curve(
      posterior_state = posterior_state,
      time_grid = time_grid,
      conf_int = conf_int
    )

  summary_df <- data.frame(
    control_label = control_label,
    n_control = nrow(control_df),
    n_evaluable_control = km_result$n_evaluable,
    observed_km_median = km_result$observed_median,
    posterior_logrank_mean = posterior_logrank$mean,
    posterior_logrank_median = posterior_logrank$median,
    bootstrap_q95 = bootstrap_result$q95,
    bootstrap_exceedances_mean = exceed_mean,
    bootstrap_replicates = bootstrap_result$n_bootstrap,
    bootstrap_p_value_mean = p_boot_mean,
    compatible = compatible
  )

  plot_data <- rbind(
    km_result$curve,
    posterior_curve
  )

  list(
    summary = summary_df,
    km_curve = km_result$curve,
    posterior_curve = posterior_curve,
    plot_data = plot_data,
    prediction_data = prediction_data,
    compatibility_stage2_result = compatibility_stage2_result,
    posterior_state = posterior_state,
    historical_potential_outcomes = historical_potential_outcomes,
    posterior_logrank = posterior_logrank,
    bootstrap = bootstrap_result
  )
}


plot_control_compatibility_check <- function(
    compatibility_result,
    x_axis_min = 0,
    x_axis_max = NULL
) {
  plot_data <- compatibility_result$plot_data
  summary_df <- compatibility_result$summary

  posterior_data <- plot_data[
    plot_data$curve == "Historical posterior predictive control",
    ,
    drop = FALSE
  ]

  km_data <- plot_data[
    plot_data$curve == "Uploaded concurrent control KM",
    ,
    drop = FALSE
  ]

  observed_median <- suppressWarnings(
    as.numeric(summary_df$observed_km_median[1])
  )

  if (is.null(x_axis_max)) {
    x_axis_max <- max(plot_data$time, na.rm = TRUE)
  }

  p <- ggplot2::ggplot() +
    ggplot2::geom_ribbon(
      data = posterior_data,
      ggplot2::aes(
        x = time,
        ymin = lower,
        ymax = upper,
        fill = curve
      ),
      alpha = 0.22,
      color = NA
    ) +
    ggplot2::geom_ribbon(
      data = km_data,
      ggplot2::aes(
        x = time,
        ymin = lower,
        ymax = upper,
        fill = curve
      ),
      alpha = 0.22,
      color = NA
    ) +
    ggplot2::geom_line(
      data = posterior_data,
      ggplot2::aes(
        x = time,
        y = survival,
        color = curve
      ),
      linewidth = 1.05
    ) +
    ggplot2::geom_step(
      data = km_data,
      ggplot2::aes(
        x = time,
        y = survival,
        color = curve
      ),
      linewidth = 1.05
    )

  if (is.finite(observed_median)) {
    p <- p +
      ggplot2::geom_segment(
        ggplot2::aes(
          x = x_axis_min,
          xend = observed_median,
          y = 0.5,
          yend = 0.5
        ),
        linetype = "dashed",
        linewidth = 0.7,
        color = "grey35"
      ) +
      ggplot2::geom_segment(
        ggplot2::aes(
          x = observed_median,
          xend = observed_median,
          y = 0,
          yend = 0.5
        ),
        linetype = "dashed",
        linewidth = 0.7,
        color = "grey35"
      )
  }

  p +
    ggplot2::coord_cartesian(
      xlim = c(x_axis_min, x_axis_max),
      ylim = c(0, 1)
    ) +
    ggplot2::labs(
      x = "Time",
      y = "Estimated Survival Probability",
      title = "Optional Concurrent-Control Compatibility Diagnostic"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        face = "bold"
      ),
      legend.position = "bottom",
      legend.title = ggplot2::element_blank()
    )
}


make_control_compatibility_summary_grob <- function(
    compatibility_result,
    base_size = 20
) {
  summary_df <- compatibility_result$summary

  compatibility_label <- if (
    isTRUE(summary_df$compatible[1])
  ) {
    "Compatible"
  } else {
    "Potential incompatibility detected"
  }

  table_df <- data.frame(
    `Control N` =
      summary_df$n_evaluable_control[1],
    `Mean posterior log-rank chi-square` =
      round(summary_df$posterior_logrank_mean[1], 3),
    `Bootstrap p-value` =
      round(summary_df$bootstrap_p_value_mean[1], 3),
    `Result` =
      compatibility_label,
    check.names = FALSE
  )

  gridExtra::tableGrob(
    table_df,
    rows = NULL,
    theme = gridExtra::ttheme_default(
      base_size = base_size
    )
  )
}


make_control_compatibility_combined_grob <- function(
    compatibility_result,
    x_axis_min = 0,
    x_axis_max = NULL,
    plot_theme = ggplot2::theme_bw(),
    table_base_size = 20
) {
  p <- plot_control_compatibility_check(
    compatibility_result = compatibility_result,
    x_axis_min = x_axis_min,
    x_axis_max = x_axis_max
  ) +
    plot_theme

  table_grob <- make_control_compatibility_summary_grob(
    compatibility_result = compatibility_result,
    base_size = table_base_size
  )

  gridExtra::arrangeGrob(
    p,
    table_grob,
    ncol = 1,
    heights = c(7, 1.15)
  )
}
