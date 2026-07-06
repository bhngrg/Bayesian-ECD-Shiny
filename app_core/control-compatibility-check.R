# Posterior predictive concurrent-control compatibility check.
#
# This optional diagnostic compares the uploaded concurrent-control arm
# with the historical posterior predictive control distribution.

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
  
  keep <- is.finite(time_values) & !is.na(status_values)
  
  if (sum(keep) < 5) {
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
  
  if (length(observed_median) == 0) {
    observed_median <- NA_real_
  }
  
  list(
    fit = fit,
    curve = km_curve,
    observed_median = observed_median,
    n_evaluable = nrow(tmp_df)
  )
}


select_control_compatibility_draws <- function(
    result,
    burnin = 200L,
    n_posterior_draws = 200L
) {
  n_iter <- dim(result$Lognormal_Mu_Cube)[1]
  
  if (burnin >= n_iter) {
    stop("burnin must be smaller than the number of stored MCMC iterations.")
  }
  
  available_draws <- seq.int(burnin + 1L, n_iter)
  
  n_posterior_draws <- min(
    as.integer(n_posterior_draws),
    length(available_draws)
  )
  
  unique(round(seq(
    from = min(available_draws),
    to = max(available_draws),
    length.out = n_posterior_draws
  )))
}


get_control_compatibility_treatment_index <- function(
    result,
    control_label = "Control"
) {
  if (!control_label %in% result$Treatment_Levels) {
    stop(
      "Control label '", control_label,
      "' was not found in the model treatment levels. Available levels: ",
      paste(result$Treatment_Levels, collapse = ", ")
    )
  }
  
  unique(result$Treatment_Indices[
    which(result$Treatment_Levels == control_label)
  ])
}


get_control_compatibility_pi_index <- function(
    result,
    treatment_index
) {
  pi_index <- result$Combined_Indices[
    which(result$Combined_Indices[, "Treatment_Indices"] == treatment_index),
    1
  ]
  
  unique(pi_index)
}


extract_draw_weights_for_uploaded_population <- function(
    result,
    iter_index,
    pi_index
) {
  weights_i <- drop(result$picube[
    iter_index,
    ,
    pi_index,
    drop = FALSE
  ])
  
  if (is.matrix(weights_i)) {
    weights_i <- rowMeans(weights_i)
  }
  
  weights_i <- as.numeric(weights_i)
  weights_i[!is.finite(weights_i)] <- 0
  
  if (sum(weights_i) <= 0) {
    stop("Posterior mixture weights are not valid for one of the selected draws.")
  }
  
  weights_i / sum(weights_i)
}


compute_historical_predictive_control_curve <- function(
    result,
    control_label = "Control",
    time_grid = seq(1, 1251, length.out = 201),
    burnin = 200L,
    n_posterior_draws = 200L,
    conf_int = 0.95
) {
  draw_indices <- select_control_compatibility_draws(
    result = result,
    burnin = burnin,
    n_posterior_draws = n_posterior_draws
  )
  
  treatment_index <- get_control_compatibility_treatment_index(
    result = result,
    control_label = control_label
  )
  
  pi_index <- get_control_compatibility_pi_index(
    result = result,
    treatment_index = treatment_index
  )
  
  log_time_grid <- log(time_grid)
  
  survival_draw_matrix <- sapply(draw_indices, function(iter_index) {
    mu_i <- drop(result$Lognormal_Mu_Cube[
      iter_index,
      ,
      treatment_index,
      drop = FALSE
    ])
    
    sig_i <- sqrt(drop(result$Lognormal_Sig_Cube[
      iter_index,
      ,
      treatment_index,
      drop = FALSE
    ]))
    
    weights_i <- extract_draw_weights_for_uploaded_population(
      result = result,
      iter_index = iter_index,
      pi_index = pi_index
    )
    
    vapply(log_time_grid, function(log_time) {
      exp(survival.fn.lognorm(
        time = log_time,
        mu = mu_i,
        sig = sig_i,
        weights = weights_i
      ))
    }, numeric(1))
  })
  
  low_quant <- (1 - conf_int) / 2
  high_quant <- 1 - low_quant
  
  data.frame(
    time = time_grid,
    survival = apply(survival_draw_matrix, 1, stats::median, na.rm = TRUE),
    lower = apply(
      survival_draw_matrix,
      1,
      stats::quantile,
      probs = low_quant,
      na.rm = TRUE
    ),
    upper = apply(
      survival_draw_matrix,
      1,
      stats::quantile,
      probs = high_quant,
      na.rm = TRUE
    ),
    curve = "Historical posterior predictive control"
  )
}


draw_historical_predictive_control_sample_medians <- function(
    result,
    n_control,
    control_label = "Control",
    burnin = 200L,
    n_posterior_draws = 200L,
    seed = 20260706L
) {
  draw_indices <- select_control_compatibility_draws(
    result = result,
    burnin = burnin,
    n_posterior_draws = n_posterior_draws
  )
  
  treatment_index <- get_control_compatibility_treatment_index(
    result = result,
    control_label = control_label
  )
  
  pi_index <- get_control_compatibility_pi_index(
    result = result,
    treatment_index = treatment_index
  )
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  medians <- vapply(draw_indices, function(iter_index) {
    mu_i <- drop(result$Lognormal_Mu_Cube[
      iter_index,
      ,
      treatment_index,
      drop = FALSE
    ])
    
    sig_i <- sqrt(drop(result$Lognormal_Sig_Cube[
      iter_index,
      ,
      treatment_index,
      drop = FALSE
    ]))
    
    weights_i <- extract_draw_weights_for_uploaded_population(
      result = result,
      iter_index = iter_index,
      pi_index = pi_index
    )
    
    selected_clusters <- sample(
      x = seq_along(weights_i),
      size = n_control,
      replace = TRUE,
      prob = weights_i
    )
    
    predicted_times <- stats::rlnorm(
      n = n_control,
      meanlog = mu_i[selected_clusters],
      sdlog = sig_i[selected_clusters]
    )
    
    stats::median(predicted_times, na.rm = TRUE)
  }, numeric(1))
  
  data.frame(
    posterior_draw = draw_indices,
    predicted_sample_median = medians
  )
}


run_control_compatibility_check <- function(
    result,
    uploaded_data,
    time_col,
    censor_col,
    trt_col,
    control_label = "Control",
    time_grid = seq(1, 1251, length.out = 201),
    burnin = 200L,
    n_posterior_draws = 200L,
    conf_int = 0.95,
    seed = 20260706L
) {
  required_cols <- c(time_col, censor_col, trt_col)
  missing_cols <- setdiff(required_cols, names(uploaded_data))
  
  if (length(missing_cols) > 0) {
    stop(
      "Uploaded data are missing required column(s): ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  control_df <- uploaded_data[
    as.character(uploaded_data[[trt_col]]) == control_label,
    ,
    drop = FALSE
  ]
  
  if (nrow(control_df) == 0) {
    stop(
      "No concurrent-control patients were found for control label '",
      control_label,
      "'."
    )
  }
  
  km_result <- compute_uploaded_control_km_curve(
    control_df = control_df,
    time_col = time_col,
    status_col = censor_col,
    conf_int = conf_int
  )
  
  posterior_curve <- compute_historical_predictive_control_curve(
    result = result,
    control_label = control_label,
    time_grid = time_grid,
    burnin = burnin,
    n_posterior_draws = n_posterior_draws,
    conf_int = conf_int
  )
  
  posterior_sample_medians <- draw_historical_predictive_control_sample_medians(
    result = result,
    n_control = km_result$n_evaluable,
    control_label = control_label,
    burnin = burnin,
    n_posterior_draws = n_posterior_draws,
    seed = seed
  )
  
  low_quant <- (1 - conf_int) / 2
  high_quant <- 1 - low_quant
  
  median_interval <- stats::quantile(
    posterior_sample_medians$predicted_sample_median,
    probs = c(low_quant, 0.5, high_quant),
    na.rm = TRUE
  )
  
  compatible <- !is.na(km_result$observed_median) &&
    km_result$observed_median >= median_interval[[1]] &&
    km_result$observed_median <= median_interval[[3]]
  
  summary_df <- data.frame(
    control_label = control_label,
    n_control = nrow(control_df),
    n_evaluable_control = km_result$n_evaluable,
    observed_km_median = km_result$observed_median,
    posterior_predictive_median_q025 = unname(median_interval[[1]]),
    posterior_predictive_median_q500 = unname(median_interval[[2]]),
    posterior_predictive_median_q975 = unname(median_interval[[3]]),
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
    posterior_sample_medians = posterior_sample_medians,
    plot_data = plot_data
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
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      legend.position = "bottom",
      legend.title = ggplot2::element_blank()
    )
}



make_control_compatibility_summary_grob <- function(
    compatibility_result,
    base_size = 20
) {
  summary_df <- compatibility_result$summary
  
  observed_median <- summary_df$observed_km_median[1]
  
  observed_median_label <- ifelse(
    is.na(observed_median),
    "Not reached",
    as.character(round(observed_median, 2))
  )
  
  compatibility_label <- ifelse(
    is.na(observed_median),
    "Observed KM median not reached; review curve-level agreement",
    ifelse(
      isTRUE(summary_df$compatible[1]),
      "Compatible",
      "Potential incompatibility detected"
    )
  )
  
  table_df <- data.frame(
    `Control N` = summary_df$n_control[1],
    `Observed KM median` = observed_median_label,
    `Predictive median 2.5%` = round(
      summary_df$posterior_predictive_median_q025[1],
      2
    ),
    `Predictive median 50%` = round(
      summary_df$posterior_predictive_median_q500[1],
      2
    ),
    `Predictive median 97.5%` = round(
      summary_df$posterior_predictive_median_q975[1],
      2
    ),
    `Result` = compatibility_label,
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

