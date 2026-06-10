# Plotting function for posterior lognormal survival and density summaries.
# Depends on survival.fn.lognorm() and density.fn.lognorm().

plots_lognormal <- function(result, cntrl, trt, timepoints, burnin = 1e3,
                            conf_int = 0.95) {
  
  suppressWarnings(suppressMessages(library(foreach)))
  suppressWarnings(suppressMessages(library(doParallel)))
  
  if (cntrl == trt) {
    stop("No plot will be seen when the control is the same as the treatment.")
  }
  
  if (!(cntrl %in% result$Treatment_Levels)) {
    stop("Make sure that this control designation exists in your Treatment Levels.")
  }
  
  if (!(trt %in% result$Treatment_Levels)) {
    stop("Make sure that the treatment designation exists in your Treatment Levels.")
  }
  
  if (burnin >= dim(result$Lognormal_Mu_Cube)[1]) {
    stop("burnin must be smaller than the number of stored MCMC iterations.")
  }
  
  trt_ind <- unique(result$Treatment_Indices[which(result$Treatment_Levels == trt)])
  cntrl_ind <- unique(result$Treatment_Indices[which(result$Treatment_Levels == cntrl)])
  
  keep_iter <- seq.int(burnin + 1, dim(result$Lognormal_Mu_Cube)[1])
  
  make_draw_list <- function(x) {
    lapply(seq_len(dim(x)[1]), function(i) {
      drop(x[i, , , drop = FALSE])
    })
  }
  
  mu1 <- make_draw_list(result$Lognormal_Mu_Cube[keep_iter, , trt_ind, drop = FALSE])
  mu2 <- make_draw_list(result$Lognormal_Mu_Cube[keep_iter, , cntrl_ind, drop = FALSE])
  
  s1 <- lapply(
    make_draw_list(result$Lognormal_Sig_Cube[keep_iter, , trt_ind, drop = FALSE]),
    sqrt
  )
  
  s2 <- lapply(
    make_draw_list(result$Lognormal_Sig_Cube[keep_iter, , cntrl_ind, drop = FALSE]),
    sqrt
  )
  
  pi1_ind <- result$Combined_Indices[
    which(result$Combined_Indices[, "Treatment_Indices"] == trt_ind),
    1
  ]
  
  pi2_ind <- result$Combined_Indices[
    which(result$Combined_Indices[, "Treatment_Indices"] == cntrl_ind),
    1
  ]
  
  pi1 <- make_draw_list(result$picube[keep_iter, , pi1_ind, drop = FALSE])
  pi2 <- make_draw_list(result$picube[keep_iter, , pi2_ind, drop = FALSE])
  
  if (is.matrix(pi1[[1]])) {
    pi1 <- lapply(pi1, rowMeans)
  }
  
  if (is.matrix(pi2[[1]])) {
    pi2 <- lapply(pi2, rowMeans)
  }
  
  log_timepoints <- log(timepoints)
  
  # ------------------------------------------------------------
  # Parallel evaluation
  # Important:
  # - Do NOT use .packages = "ExtendedCAPPMx"
  # - Export local sourced functions to workers instead
  # ------------------------------------------------------------
  ncores <- max(parallel::detectCores() - 1, 1)
  cl <- parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl)
  
  on.exit({
    try(parallel::stopCluster(cl), silent = TRUE)
  }, add = TRUE)
  
  logdens1 <- foreach(
    i = seq_along(mu1),
    .combine = "cbind",
    .export = c("density.fn.lognorm"),
    .packages = c("stats", "matrixStats")
  ) %dopar% {
    mu_i <- mu1[[i]]
    s_i <- s1[[i]]
    pi_i <- pi1[[i]]
    
    sapply(log_timepoints, function(t) {
      density.fn.lognorm(t, mu_i, s_i, pi_i)
    })
  }
  
  logsurv1 <- foreach(
    i = seq_along(mu1),
    .combine = "cbind",
    .export = c("survival.fn.lognorm"),
    .packages = c("stats", "matrixStats")
  ) %dopar% {
    mu_i <- mu1[[i]]
    s_i <- s1[[i]]
    pi_i <- pi1[[i]]
    
    sapply(log_timepoints, function(t) {
      survival.fn.lognorm(t, mu_i, s_i, pi_i)
    })
  }
  
  logdens2 <- foreach(
    i = seq_along(mu2),
    .combine = "cbind",
    .export = c("density.fn.lognorm"),
    .packages = c("stats", "matrixStats")
  ) %dopar% {
    mu_i <- mu2[[i]]
    s_i <- s2[[i]]
    pi_i <- pi2[[i]]
    
    sapply(log_timepoints, function(t) {
      density.fn.lognorm(t, mu_i, s_i, pi_i)
    })
  }
  
  logsurv2 <- foreach(
    i = seq_along(mu2),
    .combine = "cbind",
    .export = c("survival.fn.lognorm"),
    .packages = c("stats", "matrixStats")
  ) %dopar% {
    mu_i <- mu2[[i]]
    s_i <- s2[[i]]
    pi_i <- pi2[[i]]
    
    sapply(log_timepoints, function(t) {
      survival.fn.lognorm(t, mu_i, s_i, pi_i)
    })
  }
  
  try(parallel::stopCluster(cl), silent = TRUE)
  try(foreach::registerDoSEQ(), silent = TRUE)
  
  # Hazard rate plot
  log.hr1 <- apply(logdens1 - logsurv1, 1, stats::median)
  log.hr2 <- apply(logdens2 - logsurv2, 1, stats::median)
  
  hr1 <- exp(log.hr1)
  hr2 <- exp(log.hr2)
  
  hr1.dat <- data.frame(
    HR = hr1,
    Time = timepoints,
    Arm = trt
  )
  
  hr2.dat <- data.frame(
    HR = hr2,
    Time = timepoints,
    Arm = cntrl
  )
  
  hr.mean.compare <- rbind(hr1.dat, hr2.dat)
  hr.mean.compare$Arm <- factor(hr.mean.compare$Arm, levels = c(trt, cntrl))
  
  dat.hrplot <- reshape2::melt(
    hr.mean.compare,
    id.vars = c("Time", "Arm")
  )
  
  p.hazard.rate <- ggplot2::ggplot(
    dat.hrplot,
    ggplot2::aes(x = Time, y = value, color = Arm)
  ) +
    ggplot2::geom_line() +
    ggplot2::ylab("Hazard rate")
  
  # Quantiles
  low_quant <- (1 - conf_int) / 2
  high_quant <- 1 - ((1 - conf_int) / 2)
  
  # Survival curves
  surv1.median <- exp(apply(logsurv1, 1, stats::median))
  surv2.median <- exp(apply(logsurv2, 1, stats::median))
  
  surv1.qtiles <- exp(apply(
    logsurv1,
    1,
    stats::quantile,
    probs = c(low_quant, high_quant)
  ))
  
  surv2.qtiles <- exp(apply(
    logsurv2,
    1,
    stats::quantile,
    probs = c(low_quant, high_quant)
  ))
  
  surv1.dat <- data.frame(
    survival = surv1.median,
    lower = surv1.qtiles[1, ],
    upper = surv1.qtiles[2, ],
    Time = timepoints,
    Arm = trt
  )
  
  surv2.dat <- data.frame(
    survival = surv2.median,
    lower = surv2.qtiles[1, ],
    upper = surv2.qtiles[2, ],
    Time = timepoints,
    Arm = cntrl
  )
  
  surv.median.compare <- rbind(surv1.dat, surv2.dat)
  surv.median.compare$Arm <- factor(
    surv.median.compare$Arm,
    levels = c(trt, cntrl)
  )
  
  p.surv <- ggplot2::ggplot(
    surv.median.compare,
    ggplot2::aes(x = Time, y = survival, color = Arm)
  ) +
    ggplot2::geom_line() +
    ggplot2::ylim(0, 1) +
    ggplot2::ylab("Estimated Survival Probability") +
    ggplot2::geom_ribbon(
      data = surv.median.compare,
      ggplot2::aes(ymin = lower, ymax = upper, fill = Arm),
      alpha = 0.3,
      linetype = 0
    )
  
  # Hazard ratio with credible interval
  tmp.log.HR <- (logdens1 - logsurv1) - (logdens2 - logsurv2)
  
  log.HR.mean <- apply(tmp.log.HR, 1, stats::median)
  
  log.HR.qtiles <- apply(
    tmp.log.HR,
    1,
    stats::quantile,
    probs = c(low_quant, high_quant)
  )
  
  HR.mean <- exp(log.HR.mean)
  HR.qtiles <- exp(log.HR.qtiles)
  
  hazard.ratio <- data.frame(
    HR = HR.mean,
    lower = HR.qtiles[1, ],
    upper = HR.qtiles[2, ],
    Time = timepoints
  )
  
  p.hazarad.ratio2 <- ggplot2::ggplot(
    data = hazard.ratio,
    ggplot2::aes(Time, HR)
  ) +
    ggplot2::geom_line() +
    ggplot2::geom_ribbon(
      data = hazard.ratio,
      ggplot2::aes(ymin = lower, ymax = upper),
      alpha = 0.3
    ) +
    ggplot2::geom_hline(
      yintercept = 1.0,
      colour = "red",
      linetype = 2
    ) +
    ggplot2::ylab("Estimated Hazard ratio")
  
  list(
    p.surv = p.surv,
    p.hazarad = p.hazard.rate,
    p.hazarad.ratio.withCI = p.hazarad.ratio2,
    Hazard.Ratio = exp(tmp.log.HR),
    surv.data = surv.median.compare,
    HR.data = hazard.ratio
  )
}