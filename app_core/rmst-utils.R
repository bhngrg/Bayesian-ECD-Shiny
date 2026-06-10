# RMST helper and plotting functions.
# Depends on survival.fn.lognorm() and subgroup_data().
# Uses foreach + doParallel with explicit function exports for cross-platform workers.

rmst_lognormal <- function(result, cntrl, trt, t_star, burnin = 200) {
  
  suppressWarnings(suppressMessages(library(ggplot2)))
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
  
  timepoints <- log(seq(.Machine$double.neg.eps, t_star, length.out = 200))
  
  ncores <- max(parallel::detectCores() - 1, 1)
  cl <- parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl)
  
  on.exit({
    try(parallel::stopCluster(cl), silent = TRUE)
  }, add = TRUE)
  
  logsurv1 <- foreach(
    i = seq_along(mu1),
    .combine = "cbind",
    .export = c("survival.fn.lognorm"),
    .packages = c("stats", "matrixStats")
  ) %dopar% {
    mu_i <- mu1[[i]]
    s_i <- s1[[i]]
    pi_i <- pi1[[i]]
    
    sapply(timepoints, function(t) {
      survival.fn.lognorm(t, mu_i, s_i, pi_i)
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
    
    sapply(timepoints, function(t) {
      survival.fn.lognorm(t, mu_i, s_i, pi_i)
    })
  }
  
  rmst_calc <- function(x, y, xmin, xmax) {
    out <- integrate(
      f = approxfun(x, y, rule = 2),
      lower = xmin,
      upper = xmax
    )
    
    out$value
  }
  
  rmst1 <- foreach(
    i = seq_len(ncol(logsurv1)),
    .combine = "c",
    .packages = "stats",
    .export = "rmst_calc"
  ) %dopar% {
    rmst_calc(
      x = exp(timepoints),
      y = exp(logsurv1)[, i],
      xmin = 0,
      xmax = t_star
    )
  }
  
  rmst2 <- foreach(
    i = seq_len(ncol(logsurv2)),
    .combine = "c",
    .packages = "stats",
    .export = "rmst_calc"
  ) %dopar% {
    rmst_calc(
      x = exp(timepoints),
      y = exp(logsurv2)[, i],
      xmin = 0,
      xmax = t_star
    )
  }
  
  surv1.median <- exp(apply(logsurv1, 1, stats::median))
  surv2.median <- exp(apply(logsurv2, 1, stats::median))
  
  surv1.dat <- data.frame(
    survival = surv1.median,
    Time = exp(timepoints),
    Arm = trt
  )
  
  surv2.dat <- data.frame(
    survival = surv2.median,
    Time = exp(timepoints),
    Arm = cntrl
  )
  
  rmst_diff <- rmst1 - rmst2
  rmst_ratio <- rmst1 / rmst2
  
  rmst.result <- data.frame(
    Drug.Type = c(trt, cntrl),
    RMST.Days = c(stats::median(rmst1), stats::median(rmst2))
  )
  
  names(rmst.result) <- c("Drug Type", "RMST (Days)")
  
  rmst.diff <- data.frame(
    quantile = c("2.5%", "50%", "97.5%"),
    Difference.Days = c(
      stats::quantile(rmst_diff, probs = 0.025, na.rm = TRUE),
      stats::median(rmst_diff, na.rm = TRUE),
      stats::quantile(rmst_diff, probs = 0.975, na.rm = TRUE)
    )
  )
  
  names(rmst.diff) <- c("Quantile", "Difference (Days)")
  
  rmst.ratio <- data.frame(
    quantile = c("2.5%", "50%", "97.5%"),
    Ratio = c(
      stats::quantile(rmst_ratio, probs = 0.025, na.rm = TRUE),
      stats::median(rmst_ratio, na.rm = TRUE),
      stats::quantile(rmst_ratio, probs = 0.975, na.rm = TRUE)
    )
  )
  
  names(rmst.ratio) <- c("Quantile", "Ratio")
  
  p.surv.trt <- ggplot2::ggplot(
    surv1.dat,
    ggplot2::aes(x = Time, y = survival)
  ) +
    ggplot2::geom_line(color = "#F8766D") +
    ggplot2::ylim(0, 1) +
    ggplot2::xlim(0, t_star + 100) +
    ggplot2::ylab("Estimated Survival Probability") +
    ggplot2::geom_ribbon(
      data = surv1.dat,
      ggplot2::aes(ymin = 0, ymax = survival),
      alpha = 0.2,
      linetype = 0,
      fill = "#F8766D"
    )
  
  p.surv.cntrl <- ggplot2::ggplot(
    surv2.dat,
    ggplot2::aes(x = Time, y = survival)
  ) +
    ggplot2::geom_line(color = "#00BFC4") +
    ggplot2::ylim(0, 1) +
    ggplot2::xlim(0, t_star + 100) +
    ggplot2::ylab("Estimated Survival Probability") +
    ggplot2::geom_ribbon(
      data = surv2.dat,
      ggplot2::aes(ymin = 0, ymax = survival),
      alpha = 0.2,
      linetype = 0,
      fill = "#00BFC4"
    )
  
  list(
    p.surv.trt = p.surv.trt,
    p.surv.cntrl = p.surv.cntrl,
    RMST = rmst.result,
    Difference = rmst.diff,
    Ratio = rmst.ratio
  )
}


RMST_result <- function(result_subgroup_data = list(), time_horizons = c()) {
  
  suppressWarnings(suppressMessages(library(matrixStats)))
  suppressWarnings(suppressMessages(library(foreach)))
  suppressWarnings(suppressMessages(library(doParallel)))
  
  subgroup_ind <- result_subgroup_data$subgroup_index
  thetas <- result_subgroup_data$thetas
  trt_indices <- result_subgroup_data$treatment_index
  cntrl_indices <- result_subgroup_data$control_index
  
  trt_names <- rev(gsub("mu_", "", names(thetas[1:2])))
  
  timepoints <- lapply(seq_along(time_horizons), function(i) {
    log(seq(
      from = .Machine$double.neg.eps,
      to = time_horizons[i],
      length.out = 200
    ))
  })
  
  ncores <- max(parallel::detectCores() - 1, 1)
  cl <- parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl)
  
  on.exit({
    try(parallel::stopCluster(cl), silent = TRUE)
  }, add = TRUE)
  
  RMST_results <- foreach(
    i = seq_along(timepoints),
    .packages = c("stats", "matrixStats")
  ) %dopar% {
    
    surv_trt <- lapply(seq_along(timepoints[[i]]), function(j) {
      stats::pnorm(
        q = timepoints[[i]][j],
        mean = thetas[[trt_indices[1]]][, subgroup_ind],
        sd = thetas[[trt_indices[2]]][, subgroup_ind],
        lower.tail = FALSE,
        log.p = TRUE
      )
    })
    
    surv_cntrl <- lapply(seq_along(timepoints[[i]]), function(j) {
      stats::pnorm(
        q = timepoints[[i]][j],
        mean = thetas[[cntrl_indices[1]]][, subgroup_ind],
        sd = thetas[[cntrl_indices[2]]][, subgroup_ind],
        lower.tail = FALSE,
        log.p = TRUE
      )
    })
    
    tmp_cntrl <- exp(
      sapply(surv_cntrl, matrixStats::rowLogSumExps) -
        log(length(subgroup_ind))
    )
    
    tmp_trt <- exp(
      sapply(surv_trt, matrixStats::rowLogSumExps) -
        log(length(subgroup_ind))
    )
    
    timepoints_exp <- exp(timepoints[[i]])
    
    rmst_calc <- function(x, y, xmin, xmax) {
      out <- integrate(
        f = approxfun(x, y, rule = 2),
        lower = xmin,
        upper = xmax
      )
      
      out$value
    }
    
    RMST_cntrl <- sapply(seq_len(nrow(tmp_cntrl)), function(j) {
      rmst_calc(
        x = timepoints_exp,
        y = tmp_cntrl[j, ],
        xmin = 0,
        xmax = time_horizons[i]
      )
    })
    
    RMST_trt <- sapply(seq_len(nrow(tmp_trt)), function(j) {
      rmst_calc(
        x = timepoints_exp,
        y = tmp_trt[j, ],
        xmin = 0,
        xmax = time_horizons[i]
      )
    })
    
    RMST_diff <- RMST_trt - RMST_cntrl
    RMST_ratio <- RMST_trt / RMST_cntrl
    
    matrix(
      c(
        stats::quantile(RMST_trt, probs = c(0.025, 0.5, 0.975), na.rm = TRUE),
        stats::quantile(RMST_cntrl, probs = c(0.025, 0.5, 0.975), na.rm = TRUE),
        stats::quantile(RMST_diff, probs = c(0.025, 0.5, 0.975), na.rm = TRUE),
        stats::quantile(RMST_ratio, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
      ),
      nrow = 4,
      ncol = 3,
      byrow = TRUE,
      dimnames = list(
        c(trt_names[2], trt_names[1], "Difference", "Ratio"),
        c("2.5%", "50%", "97.5%")
      )
    )
  }
  
  names(RMST_results) <- paste(time_horizons, "days", sep = " ")
  
  RMST_results
}