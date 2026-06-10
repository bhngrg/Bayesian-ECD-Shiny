# p(HR < threshold | Data)
posterior_prob_hr_lt <- function(result,
                                 input_df,
                                 cntrl,
                                 trt,
                                 timepoint_days,
                                 threshold = c(0.6, 0.8, 1.0),
                                 burnin = 200,
                                 subgroup_list = NULL,
                                 use_pred = TRUE,
                                 conf_int = 0.95) {
  
  suppressWarnings(suppressMessages(library(matrixStats)))
  suppressWarnings(suppressMessages(library(plyr)))
  
  ## ----------------------------------------------------------
  ## basic checks
  ## ----------------------------------------------------------
  if (is.null(input_df) || !is.data.frame(input_df)) {
    stop("input_df must be a data.frame.")
  }
  
  if (cntrl == trt) stop("Control and treatment must be different.")
  if (!(cntrl %in% result$Treatment_Levels)) stop("cntrl not found in Treatment_Levels.")
  if (!(trt %in% result$Treatment_Levels)) stop("trt not found in Treatment_Levels.")
  
  if (use_pred) {
    if (!("Predicted_Allocation_variables" %in% names(result))) {
      stop("Predicted_Allocation_variables not found in result.")
    }
  } else {
    if (!("Allocation_variables" %in% names(result))) {
      stop("Allocation_variables not found in result.")
    }
  }
  
  if (!("Lognormal_Mu_Cube" %in% names(result))) {
    stop("Lognormal_Mu_Cube not found in result.")
  }
  if (!("Lognormal_Sig_Cube" %in% names(result))) {
    stop("Lognormal_Sig_Cube not found in result.")
  }
  
  if (burnin < 0) stop("burnin must be >= 0.")
  
  ## ----------------------------------------------------------
  ## treatment indices
  ## ----------------------------------------------------------
  trt_ind   <- unique(result$Treatment_Indices[result$Treatment_Levels == trt])
  cntrl_ind <- unique(result$Treatment_Indices[result$Treatment_Levels == cntrl])
  
  ## ----------------------------------------------------------
  ## allocation variables
  ## ----------------------------------------------------------
  if (use_pred) {
    alloc_source <- result$Predicted_Allocation_variables
  } else {
    alloc_source <- result$Allocation_variables
  }
  
  if (nrow(alloc_source) <= burnin) {
    stop("burnin is too large relative to the number of MCMC iterations.")
  }
  
  alloc_mat <- alloc_source[-seq_len(burnin), , drop = FALSE] + 1
  
  pred_alloc_trt   <- alply(alloc_mat, 1, .fun = identity)
  pred_alloc_cntrl <- alply(alloc_mat, 1, .fun = identity)
  
  ## ----------------------------------------------------------
  ## parameter cubes
  ## ----------------------------------------------------------
  mu_trt_list <- alply(result$Lognormal_Mu_Cube[-seq_len(burnin), , trt_ind], 1, .fun = identity)
  mu_ctl_list <- alply(result$Lognormal_Mu_Cube[-seq_len(burnin), , cntrl_ind], 1, .fun = identity)
  
  sd_trt_list <- alply(result$Lognormal_Sig_Cube[-seq_len(burnin), , trt_ind], 1,
                       .fun = purrr::compose(sqrt, identity))
  sd_ctl_list <- alply(result$Lognormal_Sig_Cube[-seq_len(burnin), , cntrl_ind], 1,
                       .fun = purrr::compose(sqrt, identity))
  
  ## ----------------------------------------------------------
  ## select patient-specific cluster parameters at each MCMC draw
  ## rows = posterior draws, cols = patients
  ## ----------------------------------------------------------
  mu_trt <- t(sapply(seq_along(mu_trt_list), function(m) mu_trt_list[[m]][pred_alloc_trt[[m]]]))
  mu_ctl <- t(sapply(seq_along(mu_ctl_list), function(m) mu_ctl_list[[m]][pred_alloc_cntrl[[m]]]))
  sd_trt <- t(sapply(seq_along(sd_trt_list), function(m) sd_trt_list[[m]][pred_alloc_trt[[m]]]))
  sd_ctl <- t(sapply(seq_along(sd_ctl_list), function(m) sd_ctl_list[[m]][pred_alloc_cntrl[[m]]]))
  
  ## ----------------------------------------------------------
  ## subgroup logic
  ## subgroup_list format examples:
  ##   list(sex = "0", eor = "01", kps = NA, age = c(48, 65))
  ## categorical values refer to factor/integer level codes
  ## ----------------------------------------------------------
  subgroup_ind <- seq_len(nrow(input_df))
  
  if (!is.null(subgroup_list) && !all(is.na(subgroup_list))) {
    
    if (!is.list(subgroup_list)) {
      stop("subgroup_list must be a list or NULL.")
    }
    
    missing_vars <- setdiff(names(subgroup_list), names(input_df))
    if (length(missing_vars) > 0) {
      stop("These subgroup variables are not in input_df: ",
           paste(missing_vars, collapse = ", "))
    }
    
    nonempty_names <- names(subgroup_list)[!vapply(subgroup_list, function(x) all(is.na(x)), logical(1))]
    
    cat_vars  <- nonempty_names[vapply(subgroup_list[nonempty_names], is.character, logical(1))]
    cont_vars <- nonempty_names[!vapply(subgroup_list[nonempty_names], is.character, logical(1))]
    
    ## ----- categorical subgroup filters -----
    if (length(cat_vars) > 0) {
      
      tmp_cat <- lapply(cat_vars, function(v) {
        spec <- subgroup_list[[v]]
        
        ## split "01" -> c("0","1")
        lev_codes <- strsplit(spec, split = "")[[1]]
        
        if (length(lev_codes) == 0) {
          return(integer(0))
        }
        
        lev_codes <- as.integer(lev_codes)
        if (any(is.na(lev_codes))) {
          stop("Categorical subgroup specification for ", v,
               " must contain only digit codes like '0', '01', '012'.")
        }
        
        x <- input_df[[v]]
        
        ## interpret subgroup specification on the underlying factor/integer codes
        if (is.factor(x)) {
          x_code <- as.integer(x) - 1L
        } else if (is.numeric(x) || is.integer(x)) {
          x_code <- as.integer(x)
        } else {
          stop("Categorical subgroup variable ", v,
               " must be stored as factor, integer, or numeric.")
        }
        
        which(x_code %in% lev_codes)
      })
      
      subgroup_cat <- Reduce(intersect, tmp_cat)
    } else {
      subgroup_cat <- seq_len(nrow(input_df))
    }
    
    ## ----- continuous subgroup filters -----
    if (length(cont_vars) > 0) {
      
      tmp_cont <- lapply(cont_vars, function(v) {
        rng <- subgroup_list[[v]]
        
        if (!is.numeric(rng) || length(rng) != 2 || any(!is.finite(rng))) {
          stop("Continuous subgroup specification for ", v,
               " must be a numeric vector of length 2, e.g. c(lower, upper).")
        }
        
        x <- input_df[[v]]
        if (!is.numeric(x)) {
          stop("Continuous subgroup variable ", v, " must be numeric.")
        }
        
        which(x >= rng[1] & x < rng[2])
      })
      
      subgroup_cont <- Reduce(intersect, tmp_cont)
    } else {
      subgroup_cont <- seq_len(nrow(input_df))
    }
    
    subgroup_ind <- intersect(subgroup_cat, subgroup_cont)
  }
  
  if (length(subgroup_ind) <= 1) {
    stop("This subgroup is too small for a meaningful analysis.")
  }
  
  ## ----------------------------------------------------------
  ## one timepoint only
  ## ----------------------------------------------------------
  log_time <- log(timepoint_days)
  
  ## log survival for each patient, each posterior draw
  surv_trt_log <- pnorm(q = log_time,
                        mean = mu_trt[, subgroup_ind, drop = FALSE],
                        sd   = sd_trt[, subgroup_ind, drop = FALSE],
                        lower.tail = FALSE,
                        log.p = TRUE)
  
  surv_ctl_log <- pnorm(q = log_time,
                        mean = mu_ctl[, subgroup_ind, drop = FALSE],
                        sd   = sd_ctl[, subgroup_ind, drop = FALSE],
                        lower.tail = FALSE,
                        log.p = TRUE)
  
  ## log density for each patient, each posterior draw
  f_trt_log <- dnorm(x = log_time,
                     mean = mu_trt[, subgroup_ind, drop = FALSE],
                     sd   = sd_trt[, subgroup_ind, drop = FALSE],
                     log = TRUE)
  
  f_ctl_log <- dnorm(x = log_time,
                     mean = mu_ctl[, subgroup_ind, drop = FALSE],
                     sd   = sd_ctl[, subgroup_ind, drop = FALSE],
                     log = TRUE)
  
  ## subgroup-average log hazard for each posterior draw
  log_h_trt <- rowLogSumExps(f_trt_log) - rowLogSumExps(surv_trt_log)
  log_h_ctl <- rowLogSumExps(f_ctl_log) - rowLogSumExps(surv_ctl_log)
  
  ## posterior draws of HR(t)
  hr_draws <- exp(log_h_trt - log_h_ctl)
  
  ## posterior probabilities for multiple thresholds
  threshold <- as.numeric(threshold)
  
  post_prob <- sapply(threshold, function(th) {
    mean(hr_draws < th, na.rm = TRUE)
  })
  
  post_prob <- matrix(post_prob, nrow = 1)
  colnames(post_prob) <- paste0("HR_lt_", threshold)
  rownames(post_prob) <- paste0("t_", timepoint_days)
  
  ## posterior summary
  low_quant  <- (1 - conf_int) / 2
  high_quant <- 1 - low_quant
  hr_summary <- quantile(hr_draws,
                         probs = c(low_quant, 0.5, high_quant),
                         na.rm = TRUE)
  
  out <- list(
    timepoint_days = timepoint_days,
    threshold = threshold,
    posterior_probability = post_prob,
    hr_draws = hr_draws,
    hr_summary = hr_summary,
    subgroup_index = subgroup_ind
  )
  
  return(out)
}

# p(RMST_ratio > threshold | Data)
# p(RMST_difference > threshold | Data)
posterior_prob_rmst <- function(result,
                                input_df,
                                cntrl,
                                trt,
                                time_horizons,
                                ratio_threshold = c(1.0, 1.2, 1.4),
                                diff_threshold  = c(0, 60, 120),
                                burnin = 200,
                                subgroup_list = NULL,
                                use_pred = TRUE,
                                conf_int = 0.95,
                                n_time_grid = 200) {
  
  suppressWarnings(suppressMessages(library(matrixStats)))
  suppressWarnings(suppressMessages(library(plyr)))
  
  ## ----------------------------------------------------------
  ## basic checks
  ## ----------------------------------------------------------
  if (is.null(input_df) || !is.data.frame(input_df)) {
    stop("input_df must be a data.frame.")
  }
  
  if (cntrl == trt) stop("Control and treatment must be different.")
  if (!(cntrl %in% result$Treatment_Levels)) stop("cntrl not found in Treatment_Levels.")
  if (!(trt %in% result$Treatment_Levels)) stop("trt not found in Treatment_Levels.")
  
  if (length(time_horizons) < 1) {
    stop("Please provide at least one time horizon.")
  }
  
  time_horizons <- as.numeric(time_horizons)
  if (any(!is.finite(time_horizons)) || any(time_horizons <= 0)) {
    stop("All time_horizons must be positive finite numbers.")
  }
  time_horizons <- sort(unique(time_horizons))
  
  if (!is.numeric(n_time_grid) || length(n_time_grid) != 1 || n_time_grid < 2) {
    stop("n_time_grid must be a single integer >= 2.")
  }
  n_time_grid <- as.integer(n_time_grid)
  
  if (burnin < 0) stop("burnin must be >= 0.")
  
  if (use_pred) {
    if (!("Predicted_Allocation_variables" %in% names(result))) {
      stop("Predicted_Allocation_variables not found in result.")
    }
  } else {
    if (!("Allocation_variables" %in% names(result))) {
      stop("Allocation_variables not found in result.")
    }
  }
  
  if (!("Lognormal_Mu_Cube" %in% names(result))) {
    stop("Lognormal_Mu_Cube not found in result.")
  }
  if (!("Lognormal_Sig_Cube" %in% names(result))) {
    stop("Lognormal_Sig_Cube not found in result.")
  }
  
  ## ----------------------------------------------------------
  ## treatment indices
  ## ----------------------------------------------------------
  trt_ind   <- unique(result$Treatment_Indices[result$Treatment_Levels == trt])
  cntrl_ind <- unique(result$Treatment_Indices[result$Treatment_Levels == cntrl])
  
  ## ----------------------------------------------------------
  ## allocation variables
  ## ----------------------------------------------------------
  if (use_pred) {
    alloc_source <- result$Predicted_Allocation_variables
  } else {
    alloc_source <- result$Allocation_variables
  }
  
  if (nrow(alloc_source) <= burnin) {
    stop("burnin is too large relative to the number of MCMC iterations.")
  }
  
  if (burnin == 0) {
    alloc_mat <- alloc_source[, , drop = FALSE] + 1
  } else {
    alloc_mat <- alloc_source[-seq_len(burnin), , drop = FALSE] + 1
  }
  
  pred_alloc_trt   <- alply(alloc_mat, 1, .fun = identity)
  pred_alloc_cntrl <- alply(alloc_mat, 1, .fun = identity)
  
  ## ----------------------------------------------------------
  ## parameter cubes
  ## ----------------------------------------------------------
  if (burnin == 0) {
    mu_trt_list <- alply(result$Lognormal_Mu_Cube[, , trt_ind], 1, .fun = identity)
    mu_ctl_list <- alply(result$Lognormal_Mu_Cube[, , cntrl_ind], 1, .fun = identity)
    
    sd_trt_list <- alply(result$Lognormal_Sig_Cube[, , trt_ind], 1,
                         .fun = purrr::compose(sqrt, identity))
    sd_ctl_list <- alply(result$Lognormal_Sig_Cube[, , cntrl_ind], 1,
                         .fun = purrr::compose(sqrt, identity))
  } else {
    mu_trt_list <- alply(result$Lognormal_Mu_Cube[-seq_len(burnin), , trt_ind], 1, .fun = identity)
    mu_ctl_list <- alply(result$Lognormal_Mu_Cube[-seq_len(burnin), , cntrl_ind], 1, .fun = identity)
    
    sd_trt_list <- alply(result$Lognormal_Sig_Cube[-seq_len(burnin), , trt_ind], 1,
                         .fun = purrr::compose(sqrt, identity))
    sd_ctl_list <- alply(result$Lognormal_Sig_Cube[-seq_len(burnin), , cntrl_ind], 1,
                         .fun = purrr::compose(sqrt, identity))
  }
  
  ## ----------------------------------------------------------
  ## select patient-specific cluster parameters at each MCMC draw
  ## rows = posterior draws, cols = patients
  ## ----------------------------------------------------------
  mu_trt <- t(sapply(seq_along(mu_trt_list), function(m) mu_trt_list[[m]][pred_alloc_trt[[m]]]))
  mu_ctl <- t(sapply(seq_along(mu_ctl_list), function(m) mu_ctl_list[[m]][pred_alloc_cntrl[[m]]]))
  sd_trt <- t(sapply(seq_along(sd_trt_list), function(m) sd_trt_list[[m]][pred_alloc_trt[[m]]]))
  sd_ctl <- t(sapply(seq_along(sd_ctl_list), function(m) sd_ctl_list[[m]][pred_alloc_cntrl[[m]]]))
  
  ## ----------------------------------------------------------
  ## subgroup logic
  ## subgroup_list format examples:
  ##   list(sex = "0", eor = "01", kps = NA, age = c(48, 65))
  ## categorical values refer to factor/integer level codes
  ## ----------------------------------------------------------
  subgroup_ind <- seq_len(nrow(input_df))
  
  if (!is.null(subgroup_list) && !all(is.na(subgroup_list))) {
    
    if (!is.list(subgroup_list)) {
      stop("subgroup_list must be a list or NULL.")
    }
    
    missing_vars <- setdiff(names(subgroup_list), names(input_df))
    if (length(missing_vars) > 0) {
      stop("These subgroup variables are not in input_df: ",
           paste(missing_vars, collapse = ", "))
    }
    
    nonempty_names <- names(subgroup_list)[!vapply(subgroup_list, function(x) all(is.na(x)), logical(1))]
    
    cat_vars  <- nonempty_names[vapply(subgroup_list[nonempty_names], is.character, logical(1))]
    cont_vars <- nonempty_names[!vapply(subgroup_list[nonempty_names], is.character, logical(1))]
    
    ## ----- categorical subgroup filters -----
    if (length(cat_vars) > 0) {
      
      tmp_cat <- lapply(cat_vars, function(v) {
        spec <- subgroup_list[[v]]
        
        lev_codes <- strsplit(spec, split = "")[[1]]
        if (length(lev_codes) == 0) {
          return(integer(0))
        }
        
        lev_codes <- as.integer(lev_codes)
        if (any(is.na(lev_codes))) {
          stop("Categorical subgroup specification for ", v,
               " must contain only digit codes like '0', '01', '012'.")
        }
        
        x <- input_df[[v]]
        
        if (is.factor(x)) {
          x_code <- as.integer(x) - 1L
        } else if (is.numeric(x) || is.integer(x)) {
          x_code <- as.integer(x)
        } else {
          stop("Categorical subgroup variable ", v,
               " must be stored as factor, integer, or numeric.")
        }
        
        which(x_code %in% lev_codes)
      })
      
      subgroup_cat <- Reduce(intersect, tmp_cat)
    } else {
      subgroup_cat <- seq_len(nrow(input_df))
    }
    
    ## ----- continuous subgroup filters -----
    if (length(cont_vars) > 0) {
      
      tmp_cont <- lapply(cont_vars, function(v) {
        rng <- subgroup_list[[v]]
        
        if (!is.numeric(rng) || length(rng) != 2 || any(!is.finite(rng))) {
          stop("Continuous subgroup specification for ", v,
               " must be a numeric vector of length 2, e.g. c(lower, upper).")
        }
        
        x <- input_df[[v]]
        if (!is.numeric(x)) {
          stop("Continuous subgroup variable ", v, " must be numeric.")
        }
        
        which(x >= rng[1] & x < rng[2])
      })
      
      subgroup_cont <- Reduce(intersect, tmp_cont)
    } else {
      subgroup_cont <- seq_len(nrow(input_df))
    }
    
    subgroup_ind <- intersect(subgroup_cat, subgroup_cont)
  }
  
  if (length(subgroup_ind) <= 1) {
    stop("This subgroup is too small for a meaningful analysis.")
  }
  
  ## ----------------------------------------------------------
  ## helper: trapezoidal RMST per draw
  ## S_mat: rows = posterior draws, cols = time grid
  ## ----------------------------------------------------------
  rmst_trapz_draws <- function(time_vec, S_mat) {
    dt <- diff(time_vec)
    left  <- S_mat[, -ncol(S_mat), drop = FALSE]
    right <- S_mat[, -1,           drop = FALSE]
    rowSums(((left + right) / 2) *
              matrix(dt, nrow = nrow(S_mat), ncol = length(dt), byrow = TRUE))
  }
  
  low_quant  <- (1 - conf_int) / 2
  high_quant <- 1 - low_quant
  
  ratio_threshold <- as.numeric(ratio_threshold)
  diff_threshold  <- as.numeric(diff_threshold)
  
  ## ----------------------------------------------------------
  ## loop over horizons
  ## ----------------------------------------------------------
  out_by_horizon <- lapply(time_horizons, function(tau) {
    
    time_grid <- seq(from = .Machine$double.eps, to = tau, length.out = n_time_grid)
    log_time_grid <- log(time_grid)
    
    ## log survival matrices by time point
    surv_trt_log_list <- lapply(seq_along(log_time_grid), function(k) {
      pnorm(q = log_time_grid[k],
            mean = mu_trt[, subgroup_ind, drop = FALSE],
            sd   = sd_trt[, subgroup_ind, drop = FALSE],
            lower.tail = FALSE,
            log.p = TRUE)
    })
    
    surv_ctl_log_list <- lapply(seq_along(log_time_grid), function(k) {
      pnorm(q = log_time_grid[k],
            mean = mu_ctl[, subgroup_ind, drop = FALSE],
            sd   = sd_ctl[, subgroup_ind, drop = FALSE],
            lower.tail = FALSE,
            log.p = TRUE)
    })
    
    ## subgroup-average survival curves for each posterior draw
    S_trt_draws <- exp(sapply(surv_trt_log_list, rowLogSumExps) - log(length(subgroup_ind)))
    S_ctl_draws <- exp(sapply(surv_ctl_log_list, rowLogSumExps) - log(length(subgroup_ind)))
    
    S_trt_draws <- as.matrix(S_trt_draws)
    S_ctl_draws <- as.matrix(S_ctl_draws)
    
    ## RMST draws
    rmst_trt_draws <- rmst_trapz_draws(time_grid, S_trt_draws)
    rmst_ctl_draws <- rmst_trapz_draws(time_grid, S_ctl_draws)
    
    rmst_diff_draws  <- rmst_trt_draws - rmst_ctl_draws
    rmst_ratio_draws <- rmst_trt_draws / rmst_ctl_draws
    
    ## posterior probabilities
    post_prob_diff <- sapply(diff_threshold, function(th) {
      mean(rmst_diff_draws > th, na.rm = TRUE)
    })
    
    post_prob_ratio <- sapply(ratio_threshold, function(th) {
      mean(rmst_ratio_draws > th, na.rm = TRUE)
    })
    
    post_prob_diff <- matrix(post_prob_diff, nrow = 1)
    colnames(post_prob_diff) <- paste0("RMST_diff_gt_", diff_threshold)
    rownames(post_prob_diff) <- paste0("tau_", tau)
    
    post_prob_ratio <- matrix(post_prob_ratio, nrow = 1)
    colnames(post_prob_ratio) <- paste0("RMST_ratio_gt_", ratio_threshold)
    rownames(post_prob_ratio) <- paste0("tau_", tau)
    
    ## posterior summaries
    rmst_trt_summary <- quantile(rmst_trt_draws,
                                 probs = c(low_quant, 0.5, high_quant),
                                 na.rm = TRUE)
    
    rmst_ctl_summary <- quantile(rmst_ctl_draws,
                                 probs = c(low_quant, 0.5, high_quant),
                                 na.rm = TRUE)
    
    rmst_diff_summary <- quantile(rmst_diff_draws,
                                  probs = c(low_quant, 0.5, high_quant),
                                  na.rm = TRUE)
    
    rmst_ratio_summary <- quantile(rmst_ratio_draws,
                                   probs = c(low_quant, 0.5, high_quant),
                                   na.rm = TRUE)
    
    summary_table <- rbind(
      Treatment  = rmst_trt_summary,
      Control    = rmst_ctl_summary,
      Difference = rmst_diff_summary,
      Ratio      = rmst_ratio_summary
    )
    colnames(summary_table) <- c(paste0(100 * low_quant, "%"),
                                 "50%",
                                 paste0(100 * high_quant, "%"))
    
    list(
      tau_days = tau,
      time_grid = time_grid,
      posterior_probability_diff = post_prob_diff,
      posterior_probability_ratio = post_prob_ratio,
      rmst_trt_draws = rmst_trt_draws,
      rmst_ctl_draws = rmst_ctl_draws,
      rmst_diff_draws = rmst_diff_draws,
      rmst_ratio_draws = rmst_ratio_draws,
      rmst_trt_summary = rmst_trt_summary,
      rmst_ctl_summary = rmst_ctl_summary,
      rmst_diff_summary = rmst_diff_summary,
      rmst_ratio_summary = rmst_ratio_summary,
      summary_table = summary_table
    )
  })
  
  names(out_by_horizon) <- paste0(time_horizons, "_days")
  
  out <- list(
    treatment = trt,
    control = cntrl,
    time_horizons = time_horizons,
    ratio_threshold = ratio_threshold,
    diff_threshold = diff_threshold,
    subgroup_index = subgroup_ind,
    results_by_horizon = out_by_horizon
  )
  
  return(out)
}