# Subgroup posterior survival and density summaries.
# Depends on survival.fn.lognorm() and density.fn.lognorm().

subgroup_data <- function(result,
                          input_df = NULL, cntrl, trt, burnin = 200, timepoints,
                          subgroup_list = NULL, use_pred = TRUE, conf_int = 0.95)
{
  suppressWarnings(suppressMessages(library(ggplot2)))
  suppressWarnings(suppressMessages(library(matrixStats)))
  suppressWarnings(suppressMessages(library(plyr)))
  suppressWarnings(suppressMessages(library(foreach)))
  suppressWarnings(suppressMessages(library(doParallel)))
  #### error messages plus theta extraction ----

  # make sure that there at least 5 columns in the dataframe
  if(ncol(input_df) < 5)
  {
    stop("Please check if your input dataframe has at least: \n
    a response column,
    a censoring indicator column,
    a categorical variable(s) and/or continuous variable(s) column(s),
    a treatment type column,
    a data type column")
  }

  thetas <- list()

  if(cntrl == trt)
    stop("No plot will be seen when the control is the same as the treatment")
  if(!(cntrl %in% result$Treatment_Levels))
    stop("Make sure that this control designation exists in your Treatment Levels")
  if(!(trt %in% result$Treatment_Levels))
    stop("Make sure that the treatment designation exists in your Treatment Levels")


  if(use_pred)
  {
    if(!any(names(result) %in% "Predicted_Allocation_variables"))
    {
      stop("Please make sure that you have predicted allocation variables")
    }
  } else
  {
    if(!any(names(result) %in% "Allocation_variables"))
    {
      stop("Please make sure that you have allocation variables")
    }
  }

  trt_ind <- unique(result$Treatment_Indices[which(result$Treatment_Levels == trt)])
  rwd_ind <- unique(result$Treatment_Indices[which(result$Treatment_Levels == cntrl)])


  if (use_pred) {
    alloc_mat_name <- "Predicted_Allocation_variables"
  } else {
    alloc_mat_name <- "Allocation_variables"
  }

  if (!(alloc_mat_name %in% names(result))) {
    stop("Missing ", alloc_mat_name, " in result object.")
  }

  alloc_mat <- result[[alloc_mat_name]]

  if (burnin >= nrow(alloc_mat)) {
    stop("burnin must be smaller than the number of stored MCMC iterations.")
  }

  n_pred <- ncol(alloc_mat)
  dim.mcmc <- dim(alloc_mat)

  pred_alloc_trt <- plyr::alply(
    alloc_mat[-seq_len(burnin), , drop = FALSE] + 1,
    1,
    .fun = identity
  )

  pred_alloc_rwd <- pred_alloc_trt

  # Extract posterior treatment-specific parameters for each allocated patient.
  mu1=alply( result$Lognormal_Mu_Cube[-(1:burnin), , trt_ind], 1,.fun=identity)
  mu2=alply( result$Lognormal_Mu_Cube[-(1:burnin), , rwd_ind], 1,.fun=identity)
  s1=alply( result$Lognormal_Sig_Cube[-(1:burnin), , trt_ind], 1,.fun=purrr::compose(sqrt, identity))
  s2=alply( result$Lognormal_Sig_Cube[-(1:burnin), , rwd_ind], 1,.fun=purrr::compose(sqrt, identity))
  #

  # actual values we want
  mu1 <- t(sapply(1:length(mu1), function(x) mu1[[x]][pred_alloc_trt[[x]]]))
  mu2 <- t(sapply(1:length(mu2), function(x) mu2[[x]][pred_alloc_rwd[[x]]]))
  s1 <- t(sapply(1:length(s1), function(x) s1[[x]][pred_alloc_trt[[x]]]))
  s2 <- t(sapply(1:length(s2), function(x) s2[[x]][pred_alloc_rwd[[x]]]))

  muname1 <- paste("mu", trt, sep = "_")
  muname2 <- paste("mu", cntrl, sep = "_")
  s_name1 <- paste("sd", trt, sep = "_")
  s_name2 <- paste("sd", cntrl, sep = "_")
  thetas <- list(mu1, mu2, s1, s2)
  names(thetas) <- c(muname1, muname2, s_name1, s_name2)


  if (!is.list(subgroup_list)) {
    stop(
      "subgroup_list must be a named list. For example: ",
      "list(sex = \"01\", eor = \"0\", kps = \"02\", ",
      "age = c(48, 65)). Use NA for variables that should not be restricted."
    )
  }

  expected_subgroup_vars <- c(
    result$Variable_specifications$cont_vars,
    result$Variable_specifications$cat_vars
  )

  if (
    length(subgroup_list) != length(expected_subgroup_vars) ||
    !setequal(names(subgroup_list), expected_subgroup_vars)
  ) {
    stop(
      "The names in subgroup_list must match the covariates used by the model: ",
      paste(expected_subgroup_vars, collapse = ", "),
      ". Use NA for variables that should not be restricted."
    )
  }

  if(all(is.na(subgroup_list)))
  {
    warning("Since you specified no subpopulation restrictions, we will do analysis on the whole
         data.")
  }

  ####----

  #### finding the subgroup indices ----

  # by default, start with all patients
  subgroup_ind <- seq_len(nrow(input_df))

  # if any subgroup restrictions are specified, apply them
  if (!all(is.na(subgroup_list))) {

    nonempty_names <- names(subgroup_list)[
      !vapply(subgroup_list, function(x) all(is.na(x)), logical(1))
    ]

    cat_vars <- nonempty_names[
      vapply(subgroup_list[nonempty_names], is.character, logical(1))
    ]

    cont.vars <- nonempty_names[
      !vapply(subgroup_list[nonempty_names], is.character, logical(1))
    ]

    ## ----------------------------------------------------------
    ## categorical subgroup filters
    ## subgroup_list uses original numeric codes:
    ## sex: Female = 0, Male = 1
    ## eor: GTR = 0, STR = 1, biopsy = 2
    ## kps: > 80 = 0, (60, 80] = 1, <= 60 = 2
    ##
    ## If input_df columns are factors, convert factor levels to
    ## original codes using as.integer(x) - 1.
    ## ----------------------------------------------------------
    if (length(cat_vars) > 0) {

      cat_ind_list <- lapply(cat_vars, function(v) {

        spec <- subgroup_list[[v]]
        lev_codes <- as.integer(strsplit(spec, split = "")[[1]])

        if (any(is.na(lev_codes))) {
          stop("Categorical subpopulation specification for ", v,
               " must contain only digit codes like '0', '01', or '012'.")
        }

        x <- input_df[[v]]

        if (is.factor(x) || is.character(x)) {

          if (is.null(result$Cat_Levels) || !v %in% names(result$Cat_Levels)) {
            stop(
              "Stored categorical levels were not found for variable ",
              v,
              "."
            )
          }

          expected_levels <- result$Cat_Levels[[v]]

          f <- factor(
            as.character(x),
            levels = expected_levels
          )

          if (any(is.na(f) & !is.na(x))) {
            unseen <- unique(as.character(x)[is.na(f) & !is.na(x)])

            stop(
              "Unseen categorical level(s) for ",
              v,
              ": ",
              paste(unseen, collapse = ", ")
            )
          }

          x_code <- as.integer(f) - 1L

        } else if (is.numeric(x) || is.integer(x)) {

          x_code <- as.integer(x)

        } else {

          stop(
            "Categorical variable ",
            v,
            " must be factor, character, integer, or numeric."
          )
        }

        which(x_code %in% lev_codes)
      })

      subgroup_cat <- Reduce(intersect, cat_ind_list)

    } else {
      subgroup_cat <- seq_len(nrow(input_df))
    }

    ## ----------------------------------------------------------
    ## continuous subgroup filters
    ## age interval is [min, max)
    ## ----------------------------------------------------------
    if (length(cont.vars) > 0) {

      cont_ind_list <- lapply(cont.vars, function(v) {

        rng <- as.numeric(subgroup_list[[v]])

        if (length(rng) != 2 || any(!is.finite(rng))) {
          stop("Continuous subpopulation specification for ", v,
               " must be a numeric vector of length 2, e.g. c(48, 81).")
        }

        x <- input_df[[v]]

        if (!is.numeric(x)) {
          x <- as.numeric(x)
        }

        which(x >= rng[1] & x < rng[2])
      })

      subgroup_cont <- Reduce(intersect, cont_ind_list)

    } else {
      subgroup_cont <- seq_len(nrow(input_df))
    }

    subgroup_ind <- intersect(subgroup_cat, subgroup_cont)
  }

  ####----

  #### subgroup analysis ----
  if (length(subgroup_ind) <= 1) {
    stop(
      "This subpopulation is too small for a meaningful analysis. ",
      "Selected subpopulation size = ", length(subgroup_ind), "."
    )
  }

  timepoints <- log(timepoints)

  # extract treatment indices (using grepl)
  trt_indices <- which(grepl(trt, names(thetas)) == TRUE)
  cntrl_indices <- which(grepl(cntrl, names(thetas)) == TRUE)

  cl <- ifelse((parallel::detectCores() - 1) > 1,
               parallel::detectCores() - 1,
               1)

  cl <- parallel::makeCluster(cl)
  doParallel::registerDoParallel(cl)

  on.exit({
    try(parallel::stopCluster(cl), silent = TRUE)
    try(foreach::registerDoSEQ(), silent = TRUE)
  }, add = TRUE)

  # surv_trt <- parLapply(cl = cl, 1:length(timepoints), function(i) {
  #   tmp <- pnorm(q = timepoints[i], mean = thetas[[trt_indices[1]]][ , subgroup_ind],
  #                sd = thetas[[trt_indices[2]]][, subgroup_ind], lower.tail = FALSE,
  #                log.p = TRUE)
  #
  # })


  # Parallel computation using foreach
  surv_trt <- foreach(i = 1:length(timepoints), .packages = 'stats') %dopar% {
    pnorm(q = timepoints[i],
          mean = thetas[[trt_indices[1]]][ , subgroup_ind],
          sd   = thetas[[trt_indices[2]]][ , subgroup_ind],
          lower.tail = FALSE,
          log.p = TRUE)
  }

  surv_rwd <- foreach(i = 1:length(timepoints), .packages = 'stats') %dopar% {
    pnorm(q = timepoints[i],
          mean = thetas[[cntrl_indices[1]]][ , subgroup_ind],
          sd   = thetas[[cntrl_indices[2]]][ , subgroup_ind],
          lower.tail = FALSE,
          log.p = TRUE)
  }

  # for RMST calculations
  tmp_rwd <- sapply(surv_rwd, rowLogSumExps) - log(length(subgroup_ind))
  tmp_trt <- sapply(surv_trt, rowLogSumExps) - log(length(subgroup_ind))



  f_trt <- foreach(i = 1:length(timepoints), .packages = 'stats') %dopar% {
    dnorm(x = timepoints[i],
          mean = thetas[[trt_indices[1]]][ , subgroup_ind],
          sd   = thetas[[trt_indices[2]]][ , subgroup_ind],
          log = TRUE)
  }

  f_rwd <- foreach(i = 1:length(timepoints), .packages = 'stats') %dopar% {
    dnorm(x = timepoints[i],
          mean = thetas[[cntrl_indices[1]]][ , subgroup_ind],
          sd   = thetas[[cntrl_indices[2]]][ , subgroup_ind],
          log = TRUE)
  }


  h_trt <- sapply(1:length(timepoints), function(i) {
    tmp1 <- rowLogSumExps(f_trt[[i]]) - rowLogSumExps(surv_trt[[i]]) - log(length(subgroup_ind))
  })


  h_rwd <- sapply(1:length(timepoints), function(i) {
    tmp1 <- rowLogSumExps(f_rwd[[i]]) - rowLogSumExps(surv_rwd[[i]]) - log(length(subgroup_ind))
  })

  HR <- sapply(1:length(timepoints), function(i) {
    tmp1 <- rowLogSumExps(f_trt[[i]]) - rowLogSumExps(surv_trt[[i]])
    tmp2 <- rowLogSumExps(f_rwd[[i]]) - rowLogSumExps(surv_rwd[[i]])
    tmp1 - tmp2
    # ifelse(tmp3 > max(tmp3)/4, max(tmp3)/4, tmp3)
  })

  try(parallel::stopCluster(cl), silent = TRUE)
  try(foreach::registerDoSEQ(), silent = TRUE)

  rm(f_trt, f_rwd)


  low_quant <- (1 - conf_int)/2
  high_quant <- 1 - ((1 - conf_int)/2)
  HR <- t(apply(exp(HR), 2, quantile, probs = c(low_quant, 0.5, high_quant), na.rm = TRUE))

  timepoints <- exp(timepoints)
  tmp_rwd <- exp(tmp_rwd)
  tmp_trt <- exp(tmp_trt)

  surv_rwd <- t(apply(tmp_rwd, 2, quantile,
                      probs = c(low_quant, 0.5, high_quant), na.rm = TRUE))
  surv_trt <- t(apply(tmp_trt, 2, quantile,
                      probs = c(low_quant, 0.5, high_quant), na.rm = TRUE))


  h_trt <- t(apply(h_trt, 2, quantile,
                   probs = c(low_quant, 0.5, high_quant), na.rm = TRUE))
  h_rwd <- t(apply(h_rwd, 2, quantile,
                   probs = c(low_quant, 0.5, high_quant), na.rm = TRUE))


  rm(tmp_rwd, tmp_trt)

  surv_plot_data <- data.frame(time = rep(timepoints, times = 2),
                               surv.pred.2.5 = c(surv_rwd[, 1], surv_trt[, 1]),
                               surv.pred = c(surv_rwd[, 2], surv_trt[, 2]),
                               surv.pred.97.5 = c(surv_rwd[, 3], surv_trt[, 3]))
  rm(surv_rwd, surv_trt)


  h_plot_data <- data.frame(time = rep(timepoints, times = 2),
                            h.pred.2.5 = exp(c(h_rwd[, 1], h_trt[, 1])),
                            h.pred = exp(c(h_rwd[, 2], h_trt[, 2])),
                            h.pred.97.5 = exp(c(h_rwd[, 3], h_trt[, 3])))

  # extract trt names
  trt_names <- rev(gsub("mu_", "", names(thetas[1:2])))

  # the rest of the figures (survival, density, hazard rate, hazard ratio)
  surv_plot_data$trt <- rep(trt_names, each = length(timepoints))

  surv_plot_data$trt <- factor(surv_plot_data$trt,
                               levels = c(trt_names[2], trt_names[1]))

  h_plot_data$trt <- rep(trt_names, each = length(timepoints))

  h_plot_data$trt <- factor(h_plot_data$trt,
                            levels = c(trt_names[2], trt_names[1]))


  HR_plot_data <- data.frame(time = timepoints,
                             HR.pred.2.5 = HR[, 1],
                             HR.pred = HR[, 2],
                             HR.pred.97.5 = HR[, 3])

  rm(HR, timepoints)

  p.surv <- ggplot(surv_plot_data, aes(x=time,y=surv.pred, color=trt)) +
    geom_line() + ylim(0,1) +
    ylab("Estimated Survival Probability") +
    geom_ribbon(data=surv_plot_data,aes(ymin=surv.pred.2.5,
                                        ymax=surv.pred.97.5,
                                        fill = trt),
                alpha=0.1, linetype = 0)

  p.hazarad.ratio2=ggplot(data = HR_plot_data,  aes(time,  HR.pred)) + geom_line()+
    geom_ribbon(data=HR_plot_data,aes(ymin=HR.pred.2.5,ymax=HR.pred.97.5),alpha=0.3) +
    geom_hline(yintercept = 1.0, colour = "red", linetype = 2) + ylab("Estimated Hazard ratio")

  result <- list(surv.data = surv_plot_data, p.surv = p.surv,
                 h.data = h_plot_data,
                 HR.data = HR_plot_data, p.hazarad.ratio.withCI=p.hazarad.ratio2,
                 subgroup_index = subgroup_ind,
                 thetas = thetas, control_index = cntrl_indices,
                 treatment_index = trt_indices)

  return(result)

  ####----

}
