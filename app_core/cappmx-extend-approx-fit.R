# Main approximate extended CA-PPMx fitting function.
# Depends on Rcpp functions from Robust_stage_1and2.cpp.

############# extract relevant MCMC storage ####################################
extract_MCMC_storage <- function(input_df = NULL,
                                 input_specs = NULL,
                                 ref_trt = NULL,
                                 nmix = 15,
                                 del_range_response_1 = c(0.005, 0.02) * 15,  # for t=0
                                 nleapfrog_response_1 = 3,
                                 del_range_response_2 = c(0.005, 0.02) * 15,  # for t>0
                                 nleapfrog_response_2 = 3,
                                 del_range_alp1 = c(0.1, 0.3) * 5,
                                 nleapfrog_alp1 = 4,
                                 del_range_alp2 = c(0.1, 0.3) * 3,
                                 nleapfrog_alp2 = 3,
                                 a0 = 10.1,   # clustering hyper
                                 cc = 2       # prior mean scale for sigma^2
) {
  nrun <- 59000; burn <- 1000; thin <- 5
  
  # ---- basic checks ----------------------------------------------------------
  if (is.null(ref_trt)) stop("Please specify a reference treatment level.")
  if (!is.data.frame(input_df)) {
    stop("Please input a dataframe with at least: response, censoring indicator, ",
         "categorical and/or continuous covariates, treatment, and dataset columns.")
  }
  if (ncol(input_df) < 5) {
    stop("Input dataframe must contain >= 5 columns: response, censor_ind, ",
         "covariates, trt_type, dat_type.")
  }
  if (!is.list(input_specs)) {
    stop("Please pass a list 'input_specs' like:\n",
         "list(response='os', censor_ind='os_status', trt_type='trt', dat_type='type', ",
         "cat_vars=c('sex','kps','eor')) or with cont_vars, or both.")
  }
  nm_ok <- function(nmset)
    all(names(input_specs) %in% nmset) && all(nmset %in% names(input_specs))
  if (! ( nm_ok(c("response","censor_ind","trt_type","dat_type","cat_vars")) ||
          nm_ok(c("response","censor_ind","trt_type","dat_type","cont_vars")) ||
          nm_ok(c("response","censor_ind","trt_type","dat_type","cat_vars","cont_vars")) )) {
    stop("input_specs names must be one of:\n",
         "  response, censor_ind, trt_type, dat_type, cat_vars\n",
         "  response, censor_ind, trt_type, dat_type, cont_vars\n",
         "  response, censor_ind, trt_type, dat_type, cat_vars, cont_vars")
  }
  
  # existence checks
  if (!(ref_trt %in% unique(input_df[[input_specs$trt_type]])))
    stop("Reference treatment not found in input_df.")
  if (length(unique(input_df[[input_specs$dat_type]])) < 1L)
    stop("Need at least one dataset level in input_df.")
  
  # covariate presence
  isnullcatinput  <- is.null(input_specs$cat_vars)
  isnullcontinput <- is.null(input_specs$cont_vars)
  if (isnullcatinput && isnullcontinput)
    stop("You must have at least one categorical or continuous covariate.")
  
  # survival indicator must be logical
  surv_vec <- input_df[[input_specs$censor_ind]]
  if (!is.logical(surv_vec)) {
    stop("Survival indicators must be logical (TRUE=observed event, FALSE=censored).")
  }
  
  # ---- build indices (0-based) -----------------------------------------------
  # datasets: alphabetical order
  dat_levels <- sort(unique(input_df[[input_specs$dat_type]]))
  input_df$dat.index <- as.integer(factor(input_df[[input_specs$dat_type]],
                                          levels = dat_levels)) - 1L
  
  # treatments: ref_trt first, then others alphabetically
  other_trt <- setdiff(sort(unique(input_df[[input_specs$trt_type]])), ref_trt)
  trt_levels <- c(ref_trt, other_trt)
  input_df$trt.index <- as.integer(factor(input_df[[input_specs$trt_type]],
                                          levels = trt_levels)) - 1L
  
  # ---- prepare covariates for sampler (original encodings) -------------------
  # --- Save categorical level information (for Stage 2 consistency) ---
  if (!isnullcatinput) {
    cat_levels_list <- lapply(input_specs$cat_vars, function(v) {
      x <- input_df[[v]]
      if (is.factor(x)) {
        levels(x)                 # preserve original factor order
      } else {
        # coerce to character, then preserve the *observed* order of levels
        # (you can swap to sort(unique(.)) if you prefer alphabetical)
        unique(as.character(x))
      }
    })
    names(cat_levels_list) <- input_specs$cat_vars
  } else {
    cat_levels_list <- list()
  }
  
  # --- Build covariate matrices for the sampler -----------------------
  # CATEGORICAL: 0-based numeric codes; NA stays NA_real_
  if (isnullcatinput) {
    eta.cat <- matrix(NA_real_, nrow(input_df), 1L)
  } else {
    eta.cat <- as.matrix(
      data.frame(
        lapply(input_specs$cat_vars, function(v) {
          # force the exact Stage-1 level map, then convert to 0-based codes
          f <- factor(input_df[[v]], levels = cat_levels_list[[v]])
          as.numeric(f) - 1
        }),
        check.names = FALSE
      )
    )
  }
  
  # CONTINUOUS: original (non-scaled) values
  if (isnullcontinput) {
    eta.cont <- matrix(NA_real_, nrow(input_df), 1L)
  } else {
    eta.cont <- as.matrix(input_df[, input_specs$cont_vars, drop = FALSE])
  }
  
  rownames(eta.cat)  <- NULL
  rownames(eta.cont) <- NULL
  
  # ---- INITIAL CLUSTERING: Gower + PAM on ALL rows ---------------------------
  # Build a Gower-friendly data.frame: factors for cats, scaled numerics for conts
  build_gower_df <- function(df, cat_vars, cont_vars) {
    parts <- list()
    if (!is.null(cat_vars)) {
      cat_df <- df[ , cat_vars, drop = FALSE]
      cat_df[] <- lapply(cat_df, function(x) factor(x, exclude = NULL))
      parts <- c(parts, cat_df)
    }
    if (!is.null(cont_vars)) {
      cont_df <- df[ , cont_vars, drop = FALSE]
      # scale numerics ONLY for distance; sampler still gets unscaled eta.cont
      cont_df[] <- lapply(cont_df, function(x) as.numeric(scale(x)))
      parts <- c(parts, cont_df)
    }
    as.data.frame(parts)
  }
  
  gower_df <- build_gower_df(input_df, input_specs$cat_vars, input_specs$cont_vars)
  if (nrow(gower_df) == 1L) {
    labels <- 0L
  } else {
    ds <- cluster::daisy(gower_df, metric = "gower")
    tmp <- min(nmix, 11)
    k_use <- max(1L, min(tmp, nrow(gower_df)))
    pam_fit <- cluster::pam(ds, k = k_use, diss = TRUE)
    labels <- as.integer(pam_fit$clustering) - 1L   # 0-based for C++
  }
  
  # ---- other inputs to the sampler -------------------------------------------
  # response: log time; nu: logical event indicator (TRUE=event)
  response <- log(input_df[[input_specs$response]])
  surv_ind <- input_df[[input_specs$censor_ind]]
  
  # ncats: number of levels per categorical covariate (exclude NAs)
  if (!isnullcatinput) {
    ncats <- apply(eta.cat, 2, function(x) length(setdiff(unique(x), NA)))
  } else {
    ncats <- 1L
  }
  
  # non-missing indices per row (0-based)
  non.na.inds <- apply(eta.cat, 1, function(row) which(!is.na(row)) - 1L, simplify = TRUE)
  if (!is.list(non.na.inds)) non.na.inds <- split(non.na.inds, seq_len(nrow(eta.cat)))
  non.na.inds_cont <- apply(eta.cont, 1, function(row) which(!is.na(row)) - 1L, simplify = TRUE)
  if (!is.list(non.na.inds_cont)) non.na.inds_cont <- split(non.na.inds_cont, seq_len(nrow(eta.cont)))
  
  dat.index <- input_df$dat.index
  trt.index <- input_df$trt.index
  
  ## assume you have: response (log-time), surv_ind (TRUE=event), trt.index (0-based) in Stage 1
  df0 <- 1
  m   <- cc * (a0 - 1)        # target E[β]
  s2  <- 5                    # target Var[β] on β
  
  b_v <- log1p(s2 / m^2)      # Var(log β)
  b_m <- log(m) - b_v/2       # E(log β)
  
  trt_levels_stage1 <- sort(unique(trt.index))
  n_trt <- length(trt_levels_stage1)
  
  ## per-treatment μ priors from uncensored observations
  mu_m_t <- numeric(n_trt)
  mu_v_t <- numeric(n_trt)
  
  for (tt in seq_len(n_trt)) {
    tval <- trt_levels_stage1[tt]
    idx  <- which(trt.index == tval & surv_ind == TRUE)  # TRUE/1 for failures
    # guard: handle empty/low sample arms
    if (length(idx) >= 2L) {
      mu_m_t[tt] <- mean(response[idx], na.rm = TRUE)
      mu_v_t[tt] <- stats::var(response[idx], na.rm = TRUE)
    } else if (length(idx) == 1L) {
      mu_m_t[tt] <- response[idx[1L]]
      mu_v_t[tt] <- 1e-4
    } else {
      # fallback to pooled across all failures
      pool <- which(isTRUE(surv_ind))
      mu_m_t[tt] <- mean(response[pool], na.rm = TRUE)
      mu_v_t[tt] <- stats::var(response[pool], na.rm = TRUE)
    }
    if (!is.finite(mu_v_t[tt]) || mu_v_t[tt] <= 0) mu_v_t[tt] <- 1e-4
  }
  
  ## replicate β priors across treatments to preserve your original logic
  b_m_t <- rep(b_m, n_trt)
  b_v_t <- rep(b_v, n_trt)
  
  ## pass these 4 vectors to C++ instead of scalars:
  ##   mu_m_t, mu_v_t, b_m_t, b_v_t   (each length = n_trt, in the same treatment order)
  
  
  # ---- scale continuous variable ---------------------------------------------
  # capture means and sds for reuse in Stage 2
  cont_scale_params <- list(
    means = attr(scale(eta.cont), "scaled:center"),
    sds   = attr(scale(eta.cont), "scaled:scale")
  )
  
  eta.cont <- as.matrix(scale(eta.cont))  # z-scale all continuous columns
  
  # ---- call the Rcpp sampler -------------------------------------------------
  CAPPMx.result <- background_MCMC_storage(
    dat_index = dat.index,
    trt_index = trt.index,
    st = response, nu = surv_ind,
    del = labels, eta = eta.cat,
    eta_cont = eta.cont,
    non_na_obs1 = non.na.inds,
    non_na_obs1_cont = non.na.inds_cont,
    nmix = nmix, ncat = ncats,
    a0 = a0, df0 = df0,
    mu_m_t = mu_m_t, b_v_t = b_v_t,
    b_m_t = b_m_t, mu_v_t = mu_v_t,
    del_range_lognorm_ref = del_range_response_1,
    nleapfrog_lognorm_ref = nleapfrog_response_1,
    del_range_lognorm_oth = del_range_response_2,
    nleapfrog_lognorm_oth = nleapfrog_response_2,
    alpha_hyper = c(1, 10),
    del_range_alp1 = del_range_alp1,
    nleapfrog_alp1 = nleapfrog_alp1,
    del_range_alp2 = del_range_alp2,
    nleapfrog_alp2 = nleapfrog_alp2,
    nrun = nrun, burn = burn, thin = thin
  )
  
  # ---- augment return (same fields as before) --------------------------------
  CAPPMx.result$Treatment_Levels  <- trt_levels
  CAPPMx.result$Treatment_Indices <- seq_along(trt_levels)
  
  CAPPMx.result$Data_Levels  <- dat_levels
  CAPPMx.result$Data_Indices <- seq_along(dat_levels)
  
  dat_trt.index <- cbind(input_df$dat.index, input_df$trt.index)
  colnames(dat_trt.index) <- c("Data_Indices", "Treatment_Indices")
  CAPPMx.result$Combined_Indices <- unique(dat_trt.index) + 1L
  CAPPMx.result$raw_Indices      <- input_df[ , c("dat.index","trt.index")]
  
  CAPPMx.result$ncats <- ncats
  CAPPMx.result$Variable_specifications <- input_specs
  # Store inside Stage 1 result (important for keep clustering consistent)
  
  CAPPMx.result$Cat_Levels <- cat_levels_list
  CAPPMx.result$Cont_Scale_Params <- cont_scale_params
  return(CAPPMx.result)
}


##################### backwards compatibility ##################################

cappmx_fit=function(cat_cov_trt=NULL,cont_cov_trt=NULL, response_trt, surv_ind_trt,
                    cat_cov_rwd=NULL,cont_cov_rwd=NULL, response_rwd, surv_ind_rwd,
                    nmix=15, nrun=5e3,burn=1e3,thin=5,
                    del_range_response=c(.005,.02)*15, nleapfrog_response=3,
                    del_range_alp1 = c(.1,.3)*5, nleapfrog_alp1 = 4,
                    del_range_alp2 = c(.1,.3)*3, nleapfrog_alp2 = 3){
  
  isnullcat1=is.null(cat_cov_trt); isnullcont1=is.null(cont_cov_trt); isnullcat2=is.null(cat_cov_rwd); isnullcont2=is.null(cont_cov_rwd)
  
  if(isnullcat1 & isnullcont1 & isnullcat2 & isnullcont2 )
    stop("All covariate matrices are null!")
  
  if(class(surv_ind_trt)!="logical" | class(surv_ind_rwd)!="logical" )
    stop("Survival indicators must be logical variables!")
  # if(ncol(cat_cov_trt)!=ncol(cat_cov_rwd)| ncol(cont_cov_trt)!=ncol(cont_cov_rwd)| 
  #    nrow(cat_cov_trt)!=nrow(cont_cov_trt) | )
  
  ######CHECKS FOR NULL COV MATRICES AND SETTING NULLS
  ######IN TRT ARM
  if(isnullcat1 ){
    if(isnullcont1 )
      stop("Both covariate matrices are null in the treatment arm!") else{
        cat_cov_trt=matrix(NA,nrow=nrow(cont_cov_trt),ncol=1)
      }
  }
  if(isnullcont1 ){
    if(isnullcat1)
      stop("Both covariate matrices are null in the treatment arm!") else{
        cont_cov_trt=matrix(NA,nrow=nrow(cat_cov_trt),ncol=1)
      }
  }
  ########################
  
  ######IN RWD######
  if(isnullcat2 ){
    if(isnullcont2 )
      stop("Both covariate matrices are null in the RWD!") else{
        cat_cov_rwd=matrix(NA,nrow=nrow(cont_cov_rwd),ncol=1)
      }
  }
  if(isnullcont2 ){
    if(isnullcat2)
      stop("Both covariate matrices are null in the RWD!") else{
        cont_cov_rwd=matrix(NA,nrow=nrow(cat_cov_rwd),ncol=1)
      }
  }
  ########################################################################
  
  eta.cont1=as.matrix(cont_cov_trt)
  eta.cont2=as.matrix(cont_cov_rwd)
  eta.cont=rbind(eta.cont1,eta.cont2)
  
  eta.cat1=as.matrix(cat_cov_trt); rownames(eta.cat1)=NULL ;nsamp1=nrow(eta.cat1)
  eta.cat2=as.matrix(cat_cov_rwd); rownames(eta.cat2)=NULL ;nsamp2=nrow(eta.cat2)
  eta.cat=rbind(eta.cat1,eta.cat2)
  ncats=apply(eta.cat,2, function(x) length(setdiff(unique(x),NA)))
  
  ################################################################
  #########GET THE INDICES OF THE NON-MISSING VALUES##############
  ################################################################
  
  ###############
  ###for the categorical covs 
  ###############
  non.na.inds=apply(eta.cat,1, function(eta) purrr::compose( which,"!",is.na)(eta)-1,simplify = F)
  
  ###############
  ###for the continuous covs 
  ###############
  non.na.inds_cont=apply(eta.cont,1, function(eta) purrr::compose( which,"!",is.na)(eta)-1,simplify = F)
  
  ##################################################################
  
  ###############Initial values for MCMC input###########
  ## For computing distance matrix
  
  hamm_dist=function(x,y){
    nonnax=which(!is.na(x))
    nonnay=which(!is.na(y))
    int=intersect(nonnax,nonnay)
    if(length(int)==0)
      return (length(x)) else return(e1071::hamming.distance(x[int],y[int]))
  }
  
  ds=matrix(0,nsamp2,nsamp2)
  if(!isnullcat2 ){
    for(i in 1:nsamp2){
      for(j in (i):nsamp2){
        if(j>nsamp2)
          print(j)
        ds[i,j]=hamm_dist(eta.cat2[i,],eta.cat2[j,])
        ds[j,i]=ds[i,j]
      }
    }
  }
  ds=as.dist(ds)
  
  if(!isnullcont2 ){
    x=eta.cont2; x[!is.finite(x)]=0
    ds=ds+ dist(x)
    rm(x)
  }
  
  fit <- hclust (ds , method = "ward.D2" )
  labels2 <- cutree ( fit , k =min(nmix,7)) -1
  
  #############NULL CHECKS###############
  if(isnullcat2 ){
    eta.imputed2=randomForest::na.roughfix(eta.cont2)
    eta.imputed1=randomForest::na.roughfix(eta.cont1)
  } else if(isnullcont2 ){
    eta.imputed2=randomForest::na.roughfix(eta.cat2)
    eta.imputed1=randomForest::na.roughfix(eta.cat1)
  } else if(!isnullcat1 & !isnullcont1 & !isnullcat2 & !isnullcont2 ){
    eta.imputed2=cbind(randomForest::na.roughfix(eta.cat2),randomForest::na.roughfix(eta.cont2))
    eta.imputed1=cbind(randomForest::na.roughfix(eta.cat1),randomForest::na.roughfix(eta.cont1))
  }  else stop("Check input covariates!")
  ############################################################
  
  rffit=randomForest::randomForest(x=eta.imputed2,y=as.factor(labels2), na.action=na.roughfix)
  labels1=as.numeric(predict(rffit,eta.imputed1))-1
  
  # cat_cov_trt,cont_cov_trt, response_trt, surv_ind_trt,
  # cat_cov_rwd,cont_cov_rwd, response_rwd, surv_ind_rwd
  
  # ind1=which(surv_ind_trt); ind2=which(surv_ind_rwd)
  log_observed_failure_times=(c( response_trt[surv_ind_trt], response_rwd[surv_ind_rwd]))
  
  ####log-normal hyperparameter tuning
  df0=1
  a0=10.1 ##must be greater than 1
  m= 2*(a0-1) #var(log_observed_failure_times)*(a0-1) ###mean of variance of outcomes
  s=5
  b_v=log1p(s/m^2)
  b_m=log(m)- b_v/2
  mu_m= mean(log_observed_failure_times)
  mu_v= 1
  
  common_atoms_cat_lognormal( nmix, ncat=ncats,
                              a0=a0, df0=df0, mu_m=mu_m, mu_v=mu_v, b_m=b_m, b_v=b_v,#normal and lognormal hyperparameters for the response
                              nrun=nrun, burn=burn, thin=thin,
                              eta.cat1  , eta.cat2,
                              eta.cont1, eta.cont2,
                              response_trt, surv_ind_trt,
                              response_rwd, surv_ind_rwd,
                              non.na.inds, non.na.inds_cont,
                              labels1, labels2, 
                              del_range_lognorm=del_range_response, nleapfrog_lognorm=nleapfrog_response,
                              alpha_hyper = c(1,10),del_range_alp1 = del_range_alp1, nleapfrog_alp1 = nleapfrog_alp1,
                              del_range_alp2 = del_range_alp2, nleapfrog_alp2 = nleapfrog_alp2)
}

################################################################################



cappmx_extend_approx_fit <- function(result_CAPPMx, input_df = NULL,
                                     input_specs = NULL, ref_trt = NULL,
                                     input_df_pred = NULL,
                                     del_range_response_1 = c(0.005, 0.02) * 15,  # for t=0
                                     nleapfrog_response_1 = 4,
                                     del_range_response_2 = c(0.005, 0.02) * 13,  # for t>0
                                     nleapfrog_response_2 = 3,
                                     del_range_alp1 = c(0.1, 0.3) * 5,
                                     nleapfrog_alp1 = 4,
                                     del_range_alp2 = c(0.1, 0.3) * 3,
                                     nleapfrog_alp2 = 3,
                                     a0 = 10.1,  # clustering hyper
                                     cc = 2, # prior mean scale for sigma^2
                                     freeze_control = TRUE,
                                     eb_beta_unique = TRUE,         # NEW: if TRUE, set arm-specific beta0,t for Stage-2–unique arms
                                     beta_temper_tau     = 10,      # EB weight w = n_fail / (n_fail + tau)
                                     beta_cap_q          = c(0.10, 0.90),  # clamp β to control-based quantiles
                                     beta_var_floor_mult = 1.0,     # inflate floor on log-β variance if needed
                                     min_failures_for_EB = 3)       # below this, w := 0 (fall back to pool) 
{
  # normalize inputs (avoid tibble one-index surprises)
  if (!is.null(input_df))      input_df      <- as.data.frame(input_df)
  if (!is.null(input_df_pred)) input_df_pred <- as.data.frame(input_df_pred)
  `%||%` <- function(x, y) if (is.null(x)) y else x
  
  # ============= FIT BRANCH =============
  if (!is.data.frame(input_df_pred)) {
    
    if (length(result_CAPPMx) != 30)
      stop("You might be missing some essential output from Stage 1.")
    nmix <- nrow(result_CAPPMx$nj_val_all_mat)
    
    if (is.null(ref_trt)) stop("Please specify a reference treatment level.")
    if (!is.data.frame(input_df))
      stop("input_df must be a dataframe with response, censor, trt, dat, and covariates columns.")
    # if (!is.data.frame(input_df_pred))
    #   stop("input_df_pred must be a dataframe containing covariate columns for prediction.")
    
    stage1_cat_levels <- result_CAPPMx$Cat_Levels
    isnullcatinput  <- is.null(input_specs$cat_vars)
    isnullcontinput <- is.null(input_specs$cont_vars)
    if (isnullcatinput & isnullcontinput)
      stop("You must have at least one categorical or continuous covariate.")
    
    stage1_levels <- result_CAPPMx$Data_Levels
    in_df_levels  <- unique(as.vector(input_df[[input_specs$dat_type]]))
    curr_candidates <- setdiff(in_df_levels, stage1_levels)
    if (length(curr_candidates) != 1) {
      stop("Could not infer a unique current cohort from '", input_specs$dat_type,
           "'. Exactly one level in input_df must be new (not in Stage 1).")
    }
    curr_dat_name <- curr_candidates
    rwd_dat_names <- stage1_levels
    
    is_rwd <- input_df[[input_specs$dat_type]] %in% rwd_dat_names
    
    input_df$dat.index <- as.numeric(factor(
      as.vector(input_df[[input_specs$dat_type]]),
      levels = c(curr_dat_name, rwd_dat_names),
      labels = 1:(1 + length(rwd_dat_names))
    )) - 1L
    
    ref_trt_name    <- ref_trt
    other_trt_names <- setdiff(
      unique(c(as.vector(input_df[[input_specs$trt_type]]),
               result_CAPPMx$Treatment_Levels)),
      ref_trt_name
    )
    trt_levels <- c(ref_trt_name, other_trt_names)
    dat_levels <- c(curr_dat_name, rwd_dat_names)
    
    input_df$trt.index <- as.numeric(factor(
      as.vector(input_df[[input_specs$trt_type]]),
      levels = trt_levels,
      labels = seq_along(trt_levels)
    )) - 1L
    
    dat.trt <- input_df[!is_rwd, , drop = FALSE]
    dat.rwd <- input_df[ is_rwd, , drop = FALSE]
    
    if (nrow(dat.trt) == 0L) stop("No rows found for the current cohort in input_df.")
    
    
    # --- Build trt.convert (NEW before initializer!) ---
    trt.convert <- cbind(
      result_CAPPMx$raw_Indices$trt.index,  # OLD trt index (Stage 1 space, 0-based)
      as.numeric(factor(
        factor(result_CAPPMx$raw_Indices$trt.index,
               levels = result_CAPPMx$Treatment_Indices - 1L,
               labels = result_CAPPMx$Treatment_Levels),
        levels = trt_levels,
        labels = seq_along(trt_levels)
      )) - 1L                                # NEW trt index (Stage 2 space, 0-based)
    )
    colnames(trt.convert) <- c("old_trt_index","New_trt_index")
    trt.convert <- unique(trt.convert)
    storage.mode(trt.convert) <- "integer"
    
    ## ---- rebuild covariates in Stage-1’s encoding/scaling ----
    spec_stage1 <- result_CAPPMx$Variable_specifications
    # stage1_cat_levels <- result_CAPPMx$Cat_Levels
    stage1_cats  <- spec_stage1$cat_vars  %||% character(0)
    stage1_conts <- spec_stage1$cont_vars %||% character(0)
    
    ## --- Categorical covariates (0-based, in Stage-1 order) ---
    build_cat_block <- function(df_rows) {
      if (length(stage1_cats) == 0L) {
        mat <- matrix(NA_real_, nrow(df_rows), 1L)
        colnames(mat) <- NULL
        return(mat)
      }
      cols <- lapply(stage1_cats, function(v) {
        if (!v %in% names(stage1_cat_levels)) {
          stop("Stage 1 Cat_Levels missing variable '", v, "'.")
        }
        if (v %in% names(df_rows)) {
          f <- factor(df_rows[[v]], levels = stage1_cat_levels[[v]])
          ## If new data has unseen levels (non-NA in data but NA after factoring), error:
          if (any(is.na(f) & !is.na(df_rows[[v]]))) {
            unseen <- unique(df_rows[[v]][is.na(f) & !is.na(df_rows[[v]])])
            stop("Unseen level(s) for '", v, "': ", paste(unseen, collapse = ", "))
          }
          as.numeric(f) - 1L
        } else {
          ## variable absent in new data -> keep alignment: all NA
          rep(NA_real_, nrow(df_rows))
        }
      })
      mat <- as.matrix(data.frame(cols, check.names = FALSE))
      colnames(mat) <- stage1_cats
      mat
    }
    
    cat_cov_trt <- build_cat_block(dat.trt)
    cat_cov_rwd <- build_cat_block(dat.rwd)
    
    ## --- Continuous covariates (raw scale first, Stage-1 order) ---
    build_cont_block <- function(df_rows) {
      if (length(stage1_conts) == 0L) {
        mat <- matrix(NA_real_, nrow(df_rows), 1L)
        colnames(mat) <- NULL
        return(mat)
      }
      cols <- lapply(stage1_conts, function(v) {
        if (v %in% names(df_rows)) as.numeric(df_rows[[v]]) else rep(NA_real_, nrow(df_rows))
      })
      mat <- as.matrix(data.frame(cols, check.names = FALSE))
      colnames(mat) <- stage1_conts
      mat
    }
    
    cont_cov_trt <- build_cont_block(dat.trt)
    cont_cov_rwd <- build_cont_block(dat.rwd)
    
    ## --- Stack in (current; historical) order used below ---
    eta.cat  <- rbind(cat_cov_trt,  cat_cov_rwd)
    eta.cont <- rbind(cont_cov_trt, cont_cov_rwd)
    rownames(eta.cat) <- rownames(eta.cont) <- NULL
    
    ## --- Non-missing indices per row (0-based) ---
    non.na.inds      <- apply(eta.cat,  1, function(x) which(!is.na(x)) - 1L, simplify = FALSE)
    non.na.inds_cont <- apply(eta.cont, 1, function(x) which(!is.na(x)) - 1L, simplify = FALSE)
    
    ## --- Scale continuous variables using Stage 1 summaries (aligned by name) ---
    if (length(stage1_conts) > 0L) {
      means_stage1 <- result_CAPPMx$Cont_Scale_Params$means[stage1_conts]
      sds_stage1   <- result_CAPPMx$Cont_Scale_Params$sds[stage1_conts]
      ## guard: replace 0/NA sds with 1 to avoid blow-ups
      sds_stage1[!is.finite(sds_stage1) | sds_stage1 == 0] <- 1
      
      ## sweep expects numeric; keep NAs as NAs
      eta.cont <- sweep(eta.cont, 2, means_stage1, "-")
      eta.cont <- sweep(eta.cont, 2, sds_stage1,   "/")
    } else {
      ## keep single NA column unchanged
    }
    
    ## --- ncats comes from Stage 1 and already matches Stage-1 cat order ---
    ncats <- result_CAPPMx$ncats
    
    response <- log(c(dat.trt[[input_specs$response]], dat.rwd[[input_specs$response]]))
    surv_ind <-      c(dat.trt[[input_specs$censor_ind]], dat.rwd[[input_specs$censor_ind]])
    
    map_t_old <- function(t_new) {
      hit <- which(trt.convert[, "New_trt_index"] == t_new)
      if (length(hit) == 0L) return(-1L)
      as.integer(trt.convert[hit[1L], "old_trt_index"])
    }
    
    # --- rows that belong to "current" cohort (to initialise) ---
    n_trt_rows <- nrow(dat.trt)
    n_rwd_rows <- nrow(dat.rwd)
    n_all      <- n_trt_rows + n_rwd_rows
    
    ## ===== Init via first half draws; run Stage 2 on last half draws =====
    suppressWarnings(suppressMessages({
      library(matrixStats)
      library(mcclust)   # for minbinder
    }))
    
    ## MCMC dims (we'll still use M and K from Stage 1)
    M  <- dim(result_CAPPMx$Lognormal_Mu_Cube)[1]
    K  <- dim(result_CAPPMx$Lognormal_Mu_Cube)[2]
    stopifnot(K == nrow(result_CAPPMx$nj_val_all_mat))
    
    ## choose the first half for initialisation, last half for the sampler later
    init_ids <- 0:((M/2) - 1L)               # 0-based
    # run_ids  <- (M - (M/2)):(M - 1L)         # 0-based
    
    ## dims
    N <- nrow(dat.trt)
    R <- length(init_ids)
    k_cat  <- ncol(eta.cat)
    k_cont <- ncol(eta.cont)
    
    if (length(stage1_cats) == 0L) k_cat  <- 0L
    if (length(stage1_conts) == 0L) k_cont <- 0L
    
    ## two cubes: [R x K x N]
    cov_loglik_cat  <- array(0, dim = c(R, K, N))
    cov_loglik_cont <- array(0, dim = c(R, K, N))
    
    ## Stage-1 hyper for continuous predictive
    df_x    <- 1.0
    alpha_x <- k_cont + 30.0
    beta_x  <- 1.0
    mu_x    <- 0.0
    
    ## Pull Stage-1 covariate tallies
    nobs_cube  <- result_CAPPMx$nobs_cube      # [K x k_cat  x M]
    nj_x_cube  <- result_CAPPMx$nj_x_cube      # [K x k_cont x M]
    sum_x_cube <- result_CAPPMx$sum_j_x_cube   # [K x k_cont x M]
    ss_x_cube  <- result_CAPPMx$ss_j_x_cube    # [K x k_cont x M]
    noccu_list <- result_CAPPMx$noccu_list     # length M; [[m]][[j]][[v]] vector over levels
    
    ## fill cubes
    cov_out <- stage2_cov_loglik_precompute(
      init_ids        = init_ids,
      eta_cat         = eta.cat,
      eta_cont        = eta.cont,
      non_na_obs1     = non.na.inds,
      non_na_obs1_cont= non.na.inds_cont,
      nobs_cube       = result_CAPPMx$nobs_cube,
      nj_x_cube       = result_CAPPMx$nj_x_cube,
      sum_x_cube      = result_CAPPMx$sum_j_x_cube,
      ss_x_cube       = result_CAPPMx$ss_j_x_cube,
      noccu_old       = result_CAPPMx$noccu_list,
      ncat            = ncats,
      df_x            = df_x,
      alpha_x         = alpha_x,
      mu_x            = mu_x,
      beta_x          = beta_x
    )
    
    cov_loglik_cat  <- cov_out$cov_loglik_cat   # reorder to [R,K,N]
    cov_loglik_cont <- cov_out$cov_loglik_cont
    
    ## If there are no cats/conts, the arrays above are already zeros with correct dims.
    
    ## ---- call the covariates-only initializer ----
    Z_init <- stage2_init_alloc_cov_only(
      cov_loglik_cat  = cov_loglik_cat,
      cov_loglik_cont = cov_loglik_cont,
      init_ids        = init_ids,
      nmix            = K,
      rng_seed        = 500L
    )
    ## Z_init is [R x N] with 0-based labels
    
    
    # --- consensus init via PSM + MinBinder ---
    # Z_init is [R x n_trt_rows] with 0-based labels; convert to 1-based for mcclust
    Z1 <- Z_init + 1L
    # co-clustering (PSM)
    psm <- mcclust::comp.psm(Z1)
    init_labels_trt_1b <- mcclust::minbinder(psm, max.k = min(15, K))$cl
    init_labels_trt    <- init_labels_trt_1b - 1L  # back to 0-based
    
    # assemble initial labels for full data (current + historical)
    labels1 <- init_labels_trt
    labels2 <- integer(n_rwd_rows)  # historical can be ignored for init; sampler will use Stage-1
    labels  <- c(labels1, labels2)
    
    # --- Indices for sampler ---
    dat.index <- c(dat.trt$dat.index, result_CAPPMx$raw_Indices$dat.index + 1L)
    trt.index <- c(
      dat.trt$trt.index,
      as.numeric(factor(
        factor(result_CAPPMx$raw_Indices$trt.index,
               levels = result_CAPPMx$Treatment_Indices - 1L,
               labels = result_CAPPMx$Treatment_Levels),
        levels = trt_levels,
        labels = seq_along(trt_levels)
      )) - 1L
    )
    
    
    # right before the C++ call
    storage.mode(eta.cat)  <- "integer"
    storage.mode(eta.cont) <- "double"
    storage.mode(dat.index) <- "integer"
    storage.mode(trt.index) <- "integer"
    
    # --- Run sampler ---
    burn <- 0; thin <- 5; nrun <- (nrow(result_CAPPMx$picube) - burn)/2
    
    ## discarding half
    new_M <- (M - (M/2) + 1):(M)
    
    ## --- map Stage-2 control (new index 0) back to Stage-1 old index ---
    ctl_row <- which(trt.convert[, "New_trt_index"] == 0L)
    if (length(ctl_row) != 1L) stop("Could not identify the control row in trt.convert.")
    old_ctl_0b <- trt.convert[ctl_row, "old_trt_index"]     # 0-based
    old_ctl_1b <- old_ctl_0b + 1L                           # 1-based for array indexing
    
    ## --- pull per-draw μ0,0 and β0,0 from Stage 1 (last-half only) ---
    H <- result_CAPPMx$Lognormal_hyperparams  # [M x 2 x T_old]; [,1]=mu0_t, [,2]=beta0_t
    mu0_control_draws   <- as.numeric(H[new_M, 1, old_ctl_1b])   # length = length(new_M)
    beta0_control_draws <- as.numeric(H[new_M, 2, old_ctl_1b])   # positive, natural scale
    
    ## small guardrails
    if (any(!is.finite(mu0_control_draws)))   stop("Non-finite control μ0 draws.")
    if (any(!is.finite(beta0_control_draws) | beta0_control_draws <= 0))
      stop("Non-finite or non-positive control β0 draws.")
    
    # ------------------------------------------------------------------------------
    # compute_stage2_hyperpriors_t()
    #
    # Purpose:
    #   Construct treatment-specific lognormal hyperparameters (μ₀,t, β₀,t)
    #   for Stage 2 of the approximate model.
    #
    # Behavior:
    #   • For overlapping (shared) treatment arms — i.e., arms present in both
    #     Stage 1 and Stage 2 — this function borrows empirical hyperparameters
    #     from the Stage 1 posterior (computed from the Lognormal_hyperparams cube).
    #
    #   • For Stage-2–only (unique) treatment arms — e.g., new experimental drugs
    #     not seen in Stage 1 — it rebuilds the hyperparameters from the current
    #     Stage-2 data using the same recipe as Stage 1:
    #         - μ₀,t hyperparameters (mean, variance) estimated from log-times
    #           among uncensored failures in that arm.
    #         - β₀,t hyperparameters constructed from target moments:
    #             E[β]  = cc × (a₀ − 1)
    #             Var[β] = s²_β = 5
    #           converted to log-scale mean/variance (b_m_t, b_v_t).
    #
    #   • Final values are sanitized to ensure positivity and finiteness.
    #     Returns a list with (μ_m_t, μ_v_t, b_m_t, b_v_t) and diagnostics.
    #
    # Inputs:
    #   - result_CAPPMx : Stage 1 fitted object containing Lognormal_hyperparams
    #   - trt.convert   : 0-based mapping of Stage-1 to Stage-2 treatment indices
    #   - T_new         : total number of treatments in Stage 2
    #   - input_df      : Stage-2 data with columns trt, os, os_status
    #   - ref_trt       : reference treatment label (usually "Control")
    #   - a0, cc, s2_beta : parameters for the Stage-1 prior recipe (used only for
    #                       Stage-2-unique arms)
    #
    # Output:
    #   A list with:
    #       mu_m_t, mu_v_t, b_m_t, b_v_t,
    #       diagnostics (overlap map, pooled stats, arm roles)
    #
    # ------------------------------------------------------------------------------
    compute_stage2_hyperpriors_t <- function(
    result_CAPPMx,
    trt.convert,               # 0-based [old, new] mapping
    T_new,                     # number of Stage-2 treatments
    input_df,                  # Stage-2 dataset
    input_specs,               # app-provided column names
    ref_trt = "Control",
    use_last_half = TRUE,
    a0 = 10.1,                 # used only for Stage-2-unique arms
    cc = 2,                    # used only for Stage-2-unique arms
    s2_beta = 5,               # target Var[beta] used in Stage-1 recipe
    eb_beta_unique = FALSE,
    beta_temper_tau     = 10,
    beta_cap_q          = c(0.10, 0.90),
    beta_var_floor_mult = 1.0,
    min_failures_for_EB = 3
    ) {
      if (is.null(result_CAPPMx$Lognormal_hyperparams)) {
        stop("Stage 1 output missing $Lognormal_hyperparams.")
      }
      
      H <- result_CAPPMx$Lognormal_hyperparams  # [M x 2 x T_old]; [,1]=mu0_t, [,2]=beta0_t
      
      dH <- dim(H)
      if (length(dH) != 3L || dH[2] != 2L) {
        stop("$Lognormal_hyperparams must be [M x 2 x T_old].")
      }
      
      M <- dH[1]
      T_old <- dH[3]
      
      # Required column names from input_specs
      trt_col <- input_specs$trt_type
      os_col <- input_specs$response
      censor_col <- input_specs$censor_ind
      
      required_model_cols <- c(trt_col, os_col, censor_col)
      missing_model_cols <- setdiff(required_model_cols, names(input_df))
      
      if (length(missing_model_cols) > 0) {
        stop(
          "input_df is missing required model-fitting column(s): ",
          paste(missing_model_cols, collapse = ", ")
        )
      }
      
      # Draws to use
      draw_idx <- if (use_last_half) {
        seq.int(floor(M / 2) + 1L, M)
      } else {
        seq_len(M)
      }
      
      # Tidy treatment conversion table
      tc <- as.data.frame(trt.convert)
      names(tc) <- c("old", "new")
      tc$new[is.na(tc$new)] <- -1L
      tc <- unique(tc)
      
      # Validate T_new
      stopifnot(length(T_new) == 1L, is.finite(T_new), T_new >= 1L)
      
      # Pre-allocate in Stage-2 treatment order
      mu_m_t <- mu_v_t <- b_m_t <- b_v_t <- rep(NA_real_, T_new)
      
      # ---------- SHARED / OVERLAPPING ARMS ----------
      overlap <- tc[
        tc$new >= 0L & tc$new < T_new &
          tc$old >= 0L & tc$old < T_old,
        ,
        drop = FALSE
      ]
      
      if (nrow(overlap) == 0L) {
        stop("No overlapping treatments between Stage 1 and Stage 2.")
      }
      
      per_arm <- lapply(seq_len(nrow(overlap)), function(ii) {
        t_old0 <- overlap$old[ii]   # 0-based Stage-1 index
        t_new  <- overlap$new[ii]   # 0-based Stage-2 index
        t_old  <- t_old0 + 1L       # 1-based index for R array
        
        mu_draws   <- as.numeric(H[draw_idx, 1, t_old])
        beta_draws <- as.numeric(H[draw_idx, 2, t_old])
        
        ok_mu <- is.finite(mu_draws)
        ok_b  <- is.finite(beta_draws) & beta_draws > 0
        
        if (!any(ok_mu) || !any(ok_b)) {
          return(list(t_old = t_old0, t_new = t_new, ok = FALSE))
        }
        
        list(
          t_old     = t_old0,
          t_new     = t_new,
          ok        = TRUE,
          mu_mean   = mean(mu_draws[ok_mu]),
          mu_var    = stats::var(mu_draws[ok_mu]),
          logb_mean = mean(log(beta_draws[ok_b])),
          logb_var  = stats::var(log(beta_draws[ok_b]))
        )
      })
      
      keep <- vapply(per_arm, function(x) isTRUE(x$ok), logical(1))
      per_arm <- per_arm[keep]
      
      if (length(per_arm) == 0L) {
        stop("No valid posterior draws found for overlapping treatments.")
      }
      
      for (pa in per_arm) {
        idx <- pa$t_new + 1L
        
        mu_m_t[idx] <- pa$mu_mean
        mu_v_t[idx] <- if (is.finite(pa$mu_var) && pa$mu_var > 0) {
          pa$mu_var
        } else {
          NA_real_
        }
        
        b_m_t[idx] <- pa$logb_mean
        b_v_t[idx] <- if (is.finite(pa$logb_var) && pa$logb_var > 0) {
          pa$logb_var
        } else {
          NA_real_
        }
      }
      
      # Pooled fallbacks from overlapping arms
      mu_means   <- vapply(per_arm, `[[`, numeric(1), "mu_mean")
      mu_vars    <- vapply(per_arm, `[[`, numeric(1), "mu_var")
      logb_means <- vapply(per_arm, `[[`, numeric(1), "logb_mean")
      logb_vars  <- vapply(per_arm, `[[`, numeric(1), "logb_var")
      
      w_mu <- ifelse(is.finite(mu_vars) & mu_vars > 0, 1 / mu_vars, 0)
      w_logb <- ifelse(is.finite(logb_vars) & logb_vars > 0, 1 / logb_vars, 0)
      
      mu_pool <- if (sum(w_mu) > 0) {
        sum(w_mu * mu_means) / sum(w_mu)
      } else {
        mean(mu_means)
      }
      
      logb_pool <- if (sum(w_logb) > 0) {
        sum(w_logb * logb_means) / sum(w_logb)
      } else {
        mean(logb_means)
      }
      
      base_mu_v <- mean(
        ifelse(is.finite(mu_vars) & mu_vars > 0, mu_vars, NA),
        na.rm = TRUE
      )
      
      base_logb_v <- mean(
        ifelse(is.finite(logb_vars) & logb_vars > 0, logb_vars, NA),
        na.rm = TRUE
      )
      
      base_mu_v <- if (is.finite(base_mu_v) && base_mu_v > 0) {
        base_mu_v
      } else {
        1e-6
      }
      
      base_logb_v <- if (is.finite(base_logb_v) && base_logb_v > 0) {
        base_logb_v
      } else {
        1e-6
      }
      
      bad_mu_v <- !is.finite(mu_v_t) | mu_v_t <= 0
      bad_b_v  <- !is.finite(b_v_t)  | b_v_t  <= 0
      
      if (any(bad_mu_v)) {
        mu_v_t[bad_mu_v] <- base_mu_v
      }
      
      if (any(bad_b_v)) {
        b_v_t[bad_b_v] <- base_logb_v
      }
      
      # ---------- CONTROL-BASED CAPS FOR BETA ----------
      ctl_row <- which(tc$new == 0L)
      
      if (length(ctl_row) == 1L &&
          tc$old[ctl_row] >= 0L &&
          tc$old[ctl_row] < T_old) {
        
        t_old_ctl <- tc$old[ctl_row] + 1L
        ctl_logbeta <- log(as.numeric(H[draw_idx, 2, t_old_ctl]))
        ctl_logbeta <- ctl_logbeta[is.finite(ctl_logbeta)]
        
        if (length(ctl_logbeta) >= 20L) {
          beta_cap_lo <- as.numeric(
            stats::quantile(exp(ctl_logbeta), probs = beta_cap_q[1])
          )
          beta_cap_hi <- as.numeric(
            stats::quantile(exp(ctl_logbeta), probs = beta_cap_q[2])
          )
        } else {
          beta_cap_lo <- exp(logb_pool)
          beta_cap_hi <- exp(logb_pool)
        }
      } else {
        beta_cap_lo <- exp(logb_pool)
        beta_cap_hi <- exp(logb_pool)
      }
      
      beta_logvar_floor <- max(base_logb_v * beta_var_floor_mult, 1e-6)
      
      # ---------- UNIQUE / STAGE-2-ONLY ARMS ----------
      overlap_new <- vapply(per_arm, `[[`, integer(1), "t_new")
      stage2_only_0b <- setdiff(seq_len(T_new) - 1L, unique(overlap_new))
      
      # Treatment labels in the uploaded/current dataset
      arm_levels <- levels(factor(input_df[[trt_col]]))
      
      if (length(arm_levels) < T_new) {
        arm_levels <- unique(c(arm_levels, ref_trt))
      }
      
      if (length(stage2_only_0b)) {
        
        # Original Stage-1 beta recipe fallback
        if (!isTRUE(eb_beta_unique)) {
          m <- cc * (a0 - 1)
          b_v_uni <- log1p(s2_beta / m^2)
          b_m_uni <- log(m) - 0.5 * b_v_uni
        }
        
        for (t_new in stage2_only_0b) {
          idx <- t_new + 1L
          
          t_lab <- if (!is.na(arm_levels[idx])) {
            arm_levels[idx]
          } else {
            ref_trt
          }
          
          arm_df <- input_df[
            as.character(input_df[[trt_col]]) == as.character(t_lab) &
              input_df[[censor_col]] == 1,
            ,
            drop = FALSE
          ]
          
          # ----- mu hyperparameters from Stage-2 failures -----
          if (nrow(arm_df) >= 2L) {
            mu_m_t[idx] <- mean(log(arm_df[[os_col]]), na.rm = TRUE)
            mu_v_t[idx] <- stats::var(log(arm_df[[os_col]]), na.rm = TRUE)
          } else if (nrow(arm_df) == 1L) {
            mu_m_t[idx] <- log(arm_df[[os_col]][1L])
            mu_v_t[idx] <- 1e-4
          } else {
            pool <- input_df[
              input_df[[censor_col]] == 1,
              ,
              drop = FALSE
            ]
            
            mu_m_t[idx] <- mean(log(pool[[os_col]]), na.rm = TRUE)
            mu_v_t[idx] <- stats::var(log(pool[[os_col]]), na.rm = TRUE)
          }
          
          if (!is.finite(mu_v_t[idx]) || mu_v_t[idx] <= 0) {
            mu_v_t[idx] <- 1e-4
          }
          
          # ----- beta hyperparameters -----
          if (isTRUE(eb_beta_unique)) {
            nf <- nrow(arm_df)
            
            s2_arm <- if (nf >= 2L) {
              stats::var(log(arm_df[[os_col]]), na.rm = TRUE)
            } else {
              NA_real_
            }
            
            w <- if (is.finite(s2_arm) && nf >= min_failures_for_EB) {
              nf / (nf + beta_temper_tau)
            } else {
              0
            }
            
            logb_EB <- if (is.finite(s2_arm) && s2_arm > 0) {
              log(s2_arm)
            } else {
              logb_pool
            }
            
            b_m_t[idx] <- w * logb_EB + (1 - w) * logb_pool
            
            beta_center <- exp(b_m_t[idx])
            
            if (is.finite(beta_cap_lo) &&
                is.finite(beta_cap_hi) &&
                beta_cap_hi > beta_cap_lo) {
              beta_center <- min(max(beta_center, beta_cap_lo), beta_cap_hi)
            }
            
            b_m_t[idx] <- log(beta_center)
            b_v_t[idx] <- beta_logvar_floor
          } else {
            b_m_t[idx] <- b_m_uni
            b_v_t[idx] <- b_v_uni
          }
        }
      }
      
      # ---------- FINAL SANITATION ----------
      mu_m_t[!is.finite(mu_m_t)] <- mu_pool
      b_m_t[!is.finite(b_m_t)] <- logb_pool
      
      mu_v_t[!is.finite(mu_v_t) | mu_v_t <= 0] <- base_mu_v
      b_v_t[!is.finite(b_v_t) | b_v_t <= 0] <- base_logb_v
      
      stopifnot(
        length(mu_m_t) == T_new,
        length(mu_v_t) == T_new,
        length(b_m_t)  == T_new,
        length(b_v_t)  == T_new
      )
      
      list(
        mu_m_t = mu_m_t,
        mu_v_t = mu_v_t,
        b_m_t  = b_m_t,
        b_v_t  = b_v_t,
        diagnostics = list(
          overlap_map = data.frame(
            old = vapply(per_arm, `[[`, integer(1), "t_old"),
            new = vapply(per_arm, `[[`, integer(1), "t_new")
          ),
          pooled = c(
            mu_pool = mu_pool,
            logb_pool = logb_pool,
            base_mu_v = base_mu_v,
            base_logb_v = base_logb_v
          ),
          stage2_only = stage2_only_0b,
          arm_role = {
            role <- rep("shared", T_new)
            if (length(stage2_only_0b)) {
              role[stage2_only_0b + 1L] <- "unique"
            }
            role
          }
        )
      )
    }
    
    T_new <- length(trt_levels)  # or max(trt.index) + 1L
    hp_t <- compute_stage2_hyperpriors_t(
      result_CAPPMx = result_CAPPMx,
      trt.convert   = trt.convert,
      T_new         = T_new,
      input_df      = input_df,
      input_specs   = input_specs,
      ref_trt       = ref_trt,
      use_last_half = TRUE,
      a0 = a0,
      cc = cc,
      s2_beta = 5,
      eb_beta_unique      = eb_beta_unique,
      beta_temper_tau     = beta_temper_tau,
      beta_cap_q          = beta_cap_q,
      beta_var_floor_mult = beta_var_floor_mult,
      min_failures_for_EB = min_failures_for_EB
    )
    
    # Pass these directly to C++:
    # mu_m_t = hp_t$mu_m_t, mu_v_t = hp_t$mu_v_t,
    # b_m_t  = hp_t$b_m_t,  b_v_t  = hp_t$b_v_t
    
    CAPPMx.result <- common_atoms_cat_lognormal_shared_approx(
      dat_index = dat.index,
      trt_index = trt.index,
      st = response, nu = surv_ind,
      del = labels, eta = eta.cat,
      eta_cont = eta.cont,
      non_na_obs1 = non.na.inds,
      non_na_obs1_cont = non.na.inds_cont,
      eta_pred = matrix(nrow = 0, ncol = 1),
      eta_cont_pred = matrix(nrow = 0, ncol = 1),
      non_na_obs1_pred = list(0),
      non_na_obs1_cont_pred = list(0),
      nmix = nmix, ncat = ncats,
      a0 = a0, df0 = 1,
      mu_m_t = hp_t$mu_m_t,
      b_v_t = hp_t$b_v_t,
      b_m_t = hp_t$b_m_t,
      mu_v_t = hp_t$mu_v_t,
      del_range_lognorm_ref = del_range_response_1,
      nleapfrog_lognorm_ref = nleapfrog_response_1,
      del_range_lognorm_oth = del_range_response_2,
      nleapfrog_lognorm_oth = nleapfrog_response_2,
      alpha_hyper = c(1, 10),
      del_range_alp1 = del_range_alp1,
      nleapfrog_alp1 = nleapfrog_alp1,
      del_range_alp2 = del_range_alp2,
      nleapfrog_alp2 = nleapfrog_alp2,
      nrun = nrun, burn = burn, thin = thin,
      dat_index_old = result_CAPPMx$raw_Indices$dat.index + 1L,
      trt_convert = trt.convert,
      nj_val_all_old = result_CAPPMx$nj_val_all_mat[, new_M, drop = FALSE],
      nj_val_old = result_CAPPMx$nj_val_cube[, , new_M, drop = FALSE],
      nj_val_shared_old = result_CAPPMx$nj_val_shared_cube[, , new_M, drop = FALSE],
      noccu_old = result_CAPPMx$noccu_list[new_M],
      nobs_old = result_CAPPMx$nobs_cube[, , new_M, drop = FALSE],
      ss_j_x_old = result_CAPPMx$ss_j_x_cube[, , new_M, drop = FALSE],
      sum_j_x_old = result_CAPPMx$sum_j_x_cube[, , new_M, drop = FALSE],
      nj_x_old = result_CAPPMx$nj_x_cube[, , new_M, drop = FALSE],
      ss_survtime_old = result_CAPPMx$ss_survtime_cube[, , new_M, drop = FALSE],
      survtime_old = result_CAPPMx$survtime_cube[, , new_M, drop = FALSE],
      lognormal_mu_old = result_CAPPMx$Lognormal_Mu_Cube[new_M, , , drop = FALSE],
      lognormal_sig_old = result_CAPPMx$Lognormal_Sig_Cube[new_M, , , drop = FALSE],
      dirichlet_alpha_mat_old = result_CAPPMx$Dirichlet_params[new_M, , drop = FALSE],
      pimat_old = result_CAPPMx$picube[new_M, , , drop = FALSE],
      ## NEW (freeze control hyperparameters per draw)
      freeze_control = freeze_control,
      mu0_control_draws   = mu0_control_draws,      # vector length = length(new_M)
      beta0_control_draws = beta0_control_draws     # vector length = length(new_M)
    )
    
    CAPPMx.result$Treatment_Levels  <- trt_levels
    CAPPMx.result$Treatment_Indices <- seq_along(trt_levels)
    CAPPMx.result$Data_Levels       <- dat_levels
    CAPPMx.result$Data_Indices      <- seq_along(dat_levels)
    dat_trt.index <- cbind(dat.index, trt.index)
    colnames(dat_trt.index) <- c("Data_Indices", "Treatment_Indices")
    CAPPMx.result$Combined_Indices <- unique(dat_trt.index) + 1L
    CAPPMx.result$Variable_specifications <- result_CAPPMx$Variable_specifications
    
    return(CAPPMx.result)
  } else {
    if (length(result_CAPPMx) != 30)
      stop("You might be missing some essential output from Stage 1.")
    nmix <- nrow(result_CAPPMx$nj_val_all_mat)
    
    if (is.null(ref_trt)) stop("Please specify a reference treatment level.")
    if (!is.data.frame(input_df))
      stop("input_df must be a dataframe with response, censor, trt, dat, and covariates columns.")
    if (!is.data.frame(input_df_pred))
      stop("input_df_pred must be a dataframe containing covariate columns for prediction.")
    
    stage1_cat_levels <- result_CAPPMx$Cat_Levels
    isnullcatinput  <- is.null(input_specs$cat_vars)
    isnullcontinput <- is.null(input_specs$cont_vars)
    if (isnullcatinput & isnullcontinput)
      stop("You must have at least one categorical or continuous covariate.")
    
    stage1_levels <- result_CAPPMx$Data_Levels
    in_df_levels  <- unique(as.vector(input_df[[input_specs$dat_type]]))
    curr_candidates <- setdiff(in_df_levels, stage1_levels)
    if (length(curr_candidates) != 1) {
      stop("Could not infer a unique current cohort from '", input_specs$dat_type,
           "'. Exactly one level in input_df must be new (not in Stage 1).")
    }
    curr_dat_name <- curr_candidates
    rwd_dat_names <- stage1_levels
    
    is_rwd <- input_df[[input_specs$dat_type]] %in% rwd_dat_names
    
    input_df$dat.index <- as.numeric(factor(
      as.vector(input_df[[input_specs$dat_type]]),
      levels = c(curr_dat_name, rwd_dat_names),
      labels = 1:(1 + length(rwd_dat_names))
    )) - 1L
    
    ref_trt_name    <- ref_trt
    other_trt_names <- setdiff(
      unique(c(as.vector(input_df[[input_specs$trt_type]]),
               result_CAPPMx$Treatment_Levels)),
      ref_trt_name
    )
    trt_levels <- c(ref_trt_name, other_trt_names)
    dat_levels <- c(curr_dat_name, rwd_dat_names)
    
    input_df$trt.index <- as.numeric(factor(
      as.vector(input_df[[input_specs$trt_type]]),
      levels = trt_levels,
      labels = seq_along(trt_levels)
    )) - 1L
    
    dat.trt <- input_df[!is_rwd, , drop = FALSE]
    dat.rwd <- input_df[ is_rwd, , drop = FALSE]
    
    if (nrow(dat.trt) == 0L) stop("No rows found for the current cohort in input_df.")
    
    # --- Build trt.convert (NEW before initializer!) ---
    trt.convert <- cbind(
      result_CAPPMx$raw_Indices$trt.index,  # OLD trt index (Stage 1 space, 0-based)
      as.numeric(factor(
        factor(result_CAPPMx$raw_Indices$trt.index,
               levels = result_CAPPMx$Treatment_Indices - 1L,
               labels = result_CAPPMx$Treatment_Levels),
        levels = trt_levels,
        labels = seq_along(trt_levels)
      )) - 1L                                # NEW trt index (Stage 2 space, 0-based)
    )
    colnames(trt.convert) <- c("old_trt_index","New_trt_index")
    trt.convert <- unique(trt.convert)
    storage.mode(trt.convert) <- "integer"
    
    ## ---- rebuild covariates in Stage-1’s encoding/scaling ----
    spec_stage1 <- result_CAPPMx$Variable_specifications
    # stage1_cat_levels <- result_CAPPMx$Cat_Levels
    stage1_cats  <- spec_stage1$cat_vars  %||% character(0)
    stage1_conts <- spec_stage1$cont_vars %||% character(0)
    
    ## --- Categorical covariates (0-based, in Stage-1 order) ---
    build_cat_block <- function(df_rows) {
      if (length(stage1_cats) == 0L) {
        mat <- matrix(NA_real_, nrow(df_rows), 1L)
        colnames(mat) <- NULL
        return(mat)
      }
      cols <- lapply(stage1_cats, function(v) {
        if (!v %in% names(stage1_cat_levels)) {
          stop("Stage 1 Cat_Levels missing variable '", v, "'.")
        }
        if (v %in% names(df_rows)) {
          f <- factor(df_rows[[v]], levels = stage1_cat_levels[[v]])
          ## If new data has unseen levels (non-NA in data but NA after factoring), error:
          if (any(is.na(f) & !is.na(df_rows[[v]]))) {
            unseen <- unique(df_rows[[v]][is.na(f) & !is.na(df_rows[[v]])])
            stop("Unseen level(s) for '", v, "': ", paste(unseen, collapse = ", "))
          }
          as.numeric(f) - 1L
        } else {
          ## variable absent in new data -> keep alignment: all NA
          rep(NA_real_, nrow(df_rows))
        }
      })
      mat <- as.matrix(data.frame(cols, check.names = FALSE))
      colnames(mat) <- stage1_cats
      mat
    }
    
    cat_cov_trt <- build_cat_block(dat.trt)
    cat_cov_rwd <- build_cat_block(dat.rwd)
    
    ## --- Continuous covariates (raw scale first, Stage-1 order) ---
    build_cont_block <- function(df_rows) {
      if (length(stage1_conts) == 0L) {
        mat <- matrix(NA_real_, nrow(df_rows), 1L)
        colnames(mat) <- NULL
        return(mat)
      }
      cols <- lapply(stage1_conts, function(v) {
        if (v %in% names(df_rows)) as.numeric(df_rows[[v]]) else rep(NA_real_, nrow(df_rows))
      })
      mat <- as.matrix(data.frame(cols, check.names = FALSE))
      colnames(mat) <- stage1_conts
      mat
    }
    
    cont_cov_trt <- build_cont_block(dat.trt)
    cont_cov_rwd <- build_cont_block(dat.rwd)
    
    ## --- Stack in (current; historical) order used below ---
    eta.cat  <- rbind(cat_cov_trt,  cat_cov_rwd)
    eta.cont <- rbind(cont_cov_trt, cont_cov_rwd)
    rownames(eta.cat) <- rownames(eta.cont) <- NULL
    
    ## --- Non-missing indices per row (0-based) ---
    non.na.inds      <- apply(eta.cat,  1, function(x) which(!is.na(x)) - 1L, simplify = FALSE)
    non.na.inds_cont <- apply(eta.cont, 1, function(x) which(!is.na(x)) - 1L, simplify = FALSE)
    
    ## --- Scale continuous variables using Stage 1 summaries (aligned by name) ---
    if (length(stage1_conts) > 0L) {
      means_stage1 <- result_CAPPMx$Cont_Scale_Params$means[stage1_conts]
      sds_stage1   <- result_CAPPMx$Cont_Scale_Params$sds[stage1_conts]
      ## guard: replace 0/NA sds with 1 to avoid blow-ups
      sds_stage1[!is.finite(sds_stage1) | sds_stage1 == 0] <- 1
      
      ## sweep expects numeric; keep NAs as NAs
      eta.cont <- sweep(eta.cont, 2, means_stage1, "-")
      eta.cont <- sweep(eta.cont, 2, sds_stage1,   "/")
    } else {
      ## keep single NA column unchanged
    }
    
    ## --- ncats comes from Stage 1 and already matches Stage-1 cat order ---
    ncats <- result_CAPPMx$ncats
    
    response <- log(c(dat.trt[[input_specs$response]], dat.rwd[[input_specs$response]]))
    surv_ind <-      c(dat.trt[[input_specs$censor_ind]], dat.rwd[[input_specs$censor_ind]])
    
    map_t_old <- function(t_new) {
      hit <- which(trt.convert[, "New_trt_index"] == t_new)
      if (length(hit) == 0L) return(-1L)
      as.integer(trt.convert[hit[1L], "old_trt_index"])
    }
    
    # --- rows that belong to "current" cohort (to initialise) ---
    n_trt_rows <- nrow(dat.trt)
    n_rwd_rows <- nrow(dat.rwd)
    n_all      <- n_trt_rows + n_rwd_rows
    
    ## ===== Init via first half draws; run Stage 2 on last half draws =====
    suppressWarnings(suppressMessages({
      library(matrixStats)
      library(mcclust)   # for minbinder
    }))
    
    ## MCMC dims (we'll still use M and K from Stage 1)
    M  <- dim(result_CAPPMx$Lognormal_Mu_Cube)[1]
    K  <- dim(result_CAPPMx$Lognormal_Mu_Cube)[2]
    stopifnot(K == nrow(result_CAPPMx$nj_val_all_mat))
    
    ## choose the first half for initialisation, last half for the sampler later
    init_ids <- 0:((M/2) - 1L)               # 0-based
    # run_ids  <- (M - (M/2)):(M - 1L)         # 0-based
    
    ## dims
    N <- nrow(dat.trt)
    R <- length(init_ids)
    k_cat  <- ncol(eta.cat)
    k_cont <- ncol(eta.cont)
    
    if (length(stage1_cats) == 0L) k_cat  <- 0L
    if (length(stage1_conts) == 0L) k_cont <- 0L
    
    ## two cubes: [R x K x N]
    cov_loglik_cat  <- array(0, dim = c(R, K, N))
    cov_loglik_cont <- array(0, dim = c(R, K, N))
    
    ## Stage-1 hyper for continuous predictive
    df_x    <- 1.0
    alpha_x <- k_cont + 30.0
    beta_x  <- 1.0
    mu_x    <- 0.0
    
    ## Pull Stage-1 covariate tallies
    nobs_cube  <- result_CAPPMx$nobs_cube      # [K x k_cat  x M]
    nj_x_cube  <- result_CAPPMx$nj_x_cube      # [K x k_cont x M]
    sum_x_cube <- result_CAPPMx$sum_j_x_cube   # [K x k_cont x M]
    ss_x_cube  <- result_CAPPMx$ss_j_x_cube    # [K x k_cont x M]
    noccu_list <- result_CAPPMx$noccu_list     # length M; [[m]][[j]][[v]] vector over levels
    
    ## fill cubes
    cov_out <- stage2_cov_loglik_precompute(
      init_ids        = init_ids,
      eta_cat         = eta.cat,
      eta_cont        = eta.cont,
      non_na_obs1     = non.na.inds,
      non_na_obs1_cont= non.na.inds_cont,
      nobs_cube       = result_CAPPMx$nobs_cube,
      nj_x_cube       = result_CAPPMx$nj_x_cube,
      sum_x_cube      = result_CAPPMx$sum_j_x_cube,
      ss_x_cube       = result_CAPPMx$ss_j_x_cube,
      noccu_old       = result_CAPPMx$noccu_list,
      ncat            = ncats,
      df_x            = df_x,
      alpha_x         = alpha_x,
      mu_x            = mu_x,
      beta_x          = beta_x
    )
    
    cov_loglik_cat  <- cov_out$cov_loglik_cat   # reorder to [R,K,N]
    cov_loglik_cont <- cov_out$cov_loglik_cont
    
    ## If there are no cats/conts, the arrays above are already zeros with correct dims.
    
    ## ---- call the covariates-only initializer ----
    Z_init <- stage2_init_alloc_cov_only(
      cov_loglik_cat  = cov_loglik_cat,
      cov_loglik_cont = cov_loglik_cont,
      init_ids        = init_ids,
      nmix            = K,
      rng_seed        = 500L
    )
    ## Z_init is [R x N] with 0-based labels
    
    
    # --- consensus init via PSM + MinBinder ---
    # Z_init is [R x n_trt_rows] with 0-based labels; convert to 1-based for mcclust
    Z1 <- Z_init + 1L
    # co-clustering (PSM)
    psm <- mcclust::comp.psm(Z1)
    init_labels_trt_1b <- mcclust::minbinder(psm, max.k = min(15, K))$cl
    init_labels_trt    <- init_labels_trt_1b - 1L  # back to 0-based
    
    # assemble initial labels for full data (current + historical)
    labels1 <- init_labels_trt
    labels2 <- integer(n_rwd_rows)  # historical can be ignored for init; sampler will use Stage-1
    labels  <- c(labels1, labels2)
    
    ## --- Prediction covariates (0-based cat, Stage-1 order; cont raw then scaled) ---
    
    spec_stage1 <- result_CAPPMx$Variable_specifications
    stage1_cat_levels <- result_CAPPMx$Cat_Levels
    stage1_cats  <- spec_stage1$cat_vars  %||% character(0)
    stage1_conts <- spec_stage1$cont_vars %||% character(0)
    
    ## Categorical (0-based codes, Stage-1 level map)
    if (length(stage1_cats) == 0L) {
      cat_cov_pred <- matrix(NA_real_, nrow(input_df_pred), 1L)
      colnames(cat_cov_pred) <- NULL
    } else {
      cat_cov_pred <- as.matrix(data.frame(lapply(stage1_cats, function(v) {
        if (!v %in% names(stage1_cat_levels)) {
          stop("Prediction variable '", v, "' missing from Stage 1 Cat_Levels.")
        }
        if (v %in% names(input_df_pred)) {
          f <- factor(input_df_pred[[v]], levels = stage1_cat_levels[[v]])
          if (any(is.na(f) & !is.na(input_df_pred[[v]]))) {
            unseen <- unique(input_df_pred[[v]][is.na(f) & !is.na(input_df_pred[[v]])])
            stop("Unseen level(s) in prediction for '", v, "': ", paste(unseen, collapse = ", "))
          }
          as.numeric(f) - 1L
        } else {
          rep(NA_real_, nrow(input_df_pred))  # keep alignment if var absent
        }
      }), check.names = FALSE))
      colnames(cat_cov_pred) <- stage1_cats
    }
    
    ## Continuous (Stage-1 order; raw then scaled with Stage-1 params)
    if (length(stage1_conts) == 0L) {
      cont_cov_pred <- matrix(NA_real_, nrow(input_df_pred), 1L)
      colnames(cont_cov_pred) <- NULL
    } else {
      cont_cov_pred <- as.matrix(data.frame(lapply(stage1_conts, function(v) {
        if (v %in% names(input_df_pred)) as.numeric(input_df_pred[[v]]) else rep(NA_real_, nrow(input_df_pred))
      }), check.names = FALSE))
      colnames(cont_cov_pred) <- stage1_conts
      
      means_stage1 <- result_CAPPMx$Cont_Scale_Params$means[stage1_conts]
      sds_stage1   <- result_CAPPMx$Cont_Scale_Params$sds[stage1_conts]
      sds_stage1[!is.finite(sds_stage1) | sds_stage1 == 0] <- 1  # guard
      cont_cov_pred <- sweep(cont_cov_pred, 2, means_stage1, "-")
      cont_cov_pred <- sweep(cont_cov_pred, 2, sds_stage1,   "/")
    }
    
    ## Non-missing indices per row (0-based)
    non.na.inds_pred      <- apply(cat_cov_pred,  1, function(row) which(!is.na(row)) - 1L, simplify = FALSE)
    non.na.inds_cont_pred <- apply(cont_cov_pred, 1, function(row) which(!is.na(row)) - 1L, simplify = FALSE)
    
    dat.index <- c(dat.trt$dat.index, result_CAPPMx$raw_Indices$dat.index + 1L)
    trt.index <- c(
      dat.trt$trt.index,
      as.numeric(factor(
        factor(result_CAPPMx$raw_Indices$trt.index,
               levels = result_CAPPMx$Treatment_Indices - 1L,
               labels = result_CAPPMx$Treatment_Levels),
        levels = trt_levels,
        labels = seq_along(trt_levels)
      )) - 1L
    )
    
    burn <- 0; thin <- 5; nrun <- (nrow(result_CAPPMx$picube) - burn)/2
    
    # right before the C++ call
    storage.mode(eta.cat)  <- "integer"
    storage.mode(eta.cont) <- "double"
    storage.mode(dat.index) <- "integer"
    storage.mode(trt.index) <- "integer"
    
    ## discarding half
    new_M <- (M - (M/2) + 1):(M)
    
    ## --- map Stage-2 control (new index 0) back to Stage-1 old index ---
    ctl_row <- which(trt.convert[, "New_trt_index"] == 0L)
    if (length(ctl_row) != 1L) stop("Could not identify the control row in trt.convert.")
    old_ctl_0b <- trt.convert[ctl_row, "old_trt_index"]     # 0-based
    old_ctl_1b <- old_ctl_0b + 1L                           # 1-based for array indexing
    
    ## --- pull per-draw μ0,0 and β0,0 from Stage 1 (last-half only) ---
    H <- result_CAPPMx$Lognormal_hyperparams  # [M x 2 x T_old]; [,1]=mu0_t, [,2]=beta0_t
    mu0_control_draws   <- as.numeric(H[new_M, 1, old_ctl_1b])   # length = length(new_M)
    beta0_control_draws <- as.numeric(H[new_M, 2, old_ctl_1b])   # positive, natural scale
    
    ## small guardrails
    if (any(!is.finite(mu0_control_draws)))   stop("Non-finite control μ0 draws.")
    if (any(!is.finite(beta0_control_draws) | beta0_control_draws <= 0))
      stop("Non-finite or non-positive control β0 draws.")
    
    # ------------------------------------------------------------------------------
    # compute_stage2_hyperpriors_t()
    #
    # Purpose:
    #   Construct treatment-specific lognormal hyperparameters (μ₀,t, β₀,t)
    #   for Stage 2 of the approximate model.
    #
    # Behavior:
    #   • For overlapping (shared) treatment arms — i.e., arms present in both
    #     Stage 1 and Stage 2 — this function borrows empirical hyperparameters
    #     from the Stage 1 posterior (computed from the Lognormal_hyperparams cube).
    #
    #   • For Stage-2–only (unique) treatment arms — e.g., new experimental drugs
    #     not seen in Stage 1 — it rebuilds the hyperparameters from the current
    #     Stage-2 data using the same recipe as Stage 1:
    #         - μ₀,t hyperparameters (mean, variance) estimated from log-times
    #           among uncensored failures in that arm.
    #         - β₀,t hyperparameters constructed from target moments:
    #             E[β]  = cc × (a₀ − 1)
    #             Var[β] = s²_β = 5
    #           converted to log-scale mean/variance (b_m_t, b_v_t).
    #
    #   • Final values are sanitized to ensure positivity and finiteness.
    #     Returns a list with (μ_m_t, μ_v_t, b_m_t, b_v_t) and diagnostics.
    #
    # Inputs:
    #   - result_CAPPMx : Stage 1 fitted object containing Lognormal_hyperparams
    #   - trt.convert   : 0-based mapping of Stage-1 to Stage-2 treatment indices
    #   - T_new         : total number of treatments in Stage 2
    #   - input_df      : Stage-2 data with columns trt, os, os_status
    #   - ref_trt       : reference treatment label (usually "Control")
    #   - a0, cc, s2_beta : parameters for the Stage-1 prior recipe (used only for
    #                       Stage-2-unique arms)
    #
    # Output:
    #   A list with:
    #       mu_m_t, mu_v_t, b_m_t, b_v_t,
    #       diagnostics (overlap map, pooled stats, arm roles)
    #
    # ------------------------------------------------------------------------------
    compute_stage2_hyperpriors_t <- function(
    result_CAPPMx,
    trt.convert,               # 0-based [old, new] mapping
    T_new,                     # number of Stage-2 treatments
    input_df,                  # Stage-2 dataset
    input_specs,               # app-provided column names
    ref_trt = "Control",
    use_last_half = TRUE,
    a0 = 10.1,                 # used only for Stage-2-unique arms
    cc = 2,                    # used only for Stage-2-unique arms
    s2_beta = 5,               # target Var[beta] used in Stage-1 recipe
    eb_beta_unique = FALSE,
    beta_temper_tau     = 10,
    beta_cap_q          = c(0.10, 0.90),
    beta_var_floor_mult = 1.0,
    min_failures_for_EB = 3
    ) {
      if (is.null(result_CAPPMx$Lognormal_hyperparams)) {
        stop("Stage 1 output missing $Lognormal_hyperparams.")
      }
      
      H <- result_CAPPMx$Lognormal_hyperparams  # [M x 2 x T_old]; [,1]=mu0_t, [,2]=beta0_t
      
      dH <- dim(H)
      if (length(dH) != 3L || dH[2] != 2L) {
        stop("$Lognormal_hyperparams must be [M x 2 x T_old].")
      }
      
      M <- dH[1]
      T_old <- dH[3]
      
      # Required column names from input_specs
      trt_col <- input_specs$trt_type
      os_col <- input_specs$response
      censor_col <- input_specs$censor_ind
      
      required_model_cols <- c(trt_col, os_col, censor_col)
      missing_model_cols <- setdiff(required_model_cols, names(input_df))
      
      if (length(missing_model_cols) > 0) {
        stop(
          "input_df is missing required model-fitting column(s): ",
          paste(missing_model_cols, collapse = ", ")
        )
      }
      
      # Draws to use
      draw_idx <- if (use_last_half) {
        seq.int(floor(M / 2) + 1L, M)
      } else {
        seq_len(M)
      }
      
      # Tidy treatment conversion table
      tc <- as.data.frame(trt.convert)
      names(tc) <- c("old", "new")
      tc$new[is.na(tc$new)] <- -1L
      tc <- unique(tc)
      
      # Validate T_new
      stopifnot(length(T_new) == 1L, is.finite(T_new), T_new >= 1L)
      
      # Pre-allocate in Stage-2 treatment order
      mu_m_t <- mu_v_t <- b_m_t <- b_v_t <- rep(NA_real_, T_new)
      
      # ---------- SHARED / OVERLAPPING ARMS ----------
      overlap <- tc[
        tc$new >= 0L & tc$new < T_new &
          tc$old >= 0L & tc$old < T_old,
        ,
        drop = FALSE
      ]
      
      if (nrow(overlap) == 0L) {
        stop("No overlapping treatments between Stage 1 and Stage 2.")
      }
      
      per_arm <- lapply(seq_len(nrow(overlap)), function(ii) {
        t_old0 <- overlap$old[ii]   # 0-based Stage-1 index
        t_new  <- overlap$new[ii]   # 0-based Stage-2 index
        t_old  <- t_old0 + 1L       # 1-based index for R array
        
        mu_draws   <- as.numeric(H[draw_idx, 1, t_old])
        beta_draws <- as.numeric(H[draw_idx, 2, t_old])
        
        ok_mu <- is.finite(mu_draws)
        ok_b  <- is.finite(beta_draws) & beta_draws > 0
        
        if (!any(ok_mu) || !any(ok_b)) {
          return(list(t_old = t_old0, t_new = t_new, ok = FALSE))
        }
        
        list(
          t_old     = t_old0,
          t_new     = t_new,
          ok        = TRUE,
          mu_mean   = mean(mu_draws[ok_mu]),
          mu_var    = stats::var(mu_draws[ok_mu]),
          logb_mean = mean(log(beta_draws[ok_b])),
          logb_var  = stats::var(log(beta_draws[ok_b]))
        )
      })
      
      keep <- vapply(per_arm, function(x) isTRUE(x$ok), logical(1))
      per_arm <- per_arm[keep]
      
      if (length(per_arm) == 0L) {
        stop("No valid posterior draws found for overlapping treatments.")
      }
      
      for (pa in per_arm) {
        idx <- pa$t_new + 1L
        
        mu_m_t[idx] <- pa$mu_mean
        mu_v_t[idx] <- if (is.finite(pa$mu_var) && pa$mu_var > 0) {
          pa$mu_var
        } else {
          NA_real_
        }
        
        b_m_t[idx] <- pa$logb_mean
        b_v_t[idx] <- if (is.finite(pa$logb_var) && pa$logb_var > 0) {
          pa$logb_var
        } else {
          NA_real_
        }
      }
      
      # Pooled fallbacks from overlapping arms
      mu_means   <- vapply(per_arm, `[[`, numeric(1), "mu_mean")
      mu_vars    <- vapply(per_arm, `[[`, numeric(1), "mu_var")
      logb_means <- vapply(per_arm, `[[`, numeric(1), "logb_mean")
      logb_vars  <- vapply(per_arm, `[[`, numeric(1), "logb_var")
      
      w_mu <- ifelse(is.finite(mu_vars) & mu_vars > 0, 1 / mu_vars, 0)
      w_logb <- ifelse(is.finite(logb_vars) & logb_vars > 0, 1 / logb_vars, 0)
      
      mu_pool <- if (sum(w_mu) > 0) {
        sum(w_mu * mu_means) / sum(w_mu)
      } else {
        mean(mu_means)
      }
      
      logb_pool <- if (sum(w_logb) > 0) {
        sum(w_logb * logb_means) / sum(w_logb)
      } else {
        mean(logb_means)
      }
      
      base_mu_v <- mean(
        ifelse(is.finite(mu_vars) & mu_vars > 0, mu_vars, NA),
        na.rm = TRUE
      )
      
      base_logb_v <- mean(
        ifelse(is.finite(logb_vars) & logb_vars > 0, logb_vars, NA),
        na.rm = TRUE
      )
      
      base_mu_v <- if (is.finite(base_mu_v) && base_mu_v > 0) {
        base_mu_v
      } else {
        1e-6
      }
      
      base_logb_v <- if (is.finite(base_logb_v) && base_logb_v > 0) {
        base_logb_v
      } else {
        1e-6
      }
      
      bad_mu_v <- !is.finite(mu_v_t) | mu_v_t <= 0
      bad_b_v  <- !is.finite(b_v_t)  | b_v_t  <= 0
      
      if (any(bad_mu_v)) {
        mu_v_t[bad_mu_v] <- base_mu_v
      }
      
      if (any(bad_b_v)) {
        b_v_t[bad_b_v] <- base_logb_v
      }
      
      # ---------- CONTROL-BASED CAPS FOR BETA ----------
      ctl_row <- which(tc$new == 0L)
      
      if (length(ctl_row) == 1L &&
          tc$old[ctl_row] >= 0L &&
          tc$old[ctl_row] < T_old) {
        
        t_old_ctl <- tc$old[ctl_row] + 1L
        ctl_logbeta <- log(as.numeric(H[draw_idx, 2, t_old_ctl]))
        ctl_logbeta <- ctl_logbeta[is.finite(ctl_logbeta)]
        
        if (length(ctl_logbeta) >= 20L) {
          beta_cap_lo <- as.numeric(
            stats::quantile(exp(ctl_logbeta), probs = beta_cap_q[1])
          )
          beta_cap_hi <- as.numeric(
            stats::quantile(exp(ctl_logbeta), probs = beta_cap_q[2])
          )
        } else {
          beta_cap_lo <- exp(logb_pool)
          beta_cap_hi <- exp(logb_pool)
        }
      } else {
        beta_cap_lo <- exp(logb_pool)
        beta_cap_hi <- exp(logb_pool)
      }
      
      beta_logvar_floor <- max(base_logb_v * beta_var_floor_mult, 1e-6)
      
      # ---------- UNIQUE / STAGE-2-ONLY ARMS ----------
      overlap_new <- vapply(per_arm, `[[`, integer(1), "t_new")
      stage2_only_0b <- setdiff(seq_len(T_new) - 1L, unique(overlap_new))
      
      # Treatment labels in the uploaded/current dataset
      arm_levels <- levels(factor(input_df[[trt_col]]))
      
      if (length(arm_levels) < T_new) {
        arm_levels <- unique(c(arm_levels, ref_trt))
      }
      
      if (length(stage2_only_0b)) {
        
        # Original Stage-1 beta recipe fallback
        if (!isTRUE(eb_beta_unique)) {
          m <- cc * (a0 - 1)
          b_v_uni <- log1p(s2_beta / m^2)
          b_m_uni <- log(m) - 0.5 * b_v_uni
        }
        
        for (t_new in stage2_only_0b) {
          idx <- t_new + 1L
          
          t_lab <- if (!is.na(arm_levels[idx])) {
            arm_levels[idx]
          } else {
            ref_trt
          }
          
          arm_df <- input_df[
            as.character(input_df[[trt_col]]) == as.character(t_lab) &
              input_df[[censor_col]] == 1,
            ,
            drop = FALSE
          ]
          
          # ----- mu hyperparameters from Stage-2 failures -----
          if (nrow(arm_df) >= 2L) {
            mu_m_t[idx] <- mean(log(arm_df[[os_col]]), na.rm = TRUE)
            mu_v_t[idx] <- stats::var(log(arm_df[[os_col]]), na.rm = TRUE)
          } else if (nrow(arm_df) == 1L) {
            mu_m_t[idx] <- log(arm_df[[os_col]][1L])
            mu_v_t[idx] <- 1e-4
          } else {
            pool <- input_df[
              input_df[[censor_col]] == 1,
              ,
              drop = FALSE
            ]
            
            mu_m_t[idx] <- mean(log(pool[[os_col]]), na.rm = TRUE)
            mu_v_t[idx] <- stats::var(log(pool[[os_col]]), na.rm = TRUE)
          }
          
          if (!is.finite(mu_v_t[idx]) || mu_v_t[idx] <= 0) {
            mu_v_t[idx] <- 1e-4
          }
          
          # ----- beta hyperparameters -----
          if (isTRUE(eb_beta_unique)) {
            nf <- nrow(arm_df)
            
            s2_arm <- if (nf >= 2L) {
              stats::var(log(arm_df[[os_col]]), na.rm = TRUE)
            } else {
              NA_real_
            }
            
            w <- if (is.finite(s2_arm) && nf >= min_failures_for_EB) {
              nf / (nf + beta_temper_tau)
            } else {
              0
            }
            
            logb_EB <- if (is.finite(s2_arm) && s2_arm > 0) {
              log(s2_arm)
            } else {
              logb_pool
            }
            
            b_m_t[idx] <- w * logb_EB + (1 - w) * logb_pool
            
            beta_center <- exp(b_m_t[idx])
            
            if (is.finite(beta_cap_lo) &&
                is.finite(beta_cap_hi) &&
                beta_cap_hi > beta_cap_lo) {
              beta_center <- min(max(beta_center, beta_cap_lo), beta_cap_hi)
            }
            
            b_m_t[idx] <- log(beta_center)
            b_v_t[idx] <- beta_logvar_floor
          } else {
            b_m_t[idx] <- b_m_uni
            b_v_t[idx] <- b_v_uni
          }
        }
      }
      
      # ---------- FINAL SANITATION ----------
      mu_m_t[!is.finite(mu_m_t)] <- mu_pool
      b_m_t[!is.finite(b_m_t)] <- logb_pool
      
      mu_v_t[!is.finite(mu_v_t) | mu_v_t <= 0] <- base_mu_v
      b_v_t[!is.finite(b_v_t) | b_v_t <= 0] <- base_logb_v
      
      stopifnot(
        length(mu_m_t) == T_new,
        length(mu_v_t) == T_new,
        length(b_m_t)  == T_new,
        length(b_v_t)  == T_new
      )
      
      list(
        mu_m_t = mu_m_t,
        mu_v_t = mu_v_t,
        b_m_t  = b_m_t,
        b_v_t  = b_v_t,
        diagnostics = list(
          overlap_map = data.frame(
            old = vapply(per_arm, `[[`, integer(1), "t_old"),
            new = vapply(per_arm, `[[`, integer(1), "t_new")
          ),
          pooled = c(
            mu_pool = mu_pool,
            logb_pool = logb_pool,
            base_mu_v = base_mu_v,
            base_logb_v = base_logb_v
          ),
          stage2_only = stage2_only_0b,
          arm_role = {
            role <- rep("shared", T_new)
            if (length(stage2_only_0b)) {
              role[stage2_only_0b + 1L] <- "unique"
            }
            role
          }
        )
      )
    }
    
    T_new <- length(trt_levels)  # or max(trt.index) + 1L
    hp_t <- compute_stage2_hyperpriors_t(
      result_CAPPMx = result_CAPPMx,
      trt.convert   = trt.convert,
      T_new         = T_new,
      input_df      = input_df,
      input_specs   = input_specs,
      ref_trt       = ref_trt,
      use_last_half = TRUE,
      a0 = a0,
      cc = cc,
      s2_beta = 5,
      eb_beta_unique      = eb_beta_unique,
      beta_temper_tau     = beta_temper_tau,
      beta_cap_q          = beta_cap_q,
      beta_var_floor_mult = beta_var_floor_mult,
      min_failures_for_EB = min_failures_for_EB
    )
    
    # Pass these directly to C++:
    # mu_m_t = hp_t$mu_m_t, mu_v_t = hp_t$mu_v_t,
    # b_m_t  = hp_t$b_m_t,  b_v_t  = hp_t$b_v_t
    
    
    CAPPMx.result <- common_atoms_cat_lognormal_shared_approx(
      dat_index = dat.index,
      trt_index = trt.index,
      st = response, nu = surv_ind,
      del = labels, eta = eta.cat,
      eta_cont = eta.cont,
      non_na_obs1 = non.na.inds,
      non_na_obs1_cont = non.na.inds_cont,
      eta_pred = cat_cov_pred,
      eta_cont_pred = cont_cov_pred,
      non_na_obs1_pred = non.na.inds_pred,
      non_na_obs1_cont_pred = non.na.inds_cont_pred,
      nmix = nmix, ncat = ncats,
      a0 = a0, df0 = 1,
      mu_m_t = hp_t$mu_m_t,
      b_v_t = hp_t$b_v_t,
      b_m_t = hp_t$b_m_t,
      mu_v_t = hp_t$mu_v_t,
      del_range_lognorm_ref = del_range_response_1,
      nleapfrog_lognorm_ref = nleapfrog_response_1,
      del_range_lognorm_oth = del_range_response_2,
      nleapfrog_lognorm_oth = nleapfrog_response_2,
      alpha_hyper = c(1, 10),
      del_range_alp1 = del_range_alp1,
      nleapfrog_alp1 = nleapfrog_alp1,
      del_range_alp2 = del_range_alp2,
      nleapfrog_alp2 = nleapfrog_alp2,
      nrun = nrun, burn = burn, thin = thin,
      dat_index_old = result_CAPPMx$raw_Indices$dat.index + 1L,
      trt_convert = trt.convert,
      nj_val_all_old = result_CAPPMx$nj_val_all_mat[, new_M, drop = FALSE],
      nj_val_old = result_CAPPMx$nj_val_cube[, , new_M, drop = FALSE],
      nj_val_shared_old = result_CAPPMx$nj_val_shared_cube[, , new_M, drop = FALSE],
      noccu_old = result_CAPPMx$noccu_list[new_M],
      nobs_old = result_CAPPMx$nobs_cube[, , new_M, drop = FALSE],
      ss_j_x_old = result_CAPPMx$ss_j_x_cube[, , new_M, drop = FALSE],
      sum_j_x_old = result_CAPPMx$sum_j_x_cube[, , new_M, drop = FALSE],
      nj_x_old = result_CAPPMx$nj_x_cube[, , new_M, drop = FALSE],
      ss_survtime_old = result_CAPPMx$ss_survtime_cube[, , new_M, drop = FALSE],
      survtime_old = result_CAPPMx$survtime_cube[, , new_M, drop = FALSE],
      lognormal_mu_old = result_CAPPMx$Lognormal_Mu_Cube[new_M, , , drop = FALSE],
      lognormal_sig_old = result_CAPPMx$Lognormal_Sig_Cube[new_M, , , drop = FALSE],
      dirichlet_alpha_mat_old = result_CAPPMx$Dirichlet_params[new_M, , drop = FALSE],
      pimat_old = result_CAPPMx$picube[new_M, , , drop = FALSE],
      ## NEW (freeze control hyperparameters per draw)
      freeze_control = freeze_control,
      mu0_control_draws   = mu0_control_draws,      # vector length = length(new_M)
      beta0_control_draws = beta0_control_draws     # vector length = length(new_M)
    )
    
    CAPPMx.result$Treatment_Levels  <- trt_levels
    CAPPMx.result$Treatment_Indices <- seq_along(trt_levels)
    CAPPMx.result$Data_Levels       <- dat_levels
    CAPPMx.result$Data_Indices      <- seq_along(dat_levels)
    dat_trt.index <- cbind(dat.index, trt.index)
    colnames(dat_trt.index) <- c("Data_Indices", "Treatment_Indices")
    CAPPMx.result$Combined_Indices <- unique(dat_trt.index) + 1L
    CAPPMx.result$Variable_specifications <- result_CAPPMx$Variable_specifications
    
    return(CAPPMx.result)
    # ============= PREDICTION BRANCH =============
  }
  
}
