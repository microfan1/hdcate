HDCATE.conf.set_data <- function(HDCATE_instance, new_data) {
  HDCATE_instance$data <- new_data
  HDCATE_instance$local_weight <- rep(1, nrow(new_data))
  HDCATE_instance$n_obs = nrow(new_data)
}

HDCATE.conf.set_y <- function(HDCATE_instance, var_name) {
  HDCATE_instance$y_name <- var_name
}

HDCATE.conf.set_d <- function(HDCATE_instance, var_name) {
  HDCATE_instance$d_name <- var_name
}

HDCATE.conf.set_x <- function(HDCATE_instance, var_name) {
  HDCATE_instance$x_formula <- var_name
}

HDCATE.conf.set_condition_var <- function(HDCATE_instance, name=NA, min=NA, max=NA, step=NA) {
  if (!(is.null(name) || anyNA(name))) {
    HDCATE_instance$cond_var_name <- name
    HDCATE_instance$d <- length(name)
  }
  if (!(is.null(min) || anyNA(min))) {
    HDCATE_instance$cond_var_lower <- min
  }
  if (!(is.null(max) || anyNA(max))) {
    HDCATE_instance$cond_var_upper <- max
  }
  if (!(is.null(step) || anyNA(step))) {
    HDCATE_instance$cond_var_eval_step <- step
  }

  # validation
  v1 <- (!(is.null(HDCATE_instance$cond_var_name) || anyNA(HDCATE_instance$cond_var_name)))
  v2 <- (!(is.null(HDCATE_instance$cond_var_lower) || anyNA(HDCATE_instance$cond_var_lower)))
  v3 <- (!(is.null(HDCATE_instance$cond_var_upper) || anyNA(HDCATE_instance$cond_var_upper)))
  # step=0.01 by default (defined in hdcate R6 class)
  v4 <- (!(is.null(HDCATE_instance$cond_var_eval_step) || anyNA(HDCATE_instance$cond_var_eval_step)))
  if (!(v1 && v2 && v3)) {
    warning(paste0('Failed to update conditional variable to: name=', name, ', min=', min, ', max=', max, ', step=', step, '.'))
    stop('Missing some of these params for the conditional variable: `name`, `min` and `max`.')
  } else {
    HDCATE_instance$cond_interval <- seq(HDCATE_instance$cond_var_lower, HDCATE_instance$cond_var_upper, HDCATE_instance$cond_var_eval_step)
    HDCATE_instance$length_grid <- length(HDCATE_instance$cond_interval)
    message(paste0('Updated conditional variable to: name=', name, ', min=', min, ', max=', max, ', step=', step, '.'))
  }
}

HDCATE.conf.enable_cross_fitting <- function(HDCATE_instance, k_fold=5, folds=NA) {
  HDCATE_instance$use_cross_fitting <- TRUE

  if (!(is.null(folds) || anyNA(folds))) {
    k_fold <- length(folds)
    message(paste0('Manually set folds. Detect `k_fold`=', k_fold, '.'))
  }
  HDCATE_instance$k_fold <- k_fold
  HDCATE_instance$folds <- folds
  HDCATE_instance$propensity_score_model_list <- list()
  HDCATE_instance$treated_cond_exp_model_list <- list()
  HDCATE_instance$untreated_cond_exp_model_list <- list()
  HDCATE_instance$local_reg_model_list <- list()
  HDCATE_instance$hdcate <- NA
  HDCATE_instance$sigma_hat <- NA
  HDCATE_instance$eta_hat_list <- NA
  message(paste0('Updated: cross-fitting estimater is used. Num of folds=', k_fold, '.'))
}

HDCATE.conf.disable_cross_fitting <- function(HDCATE_instance){
  HDCATE_instance$use_cross_fitting <- FALSE
  HDCATE_instance$k_fold <- 5
  HDCATE_instance$folds <- NULL
  HDCATE_instance$propensity_score_model_list <- list()
  HDCATE_instance$treated_cond_exp_model_list <- list()
  HDCATE_instance$untreated_cond_exp_model_list <- list()
  HDCATE_instance$local_reg_model_list <- list()
  HDCATE_instance$hdcate <- NA
  HDCATE_instance$sigma_hat <- NA
  HDCATE_instance$eta_hat_list <- NA
  message(paste0('Updated: full-sample estimater is used.'))
}

HDCATE.conf.set_bw <- function(HDCATE_instance, value){
  if (value=='default') {
    HDCATE_instance$bw <- 'rule-of-thumb'
  } else {
    HDCATE_instance$bw <- value
  }
}

HDCATE.conf.set_local_weight <- function(HDCATE_instance, value){
  if (value=='default') {
    HDCATE_instance$local_weight <- rep(1, nrow(HDCATE_instance$data))
  } else {
    HDCATE_instance$local_weight <- value
  }
}
