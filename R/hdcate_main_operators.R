#' High-Dimensional Conditional Average Treatment Effects (HDCATE) Estimator
#' @importFrom R6 R6Class
#' @importFrom hdm rlassologit
#' @importFrom hdm rlasso
#' @importFrom locpol locPolSmootherC
#' @importFrom KernSmooth dpill
#' @importFrom stats bw.nrd0
#' @importFrom stats rnorm
#' @importFrom stats runif
#' @importFrom stats rexp
#' @importFrom caret createFolds
#'
#' @description
#' Use a two-step procedure to estimate the conditional average treatment effects (CATE) with potentially high-dimensional covariate(s).
#' Run `browseVignettes('hdcate')` to browse the user manual of this package.
#'
#' @param data data frame of the observed data
#' @param y_name variable name of the observed outcomes
#' @param d_name variable name of the treatment indicators
#' @param x_formula formula of the covariates
#'
#' @return An initialized `HDCATE` model (object), ready for estimation.
#' @export
#'
#' @examples
#' # get simulation data
#' n_obs <- 500  # Num of observations
#' n_var <- 100  # Num of observed variables
#' n_rel_var <- 4  # Num of relevant variables
#' data <- HDCATE.get_sim_data(n_obs, n_var, n_rel_var)
#' # conditional expectation model is misspecified
#' x_formula <- paste(paste0('X', c(2:n_var)), collapse ='+')
#' # for example, and alternatively, the propensity score model is misspecified
#' # x_formula <- paste(paste0('X', c(1:(n_var-1))), collapse ='+')
#'
#' # Example 1: full-sample estimator
#' # create a new HDCATE model
#' model <- HDCATE(data=data, y_name='Y', d_name='D', x_formula=x_formula)
#'
#' # estimate HDCATE function, inference, and plot
#' HDCATE.set_condition_var(model, 'X2', min=-1, max=1, step=0.01)
#' \donttest{
#' HDCATE.fit(model)
#' HDCATE.inference(model)
#' HDCATE.plot(model)
#' }
#'
#' # Example 2: cross-fitting estimator
#' # change above estimator to cross-fitting mode, 5 folds, for example.
#' HDCATE.use_cross_fitting(model, k_fold=5)
#'
#' # estimate HDCATE function, inference, and plot
#' HDCATE.set_condition_var(model, 'X2', min=-1, max=1, step=0.01)
#' \donttest{
#' HDCATE.fit(model)
#' HDCATE.inference(model)
#' HDCATE.plot(model)
#' }
HDCATE <- function(data, y_name, d_name, x_formula) {
  HDCATE_R6Class$new(data, y_name, d_name, x_formula)
}

#' Set the conditional variable in CATE
#'
#' @param HDCATE_model an object created via [HDCATE]
#' @param name name of the conditional variable
#' @param min minimum value of the conditional variable for evaluation
#' @param max maximum value of the conditional variable for evaluation
#' @param step minimum distance between two evaluation points
#'
#' @return None. The `HDCATE_model` is ready to fit.
#' @export
#'
#' @examples
#' # get simulation data
#' n_obs <- 500  # Num of observations
#' n_var <- 100  # Num of observed variables
#' n_rel_var <- 4  # Num of relevant variables
#' data <- HDCATE.get_sim_data(n_obs, n_var, n_rel_var)
#' # conditional expectation model is misspecified
#' x_formula <- paste(paste0('X', c(2:n_var)), collapse ='+')
#' # propensity score model is misspecified
#' # x_formula <- paste(paste0('X', c(1:(n_var-1))), collapse ='+')
#'
#' # create a new HDCATE model
#' model <- HDCATE(data=data, y_name='Y', d_name='D', x_formula=x_formula)
#'
#' HDCATE.set_condition_var(model, 'X2', min=-1, max=1, step=0.01)
HDCATE.set_condition_var <- function(HDCATE_model, name=NA, min=NA, max=NA, step=NA) {
  HDCATE_model$set_condition_var(name, min, max, step)
}

#' Fit the HDCATE function
#'
#' @param HDCATE_model an object created via [HDCATE]
#' @param verbose whether the verbose message is displayed, the default is `TRUE`
#'
#' @return None. The `HDCATE_model` is fitted.
#' @export
#'
#' @examples
#' # get simulation data
#' n_obs <- 500  # Num of observations
#' n_var <- 100  # Num of observed variables
#' n_rel_var <- 4  # Num of relevant variables
#' data <- HDCATE.get_sim_data(n_obs, n_var, n_rel_var)
#' # conditional expectation model is misspecified
#' x_formula <- paste(paste0('X', c(2:n_var)), collapse ='+')
#' # propensity score model is misspecified
#' # x_formula <- paste(paste0('X', c(1:(n_var-1))), collapse ='+')
#'
#' # create a new HDCATE model
#' model <- HDCATE(data=data, y_name='Y', d_name='D', x_formula=x_formula)
#'
#' HDCATE.set_condition_var(model, 'X2', min=-1, max=1, step=0.01)
#'
#' \donttest{
#' HDCATE.fit(model)
#' }
HDCATE.fit <- function(HDCATE_model, verbose=TRUE) {
  HDCATE_model$fit(verbose=verbose)
}

#' Construct uniform confidence bands
#'
#' @param HDCATE_model an object created via [HDCATE]
#' @param sig_level a (vector of) significant level, such as 0.01, or c(0.01, 0.05, 0.10)
#' @param n_rep_boot repeat n times for bootstrap, the default is 1000
#' @param verbose whether the verbose message is displayed, the default is `FALSE`
#'
#' @return None. The HDCATE confidence bands are constructed.
#' @export
#'
#' @examples
#' # get simulation data
#' n_obs <- 500  # Num of observations
#' n_var <- 100  # Num of observed variables
#' n_rel_var <- 4  # Num of relevant variables
#' data <- HDCATE.get_sim_data(n_obs, n_var, n_rel_var)
#' # conditional expectation model is misspecified
#' x_formula <- paste(paste0('X', c(2:n_var)), collapse ='+')
#' # propensity score model is misspecified
#' # x_formula <- paste(paste0('X', c(1:(n_var-1))), collapse ='+')
#'
#' # create a new HDCATE model
#' model <- HDCATE(data=data, y_name='Y', d_name='D', x_formula=x_formula)
#'
#' HDCATE.set_condition_var(model, 'X2', min=-1, max=1, step=0.01)
#'
#' \donttest{
#' HDCATE.fit(model)
#' HDCATE.inference(model)
#' }
HDCATE.inference <- function(HDCATE_model, sig_level=0.01, n_rep_boot=1000, verbose=FALSE) {
  HDCATE_model$inference(sig_level=sig_level, boot_method='normal', n_rep_boot=n_rep_boot, verbose=verbose)
}

#' Plot HDCATE function and the uniform confidence bands
#'
#' @param HDCATE_model an object created via [HDCATE]
#' @param output_pdf if `TRUE`, the plot will be saved as a PDF file, the default is `FALSE`
#' @param pdf_name file name when `output_pdf=TRUE`
#' @param include_band if `TRUE`, plot the uniform confidence bands (need: [HDCATE.inference] was called before)
#' @param test_side `'both'`, `'left'` or `'right'`, i.e. 2-side test or one-side test
#' @param y_axis_min minimum value of the Y axis to plot in the graph, the default is `auto`
#' @param y_axis_max maximum value of the Y axis to plot in the graph, the default is `auto`
#' @param display.hdcate the name of HDCATE function in the legend, the default is 'HDCATEF'
#' @param display.ate the name of average treatment effect in the legend, the default is 'ATE'
#' @param display.siglevel the name of the significant level for confidence bands in the legend, the default is 'sig_level'
#'
#' @return None. A plot will be shown or saved as PDF.
#' @export
#'
#' @examples
#' # get simulation data
#' n_obs <- 500  # Num of observations
#' n_var <- 100  # Num of observed variables
#' n_rel_var <- 4  # Num of relevant variables
#' data <- HDCATE.get_sim_data(n_obs, n_var, n_rel_var)
#' # conditional expectation model is misspecified
#' x_formula <- paste(paste0('X', c(2:n_var)), collapse ='+')
#' # propensity score model is misspecified
#' # x_formula <- paste(paste0('X', c(1:(n_var-1))), collapse ='+')
#'
#' # create a new HDCATE model
#' model <- HDCATE(data=data, y_name='Y', d_name='D', x_formula=x_formula)
#'
#' HDCATE.set_condition_var(model, 'X2', min=-1, max=1, step=0.01)
#' \donttest{
#' HDCATE.fit(model)
#' HDCATE.inference(model)
#' HDCATE.plot(model)
#' }
HDCATE.plot <- function(HDCATE_model, output_pdf=FALSE, pdf_name='hdcate_plot.pdf', include_band=TRUE, test_side='both', y_axis_min='auto', y_axis_max='auto', display.hdcate='HDCATEF', display.ate='ATE', display.siglevel='sig_level') {
  HDCATE_model$plot(output_pdf=output_pdf, pdf_name=pdf_name, include_band=include_band, test_side=test_side, y_axis_min=y_axis_min, y_axis_max=y_axis_max, display.hdcate=display.hdcate, display.ate=display.ate, display.siglevel=display.siglevel)
}


#' Use k-fold cross-fitting estimator
#'
#' @param model an object created via [HDCATE]
#' @param k_fold number of folds
#' @param folds you can manually set the folds, should be a list of index vector
#'
#' @return None.
#' @export
#'
#' @examples
#' # get simulation data
#' n_obs <- 500  # Num of observations
#' n_var <- 100  # Num of observed variables
#' n_rel_var <- 4  # Num of relevant variables
#' data <- HDCATE.get_sim_data(n_obs, n_var, n_rel_var)
#' # conditional expectation model is misspecified
#' x_formula <- paste(paste0('X', c(2:n_var)), collapse ='+')
#' # propensity score model is misspecified
#' # x_formula <- paste(paste0('X', c(1:(n_var-1))), collapse ='+')
#'
#' # create a new HDCATE model
#' model <- HDCATE(data=data, y_name='Y', d_name='D', x_formula=x_formula)
#'
#' # for example, use 5-fold cross-fitting estimator
#' HDCATE.use_cross_fitting(model, k_fold=5)
#'
#' # alternatively, pass a list of index vector to the third argument to set the folds manually,
#' # in this case, the second argument k_fold is auto detected, you can pass any value to it.
#' HDCATE.use_cross_fitting(model, k_fold=2, folds=list(c(1:250), c(251:500)))
HDCATE.use_cross_fitting <- function(model, k_fold=5, folds=NULL) {
  HDCATE.conf.enable_cross_fitting(model, k_fold, folds)
}

#' Use full-sample estimator
#'
#' @description This is the default mode when creating a model via [HDCATE]
#'
#' @param model an object created via [HDCATE]
#'
#' @return None.
#' @export
#'
#' @examples
#' # get simulation data
#' n_obs <- 500  # Num of observations
#' n_var <- 100  # Num of observed variables
#' n_rel_var <- 4  # Num of relevant variables
#' data <- HDCATE.get_sim_data(n_obs, n_var, n_rel_var)
#' # conditional expectation model is misspecified
#' x_formula <- paste(paste0('X', c(2:n_var)), collapse ='+')
#' # propensity score model is misspecified
#' # x_formula <- paste(paste0('X', c(1:(n_var-1))), collapse ='+')
#'
#' # create a new HDCATE model
#' model <- HDCATE(data=data, y_name='Y', d_name='D', x_formula=x_formula)
#'
#' HDCATE.use_full_sample(model)
HDCATE.use_full_sample <- function(model) {
  HDCATE.conf.disable_cross_fitting(model)
}

#' Set user-defined first-stage estimating methods
#'
#' @description Set user-defined ML methods (such as random forests, elastic-net, boosting) to run the first-stage estimation.
#'
#' @param model an object created via [HDCATE]
#' @param fit.treated function that accepts a data.frame as the only argument, fits the treated expectation function, and returns a fitted object
#' @param fit.untreated function that accepts a data.frame as the only argument, fits the untreated expectation function, and returns a fitted object
#' @param fit.propensity function that accepts a data.frame as the only argument, fits the propensity function, and return a fitted object
#' @param predict.treated function that accepts the returned object of `fit.treated` and a data.frame as arguments, and returns the predicted vector of that data.frame
#' @param predict.untreated function that accepts the returned object of `fit.untreated` and a data.frame as arguments, and returns the predicted vector that data.frame
#' @param predict.propensity function that accepts the returned object of `fit.propensity` and a data.frame as arguments, and returns the predicted vector that data.frame
#'
#' @return None.
#' @export
#'
#' @examples
#' # get simulation data
#' n_obs <- 500  # Num of observations
#' n_var <- 100  # Num of observed variables
#' n_rel_var <- 4  # Num of relevant variables
#' data <- HDCATE.get_sim_data(n_obs, n_var, n_rel_var)
#' # conditional expectation model is misspecified
#' x_formula <- paste(paste0('X', c(2:n_var)), collapse ='+')
#' # propensity score model is misspecified
#' # x_formula <- paste(paste0('X', c(1:(n_var-1))), collapse ='+')
#'
#' # create a new HDCATE model
#' model <- HDCATE(data=data, y_name='Y', d_name='D', x_formula=x_formula)
#'
#' # manually define a lasso method
#' my_lasso_fit_exp <- function(df) {
#'   hdm::rlasso(as.formula(paste0('Y', "~", x_formula)), df)
#' }
#' my_lasso_predict_exp <- function(fitted_model, df) {
#'   predict(fitted_model, df)
#' }
#' my_lasso_fit_ps <- function(df) {
#'   hdm::rlassologit(as.formula(paste0('D', "~", x_formula)), df)
#' }
#' my_lasso_predict_ps <- function(fitted_model, df) {
#'   predict(fitted_model, df, type="response")
#' }
#'
#' # Apply the "my-lasso" apporach to the first stage
#' HDCATE.set_first_stage(
#'   model,
#'   my_lasso_fit_exp,
#'   my_lasso_fit_exp,
#'   my_lasso_fit_ps,
#'   my_lasso_predict_exp,
#'   my_lasso_predict_exp,
#'   my_lasso_predict_ps
#' )
#'
HDCATE.set_first_stage <- function(model,
                                   fit.treated,
                                   fit.untreated,
                                   fit.propensity,
                                   predict.treated,
                                   predict.untreated,
                                   predict.propensity) {
  model$user_defined_first_stage <- TRUE
  model$propensity_est_method <- 'user-defined'
  model$expectation_est_method <- 'user-defined'
  model$user_defined_fit_treated <- fit.treated
  model$user_defined_fit_untreated <- fit.untreated
  model$user_defined_fit_ps <- fit.propensity
  model$user_defined_predict_treated <- predict.treated
  model$user_defined_predict_untreated <- predict.untreated
  model$user_defined_predict_ps <- predict.propensity
}

#' Clear the user-defined first-stage estimating methods
#'
#' @description Inverse operation of [HDCATE.set_first_stage]
#'
#' @param model an object created via [HDCATE]
#'
#' @return None.
#' @export
#'
#' @examples
#' # get simulation data
#' n_obs <- 500  # Num of observations
#' n_var <- 100  # Num of observed variables
#' n_rel_var <- 4  # Num of relevant variables
#' data <- HDCATE.get_sim_data(n_obs, n_var, n_rel_var)
#' # conditional expectation model is misspecified
#' x_formula <- paste(paste0('X', c(2:n_var)), collapse ='+')
#' # propensity score model is misspecified
#' # x_formula <- paste(paste0('X', c(1:(n_var-1))), collapse ='+')
#'
#' # create a new HDCATE model
#' model <- HDCATE(data=data, y_name='Y', d_name='D', x_formula=x_formula)
#'
#' # ... manually set user-defined first-stage estimating methods via `HDCATE.set_first_stage`
#'
#' # Clear those user-defined methods and use the built-in method
#' HDCATE.unset_first_stage(model)
HDCATE.unset_first_stage <- function(model) {
  model$user_defined_first_stage <- FALSE
  model$propensity_est_method <- 'LASSO'
  model$expectation_est_method <- 'LASSO'
  model$user_defined_fit_treated <- NA
  model$user_defined_fit_untreated <- NA
  model$user_defined_fit_ps <- NA
  model$user_defined_predict_treated <- NA
  model$user_defined_predict_untreated <- NA
  model$user_defined_predict_ps <- NA
}


#' Set bandwidth
#'
#' @description Set user-defined bandwidth.
#'
#' @param model an object created via [HDCATE]
#' @param bandwidth the value of bandwidth
#'
#' @return None.
#' @export
#'
#' @examples
#' # get simulation data
#' n_obs <- 500  # Num of observations
#' n_var <- 100  # Num of observed variables
#' n_rel_var <- 4  # Num of relevant variables
#' data <- HDCATE.get_sim_data(n_obs, n_var, n_rel_var)
#' # conditional expectation model is misspecified
#' x_formula <- paste(paste0('X', c(2:n_var)), collapse ='+')
#' # propensity score model is misspecified
#' # x_formula <- paste(paste0('X', c(1:(n_var-1))), collapse ='+')
#'
#' # create a new HDCATE model
#' model <- HDCATE(data=data, y_name='Y', d_name='D', x_formula=x_formula)
#'
#' # Set user-defined bandwidth, e.g., 0.15.
#' HDCATE.set_bw(model, 0.15)
HDCATE.set_bw <- function(model, bandwidth='default') {
  HDCATE.conf.set_bw(model, bandwidth)
}

# Set weight
# @keywords internal
# @description When fitting the local linear regression on the second stage, the sample weight can also be manually set  before calling [HDCATE.fit].
#
# @param model an object created via [HDCATE]
# @param weight the vector of sample weights
#
# @return None.
#
# @examples
# # get simulation data
# n_obs <- 500  # Num of observations
# n_var <- 100  # Num of observed variables
# n_rel_var <- 4  # Num of relevant variables
# data <- HDCATE.get_sim_data(n_obs, n_var, n_rel_var)
# # conditional expectation model is misspecified
# x_formula <- paste(paste0('X', c(2:n_var)), collapse ='+')
# # propensity score model is misspecified
# # x_formula <- paste(paste0('X', c(1:(n_var-1))), collapse ='+')
#
# # create a new HDCATE model
# model <- HDCATE(data=data, y_name='Y', d_name='D', x_formula=x_formula)
#
# # Set user-defined sample weights
# HDCATE.set_weight(model, c(rep(0.8, 250), rep(1.2, 250)))
# # To reset to default weight, call:
# HDCATE.set_weight(model)
HDCATE.set_weight <- function(model, weight='default') {
  HDCATE.conf.set_local_weight(model, weight)
}
