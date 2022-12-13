#' Get simulation data
#'
#' @param n_obs Num of observations
#' @param n_var Num of covariates
#' @param n_rel_var Num of relevant variables, only the first `n_rel_var` covariates are actually present in the expectation function of potential outcome, and only the last `n_rel_var` covariates are present in the propensity score function.
#' @param sig_strength_propensity signal strength in propensity score functions
#' @param sig_strength_outcome signal strength in outcome functions
#' @param intercept value of intercept in outcome functions
#'
#' @return a data.frame, which is the simulated observed data.
#' @export
#'
#' @examples
#' HDCATE.get_sim_data()
#' HDCATE.get_sim_data(n_obs=50, n_var=4, n_rel_var=2)
HDCATE.get_sim_data <- function(n_obs=500, n_var=100, n_rel_var=4,
                         sig_strength_propensity=0.5,
                         sig_strength_outcome=1,
                         intercept=10){
  X <- matrix(rnorm(n_obs*n_var, mean=0, sd=1), nrow=n_obs, ncol = n_var)
  epsilon <- rnorm(n_obs, mean =0, sd = 1)
  indexp <- rowSums(X[, (n_var-n_rel_var+1):n_var])*sig_strength_propensity
  p <- exp(indexp)/ (1+exp(indexp))
  u <- runif(n_obs, min=0, max =1)
  D <- ifelse(p>u, 1, 0)
  y1 <- intercept + rowSums(X[, 1:n_rel_var]) * sig_strength_outcome + epsilon
  Y <- D*y1
  message(paste0('Actual CATE function is: CATE(x)=', intercept, '+', sig_strength_outcome, '*x'))
  data.frame(Y,D,X)
}
