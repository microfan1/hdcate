#' @title High-Dimensional Conditional Average Treatment Effects (HDCATE) Estimator
#'
#' @keywords internal
#'
#' @description
#' Use a two-step procedure to estimate the conditional average treatment effects (CATE) for all possible values of the covariate(s).
#'
#'
#' @format [R6::R6Class] object.
#'
#' @family HDCATE
HDCATE_R6Class <- R6Class(
  # different instances for different models
  "HDCATE",
  public=list(
    data=NA,
    y_name=NA,
    d_name=NA,
    x_formula=NA,
    cond_var_name=NA,
    cond_var_lower=NA,
    cond_var_upper=NA,
    cond_var_eval_step=0.01,
    use_cross_fitting=FALSE,
    user_defined_first_stage=FALSE,
    user_defined_fit_untreated=NA,
    user_defined_fit_treated=NA,
    user_defined_fit_ps=NA,
    user_defined_predict_untreated=NA,
    user_defined_predict_treated=NA,
    user_defined_predict_ps=NA,
    propensity_est_method='lasso',
    expectation_est_method='lasso',
    propensity_score_model_list=list(),
    treated_cond_exp_model_list=list(),
    untreated_cond_exp_model_list=list(),
    local_reg_model_list=list(),
    bw='rule-of-thumb',
    k_fold=NA,
    folds=NULL,
    local_weight=NA,
    hdcate=NA,
    sigma_hat=NA,
    eta_hat_func_list=NA,
    eta_hat_list=NA,
    d=NA,
    h=NA,
    n_obs=NA,
    alpha=NA,
    n_alpha=NA,
    i_two_side=NA,
    i_left_side=NA,
    i_right_side=NA,
    cond_interval=NA,
    length_grid=NA,
    ate=NA,
    vate=NA,
    initialize = function(data, y_name, d_name, x_formula){
      # necessary param
      self$data <- data
      self$y_name <- y_name
      self$d_name <- d_name
      self$x_formula <- x_formula

      # init model
      # equal weight
      self$local_weight <- rep(1, nrow(data))
      self$n_obs = nrow(data)
    },
    propensity_hd_estimate = function(data=NA, verbose=F) {
      data <- private$check_data(data)
      method <- self$propensity_est_method
      d_name <- self$d_name
      x_formula <- self$x_formula

      if (verbose) {message(paste0('Start estimating model for propensity score, method=', method))}

      if (!self$user_defined_first_stage){
        if (method == "lasso") {
          model <- hdm::rlassologit(as.formula(paste0(d_name, "~", x_formula)), data=data)
        } else {
          stop('When estimating propensity score, a valid method should be provided, such as "lasso".')
        }
      } else {
        model <- self$user_defined_fit_ps(model.frame(as.formula(paste0(d_name, "~", x_formula)), data))
      }

      if (verbose) {message(paste0('Finish estimating propensity score.'))}

      self$propensity_score_model_list <- append(self$propensity_score_model_list, list(model))
      return(model)
    },
    conditional_expectations_hd_estimate = function(data=NA, verbose=F) {
      data <- private$check_data(data)
      method <- self$expectation_est_method
      y_name <- self$y_name
      x_formula <- self$x_formula

      if (verbose) {message(paste0('Start estimating conditional expectations, method=', method))}
      filter_treated<-data[,self$d_name]==1
      filter_untreated<-data[,self$d_name]==0

      if (!self$user_defined_first_stage){
        if (method == 'lasso') {
          model_untreated<-hdm::rlasso(as.formula(paste0(y_name, "~", x_formula)),
                                       data=data[filter_untreated,])
          model_treated<-hdm::rlasso(as.formula(paste0(y_name, "~", x_formula)),
                                     data=data[filter_treated,])
        } else {
          stop('When estimating conditional expectations, a valid method should be provided, such as "lasso".')
        }
      } else {
        model_untreated<-self$user_defined_fit_untreated(model.frame(as.formula(paste0(y_name, "~", x_formula)), data[filter_untreated,]))
        model_treated<-self$user_defined_fit_treated(model.frame(as.formula(paste0(y_name, "~", x_formula)), data[filter_treated,]))
      }

      if (verbose) {message(paste0('Finish estimating conditional expectations.'))}


      self$untreated_cond_exp_model_list <- append(self$untreated_cond_exp_model_list, list(model_untreated))
      self$treated_cond_exp_model_list <- append(self$treated_cond_exp_model_list, list(model_treated))

      return(list(predictor_untreated=model_untreated, predictor_treated=model_treated))
    },
    first_stage = function(data=NA, verbose=F) {
      data <- private$check_data(data)
      predictor_propensity <- list(predictor_propensity=self$propensity_hd_estimate(data, verbose=verbose))
      predictor_cond_expectations <- self$conditional_expectations_hd_estimate(data, verbose=verbose)
      predictor_eta_hat <- append(predictor_cond_expectations, predictor_propensity)
      return(predictor_eta_hat)
    },
    second_stage = function(predictor_eta_hat=NA, eta_hat=NA, subsample_idx=NULL, local_weight=NULL, estimate_std=TRUE, verbose=FALSE, save_model=TRUE) {
      y_name <- self$y_name
      d_name <- self$d_name
      x_formula <- self$x_formula

      if (verbose) {
        if (self$use_cross_fitting) {
          message(paste0('Start estimating HDCATE in a fold.'))
        } else {
          message(paste0('Start estimating HDCATE.'))
        }
      }

      if (is.null(subsample_idx) || anyNA(subsample_idx)) {
        # full-sample
        subsample_idx <- c(1:nrow(self$data))
      }
      if (is.null(local_weight) || anyNA(local_weight)) {
        # equal weight (out of bootstrap)
        local_weight <- self$local_weight[subsample_idx]
      }

      # Compute phi
      if (is.null(eta_hat) || anyNA(eta_hat)) {
        # for fitting case
        if (!self$user_defined_first_stage) {
          mu0_hat <- predict(predictor_eta_hat$predictor_untreated, self$data[subsample_idx,])
          mu1_hat <- predict(predictor_eta_hat$predictor_treated, self$data[subsample_idx,])
          pi_hat <- predict(predictor_eta_hat$predictor_propensity, self$data[subsample_idx,], type="response")
        } else {
          mu0_hat <- self$user_defined_predict_untreated(predictor_eta_hat$predictor_untreated, self$data[subsample_idx,])
          mu1_hat <- self$user_defined_predict_treated(predictor_eta_hat$predictor_treated, self$data[subsample_idx,])
          pi_hat <- self$user_defined_predict_ps(predictor_eta_hat$predictor_propensity, self$data[subsample_idx,])
        }
        eta_hat <- list(untreated=mu0_hat, treated=mu1_hat, propensity=pi_hat)
      } else {
        # for inference case
        mu0_hat<-eta_hat$untreated
        mu1_hat<-eta_hat$treated
        pi_hat<- eta_hat$propensity
      }
      D <- self$data[subsample_idx, self$d_name]
      Y <- self$data[subsample_idx, self$y_name]
      phi <- D * (Y - mu1_hat) / pi_hat + mu1_hat - (1 - D) * (Y - mu0_hat) / (1 - pi_hat) - mu0_hat

      # Local linear regression
      # Gaussian kernel
      ker <- locpol::gaussK
      cond_var <- self$data[subsample_idx, self$cond_var_name]
      cond_interval <- self$cond_interval

      # ATE
      ate <- mean(phi)
      vate <- ate*rep(1,length(cond_interval))

      # bandwidth
      h <- self$get_bw(
        phi=phi,
        use_sample_idx=use_sample_idx
      )

      model <- locpol::locPolSmootherC(x=cond_var, y=phi, xeval=cond_interval, bw=h,
                               weig=local_weight, deg=1, kernel=ker)
      if (save_model) {
        self$local_reg_model_list <- append(self$local_reg_model_list, list(model))
      }

      # get HDCATE estimator, i,e. tau_hat in the paper
      hdcate <- model$beta0

      if (verbose) {
        message(paste0('Finish estimating HDCATE.'))
      }

      if (estimate_std) {
        if (verbose) {
          message(paste0('Start estimating standard errors.'))
        }

        length_grid <- self$length_grid
        fX_hat <- numeric(length_grid)
        sigmasq_hat <- numeric(length_grid)
        d <- self$d
        for(i in 1:length_grid)
        {
          fX_hat[i]<-mean(ker((cond_var-cond_interval[i])/h))/h
          sigmasq_hat[i]<-mean((phi-hdcate[i])^2*(ker((cond_var-cond_interval[i])/h))^2)/(fX_hat[i]^2)/(h^d)
        }
        sigma_hat <- sqrt(sigmasq_hat)

        if (verbose) {
          message(paste0('Finish estimating standard errors.'))
        }
      } else {
        sigma_hat <- NULL
      }

      self$h <- h
      return(list(hdcate=hdcate, sigma_hat=sigma_hat, vate=vate, ate=ate, eta_hat=eta_hat))
    },
    get_bw = function(phi, use_sample_idx) {
      if (self$bw == 'plug-in'){
        # use subsample to choose bandwidth
        use_full_sample <- FALSE
      } else {
        # use full sample to choose bandwidth
        use_full_sample <- TRUE
        # update index
        use_sample_idx <- c(1:nrow(data))
      }

      # bandwidth
      cond_var <- self$data[use_sample_idx, self$cond_var_name]
      sample_size <- nrow(self$data[use_sample_idx, ])
      if (self$bw == 'rule-of-thumb'){
        # Silverman (1986), rule-of-thumb
        h_hat <- stats::bw.nrd0(cond_var)
        h <- h_hat * sample_size^(1/5) * sample_size^(-2/7)
      } else if (self$bw == 'plug-in'){
        # Ruppert,Sheather and Wand (1995), the plug-in method
        h_hat <- KernSmooth::dpill(x=cond_var, y=phi)
        h <- h_hat * sample_size^(1/5) * sample_size^(-2/7)
      } else if (is.numeric(self$bw)){
        # bandwidth can be specified by user
        h <- self$bw
      }

      # avoid h=NaN
      if (anyNA(h) || is.null(h)) {
        h <- 0.01
      }
      return(h)
    },
    #' @description Fit the HDCATE function
    #'
    #' @return estimated HDCATE
    fit = function(verbose=FALSE){
      if (!self$use_cross_fitting){
        # Full-sample estimator
        message(paste0('Use full-sample estimator.'))
        # First stage
        predictor_eta_hat <- self$first_stage(verbose=verbose)
        # return(predictor_eta_hat)
        # Second stage
        res_second_stage <- self$second_stage(predictor_eta_hat=predictor_eta_hat, verbose=verbose)

        self$hdcate <- res_second_stage$hdcate
        self$ate <- res_second_stage$ate
        self$vate <- res_second_stage$vate
        self$sigma_hat <- res_second_stage$sigma_hat
        self$eta_hat_list <- list(full_sample=res_second_stage$eta_hat)
      } else {
        if (is.numeric(self$k_fold) && (self$k_fold %% 1 == 0) && (self$k_fold >= 2)) {
          # K-fold cross-fitting estimator
          message(paste0('Use ', self$k_fold, '-fold cross-fitting estimator.'))
          if (is.null(self$folds) || anyNA(self$folds)) {
            self$folds <- caret::createFolds(self$data[, self$y_name], k=self$k_fold)
          }
          cv_res <- lapply(
            self$folds,
            function(idx) {
              data_kc <- self$data[-idx, ]
              # First stage
              predictor_eta_hat_kc <- self$first_stage(data=data_kc, verbose=verbose)
              # Second stage
              res <- self$second_stage(
                predictor_eta_hat=predictor_eta_hat_kc,
                subsample_idx=idx,
                verbose=verbose
              )
              return(res)
            }
          )

          self$hdcate <- as.matrix(rowMeans(sapply(cv_res, function(obj) obj$hdcate)))
          self$ate <- mean(sapply(cv_res, function(obj) obj$ate))
          self$vate <- as.matrix(rowMeans(sapply(cv_res, function(obj) obj$vate)))
          self$sigma_hat <- as.matrix(rowMeans(sapply(cv_res, function(obj) obj$sigma_hat)))
          self$eta_hat_list <- lapply(cv_res, function(obj) obj$eta_hat)
        } else {
          stop(paste0('`k_fold`=', self$k_fold, ' is not supported for cross-fitting method'))
        }
        # return(list(cv_res=cv_res))
      }
    },
    inference = function(sig_level=0.01, boot_method='normal', n_rep_boot=1000, verbose=FALSE){
      # uniform confidence bands
      message(paste0('Constructing uniform confidence bands...'))

      # validation
      v1 <- !(is.null(self$hdcate) || anyNA(self$hdcate))
      v2 <- !(is.null(self$sigma_hat) || anyNA(self$sigma_hat))
      v3 <- !(is.null(self$eta_hat_list) || anyNA(self$eta_hat_list))
      if (!(v1 && v2 && v3)) {
        stop('Missing fitted models. You should call `HDCATE.fit()` before calling `HDCATE.inference()`.')
      }

      # bootstrap
      # n_rep_boot<-2  # for test, just delete this line
      weights = self$draw_weights(boot_method, n_rep_boot, self$n_obs)  # one row for one boot
      hdcate_stat_boot_one_sided <- numeric(n_rep_boot)
      hdcate_stat_boot_two_sided <- numeric(n_rep_boot)
      for (b in 1:n_rep_boot){
        if (is.null(self$folds) || anyNA(self$folds)) {
          # full_sample
          boot_loc_res <- self$second_stage(
            eta_hat=self$eta_hat_list[[1]],
            subsample_idx=NULL,
            estimate_std=FALSE,
            local_weight=weights[,b],
            verbose=verbose,
            save_model=FALSE
          )
          hdcate_boot <- boot_loc_res$hdcate
        } else {
          cv_res <- mapply(
            function(idx, list_idx) {
              boot_loc_res <- self$second_stage(
                eta_hat=self$eta_hat_list[[list_idx]],
                subsample_idx=idx,
                estimate_std=FALSE,
                local_weight=weights[idx,b],
                verbose=verbose,
                save_model=FALSE
              )
              return(boot_loc_res)
            },
            self$folds,
            seq_along(self$folds),
            SIMPLIFY=FALSE
          )
          hdcate_boot <- as.matrix(rowMeans(sapply(cv_res, function(obj) obj$hdcate)))
        }
        hdcate_stat_boot_one_sided[b] <- max(sqrt(self$n_obs * (self$h ^ self$d)) * (hdcate_boot - self$hdcate) / self$sigma_hat)
        hdcate_stat_boot_two_sided[b] <- max(sqrt(self$n_obs * (self$h ^ self$d)) * abs(hdcate_boot - self$hdcate) / self$sigma_hat)
      }
      # return(list(hdcate_stat_boot_one_sided, hdcate_stat_boot_two_sided, hdcate_boot, weights))
      # CI
      # alpha can be a vector like c(0.10, 0.05, 0.01)
      alpha <- sig_level
      n_alpha <- length(alpha)
      c_alpha_one_sided <- quantile(hdcate_stat_boot_one_sided, 1 - alpha, na.rm = TRUE)
      c_alpha_two_sided <- quantile(hdcate_stat_boot_two_sided, 1 - alpha, na.rm = TRUE)

      # Uniform confidence bands
      i_left_side <- matrix(0, self$length_grid, n_alpha,
                            dimnames = list(paste0('CATE on ', self$cond_var_name, '=', self$cond_interval), paste0('Sigificant level ', alpha)))
      i_right_side <- matrix(0, self$length_grid, n_alpha,
                             dimnames = list(paste0('CATE on ', self$cond_var_name, '=', self$cond_interval), paste0('Sigificant level ', alpha)))
      i_two_side <- matrix(0, self$length_grid, n_alpha * 2,
                           dimnames = list(
                             paste0('CATE on ', self$cond_var_name, '=', self$cond_interval),
                             rep(paste0('Sigificant level ', alpha), 2))
      )
      for (idx in 1:n_alpha){
        i_left_side[,idx] <- self$hdcate - c_alpha_one_sided[idx] * self$sigma_hat / sqrt(self$n_obs * (self$h ^ self$d))
        i_right_side[,idx] <- self$hdcate + c_alpha_one_sided[idx] * self$sigma_hat / sqrt(self$n_obs * (self$h ^ self$d))
        i_two_side[,idx] <- i_left_side[,idx]
        i_two_side[,idx+n_alpha] <- i_right_side[,idx]
      }

      self$alpha <- alpha
      self$n_alpha <- n_alpha
      self$i_two_side <- i_two_side
      self$i_left_side <- i_left_side
      self$i_right_side <- i_right_side
      message(paste0('Uniform confidence bands are constructed.'))
      # res <- list(i_two_side, i_left_side, i_right_side, hdcate, vate)
      # return(list(cv_res=cv_res))
    },
    #' @description Plot the results.
    #'
    #' @param output_pdf if `TRUE`, save image to a pdf file named as `pdf_name`
    #' @param pdf_name the name of the output PDF file
    #' @param include_band if `TRUE`, plot uniform confidence bands as well.
    #' @param test_side `'both'` for a 2-sided test, `'left'` for a left-sided test or `'right'` for a right-sided test
    #' @param y_axis_min the lowest value plotted in Y axis, the default is `'auto'`
    #' @param y_axis_max the largest value plotted in Y axis, the default is `'auto'`
    #' @param display.hdcate the name of HDCATE function in the legend, the default is 'HDCATEF'
    #' @param display.ate the name of average treatment effect in the legend, the default is 'ATE'
    #' @param display.siglevel the name of the significant level for confidence bands in the legend, the default is 'sig_level'
    plot=function(output_pdf=FALSE, pdf_name='hdcate_plot.pdf', include_band=TRUE, test_side='both', y_axis_min='auto', y_axis_max='auto', display.hdcate='HDCATEF', display.ate='ATE', display.siglevel='sig_level'){
      min_value <- min(self$hdcate, self$vate)
      max_value <- max(self$hdcate, self$vate)

      if (is.null(self$i_two_side) || anyNA(self$i_two_side)) {
        include_band <- FALSE
      }

      if (include_band) {
        cb_matrix <- self$get_confidence_bands(test_side)
        min_value <- min(min_value, cb_matrix)
        max_value <- max(max_value, cb_matrix)
      }

      # set y axis
      distance <- max_value - min_value
      if (y_axis_min=='auto'){
        ylower <- min_value - distance * 0.25
      }
      if (y_axis_max=='auto'){
        yupper <- max_value + distance * 0.25
      }

      # reset user's settings when the function is exited
      oldpar <- par(no.readonly = TRUE)
      on.exit(par(oldpar))

      # plot and save
      if (output_pdf) {
        pdf(pdf_name)
      }

      plot(self$cond_interval, self$hdcate, ylim=c(ylower,yupper), xlab = self$cond_var_name, ylab="", type="l", lty=1, lwd = 2)
      par(new = T)
      plot(self$cond_interval, self$vate, ylim=c(ylower,yupper), xlab ="",ylab="", type="l", lty=4, lwd = 1.25)
      lty_set <- c(1, 4)
      lwd_set <- c(2, 1.25)
      if (include_band) {
        for (idx in 1:self$n_alpha){
          # 2, 3, 5, 6
          idx_lty <- c(3, 2, 5, 6)[idx %% 4 + 1]
          # if (idx_lty == 0) {
          #   idx_lty = 6
          # }
          if (test_side=='left' || test_side=='right') {
            par(new = T)
            plot(self$cond_interval, cb_matrix[,idx], ylim=c(ylower,yupper), xlab ="", ylab="", type="l", lty=idx_lty)
            lty_set <- append(lty_set, idx_lty)
            lwd_set <- append(lwd_set, 1)
          } else {
            # two-side test
            par(new = T)
            plot(self$cond_interval, cb_matrix[,idx], ylim=c(ylower,yupper), xlab ="", ylab="", type="l", lty=idx_lty)
            par(new = T)
            plot(self$cond_interval, cb_matrix[,idx+self$n_alpha], ylim=c(ylower,yupper), xlab ="", ylab="", type="l", lty=idx_lty)

            lty_set <- append(lty_set, rep(idx_lty, 2))
            lwd_set <- append(lwd_set, rep(1,2))
          }
        }
        legend("topright", c(display.hdcate, display.ate, paste0(display.siglevel, '=', self$alpha)), lty=lty_set, lwd=lwd_set, cex=1.2)
      } else {
        legend("topright", c(display.hdcate, display.ate), lty=lty_set, lwd=lwd_set, cex=1.2)
      }
      if (output_pdf) {
        dev.off()
      }
    },
    get_confidence_bands=function(test_side='both') {
      if ((is.null(self$i_two_side) || anyNA(self$i_two_side)) && check_exist) {
        stop('Missing results. Make sure you have called `HDCATE.fit()` and `HDCATE.inference()` before calling this method.')
      }
      if (test_side=='left'){
        cb_matrix <- self$i_left_side
      } else if (test_side=='right') {
        cb_matrix <- self$i_right_side
      } else {
        # 2-side test
        cb_matrix <- self$i_two_side
      }
      return(cb_matrix)
    },
    draw_weights=function(method, n_rep_boot, n_obs){
      if (method == "bayes") {
        weights = stats::rexp(n_rep_boot * n_obs, rate = 1)
      } else if (method == "normal") {
        weights = stats::rnorm(n_rep_boot * n_obs, mean = 1, sd = 1)
      } else if (method == "wild") {
        weights = stats::rnorm(n_rep_boot * n_obs) / sqrt(2) +
          (stats::rnorm(n_rep_boot * n_obs)^2 - 1) / 2 + 1
      } else {
        stop("invalid boot method")
      }
      weights = matrix(weights, nrow = n_rep_boot, ncol = n_obs, byrow = TRUE)
      weights = t(weights)
      return(weights)
    },
    set_condition_var=function(name=NA, min=NA, max=NA, step=NA) {
      if (!(is.null(name) || anyNA(name))) {
        self$cond_var_name <- name
        self$d <- length(name)
      }
      if (!(is.null(min) || anyNA(min))) {
        self$cond_var_lower <- min
      }
      if (!(is.null(max) || anyNA(max))) {
        self$cond_var_upper <- max
      }
      if (!(is.null(step) || anyNA(step))) {
        self$cond_var_eval_step <- step
      }

      # validation
      v1 <- (!(is.null(self$cond_var_name) || anyNA(self$cond_var_name)))
      v2 <- (!(is.null(self$cond_var_lower) || anyNA(self$cond_var_lower)))
      v3 <- (!(is.null(self$cond_var_upper) || anyNA(self$cond_var_upper)))
      # step=0.01 by default (defined in hdcate R6 class)
      v4 <- (!(is.null(self$cond_var_eval_step) || anyNA(self$cond_var_eval_step)))
      if (!(v1 && v2 && v3)) {
        warning(paste0('Failed to update conditional variable to: name=', name, ', min=', min, ', max=', max, ', step=', step, '.'))
        stop('Missing some of these params for the conditional variable: `name`, `min` and `max`.')
      } else {
        self$cond_interval <- seq(self$cond_var_lower, self$cond_var_upper, self$cond_var_eval_step)
        self$length_grid <- length(self$cond_interval)
        message(paste0('Updated conditional variable to: name=', name, ', min=', min, ', max=', max, ', step=', step, '.'))
      }
    }
  ),
  private = list(
    check_data = function(data) {
      if (!(is.data.frame(data))) {
        # for full-sample estimator
        return(self$data)
      } else {
        # for split-sample estimator
        return(data)
      }
    }
  ),
)
