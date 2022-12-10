## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  warning = FALSE,
  message = FALSE
)

# save built-in output hook
hook_output = knitr::knit_hooks$get("output")

# set a new output hook to cut long output texts
knitr::knit_hooks$set(output = function(x, options) {
  if (!is.null(n <- options$out.lines)) {
    x = xfun::split_lines(x)
    if (length(x) > n) {
      # cut
      x = c(head(x, n), '....\n')
    }
    x = paste(x, collapse = '\n')
  }
  hook_output(x, options)
})

## ----installation, eval=FALSE-------------------------------------------------
#  install.packages('devtools')
#  devtools::install_github('microfan1/hdcate', build_vignettes = TRUE)

## ----load package-------------------------------------------------------------
library(hdcate)

## -----------------------------------------------------------------------------
set.seed(1)  # set an arbitrary random seed to ensure the reproducibility
data <- HDCATE.get_sim_data(n_obs=500, n_var=100, n_rel_var=4, intercept=10)  # generate data

## ----data, out.lines=6--------------------------------------------------------
data

## ----x_formula----------------------------------------------------------------
x_formula <- "X1+X2+X3+X4+X5+X6"
x_formula

## ----model--------------------------------------------------------------------
model <- HDCATE(data=data, y_name='Y', d_name='D', x_formula=x_formula)

## ----set conditional variable-------------------------------------------------
HDCATE.set_condition_var(model, 'X2', min=-1, max=1, step=0.01)

## ----fit----------------------------------------------------------------------
HDCATE.fit(model)

## ----get fitted ATE-----------------------------------------------------------
# ATE
model$ate

## ----get fitted CATE, out.lines=6---------------------------------------------
# CATE
model$hdcate

## ----inferences---------------------------------------------------------------
HDCATE.inference(model)

## ----two-sided test, out.lines=6----------------------------------------------
# two-sided bands
model$i_two_side

## ----left-sided test, out.lines=6---------------------------------------------
# left-sided bands
model$i_left_side

## ----right-sided test, out.lines=6--------------------------------------------
# right-sided bands
model$i_right_side

## ----plot, fig.show='hold' , fig.width=6, fig.height=6, fig.align='center', fig.cap='HDCATE Estimation (propensity function is misspecified)', dpi=200, out.width="80%"----
HDCATE.plot(model)

# Alternatively, if the model is not inferenced, use the code below:
# HDCATE.plot(model, include_band=FALSE)

# Alternatively, we may want to save the full plot as a high-definition PDF file:
# HDCATE.plot(model, output_pdf=TRUE, pdf_name='hdcate_plot.pdf', 
#   include_band=TRUE, test_side='both', y_axis_min='auto', y_axis_max='auto',
#   display.hdcate='HDCATEF', display.ate='ATE', display.siglevel='sig_level'
# )

## ----change formula, fig.show='hold' , fig.width=6, fig.height=6, fig.align='center', fig.cap='HDCATE Estimation (expectation function is misspecified)', dpi=200, out.width="80%"----
x_formula = paste(paste0('X', c(90:100)), collapse ='+')  # x_formula="X90+...+X100"

model <- HDCATE(data=data, y_name='Y', d_name='D', x_formula=x_formula)

HDCATE.set_condition_var(model, 'X2', min=-1, max=1, step=0.01)
HDCATE.fit(model)
HDCATE.inference(model)
HDCATE.plot(model)

## ----cross-fitting------------------------------------------------------------
# create the model first, or use an existing `model` and skip this line
model <- HDCATE(data=data, y_name='Y', d_name='D', x_formula=x_formula)

# suppose we want to use the 5-fold cross-fitting approach:
HDCATE.use_cross_fitting(model, k_fold=5)

## ----manual fold--------------------------------------------------------------
# In this case, the param k_fold is auto detected, you can pass any value to it.
HDCATE.use_cross_fitting(model, k_fold=1, folds=list(c(1:250), c(251:500)))

## ----change to full sample----------------------------------------------------
HDCATE.use_full_sample(model)

## ----temp function------------------------------------------------------------
my_method_to_fit_expectation <- function(df) {
  # fit the input data.frame (df), where the first column is outcome value (Y)
  #   and the rest are covariates in x_formula
  # ...
  # return a fitted object
}
my_method_to_predict_expectation <- function(fitted_model, df) {
  # use the fitted object (returned by your function above) to predict the conditional expectation
  #   fucntion on the input data.frame (df)
  # ...
  # return a vector of predicted values
}
my_method_to_fit_propensity <- function(df) {
  # fit the input data.frame (df), where the first column is dummy treatment 
  #   value (D) and the rest are covariates in x_formula
  # ...
  # return a fitted object
}
my_method_to_predict_propensity <- function(fitted_model, df) {
  # use the fitted object (returned by your function above) to predict the propensity score 
  #   fucntion on the input data.frame (df)
  # ...
  # return a vector of predicted values
}

## ----lasso example------------------------------------------------------------
x_formula <- "X1+X2+X3+X4+X5+X6"  # propensity function is misspecified
# create a model
model <- HDCATE(data=data, y_name='Y', d_name='D', x_formula=x_formula)

# using the template above, we can write:
my_method_to_fit_expectation <- function(df) {
  # fit the input data.frame (df), where the first column is outcome value (Y)
  #   and the rest are covariates in x_formula
  hdm::rlasso(as.formula(paste0('Y', "~", x_formula)), df)
  # return a fitted object
}
my_method_to_predict_expectation <- function(fitted_model, df) {
  # use the fitted object (fitted_model) to predict the conditional expectation
  #   fucntion on the input data.frame (df)
  predict(fitted_model, df)
  # return a vector of predicted values
}
my_method_to_fit_propensity <- function(df) {
  # fit the input data.frame (df), where the first column is dummy treatment 
  #   variable (D) and the rest are covariates in x_formula
  hdm::rlassologit(as.formula(paste0('D', "~", x_formula)), df)
  # return a fitted object
}
my_method_to_predict_propensity <- function(fitted_model, df) {
  # use the fitted object (fitted_model) to predict the propensity score 
  #   fucntion on the input data.frame (df)
  predict(fitted_model, df, type="response")
  # return a vector of predicted values
}

## ----reset1, include=FALSE----------------------------------------------------
model <- HDCATE(data=data, y_name='Y', d_name='D', x_formula=x_formula)

## ----plug in------------------------------------------------------------------
HDCATE.set_first_stage(
  model,
  fit.treated=my_method_to_fit_expectation,
  fit.untreated=my_method_to_fit_expectation,
  fit.propensity=my_method_to_fit_propensity,
  predict.treated=my_method_to_predict_expectation,
  predict.untreated=my_method_to_predict_expectation,
  predict.propensity=my_method_to_predict_propensity
)

## ----re estimate 1, fig.show='hold' , fig.width=6, fig.height=6, fig.align='center', fig.cap='HDCATE Estimation (using user-defined LASSO on the first stage)', dpi=200, out.width="80%"----
HDCATE.set_condition_var(model, 'X2', min=-1, max=1, step=0.01)
HDCATE.fit(model)
HDCATE.inference(model)
HDCATE.plot(model)

## ----recover------------------------------------------------------------------
HDCATE.unset_first_stage(model)

## ----reset before rf, include=FALSE-------------------------------------------
model <- HDCATE(data=data, y_name='Y', d_name='D', x_formula=x_formula)

## ----ramdom forest example, fig.show='hold' , fig.width=6, fig.height=6, fig.align='center', fig.cap='HDCATE Estimation (using user-defined Random Forest model on the first stage)', dpi=200, out.width="80%"----
# propensity function is misspecified, for example
x_formula <- "X1+X2+X3+X4+X5+X6"
# create a model
model <- HDCATE(data=data, y_name='Y', d_name='D', x_formula=x_formula)

# =========== User-defined ML approach start ===========

# using the template above, we can write:
my_method_to_fit_expectation <- function(df) {
  # fit the input data.frame (df), where the first column is outcome value (Y)
  #   and the rest are covariates in x_formula
  randomForest::randomForest(as.formula(paste0('Y', "~", x_formula)), data = df)
  # return a fitted object
}
my_method_to_predict_expectation <- function(fitted_model, df) {
  # use the fitted object (fitted_model) to predict the conditional expectation
  #   fucntion on the input data.frame (df)
  predict(fitted_model, df)
  # return a vector of predicted values
}
my_method_to_fit_propensity <- function(df) {
  # fit the input data.frame (df), where the first column is dummy treatment 
  #   variable (D) and the rest are covariates in x_formula
  randomForest::randomForest(as.formula(paste0('D', "~", x_formula)), data = df)
  # return a fitted object
}
my_method_to_predict_propensity <- function(fitted_model, df) {
  # use the fitted object (fitted_model) to predict the propensity score 
  #   fucntion on the input data.frame (df)
  predict(fitted_model, df)
  # return a vector of predicted values
}

HDCATE.set_first_stage(
  model,
  fit.treated=my_method_to_fit_expectation,
  fit.untreated=my_method_to_fit_expectation,
  fit.propensity=my_method_to_fit_propensity,
  predict.treated=my_method_to_predict_expectation,
  predict.untreated=my_method_to_predict_expectation,
  predict.propensity=my_method_to_predict_propensity
)

# =========== User-defined ML approach end ===========

# set conditional variable in CATE, same as above
HDCATE.set_condition_var(model, 'X2', min=-1, max=1, step=0.01)

# fit, inference and plot
HDCATE.fit(model)
HDCATE.inference(model)
HDCATE.plot(model)

## ----var importance in rf, fig.show='hold' , fig.width=6, fig.height=6, fig.align='center', fig.cap='Variable Importance Extracted from the User-defined Random Forest Model'----
# the list have only one element since we use full-sample estimator
rf_model <- model$treated_cond_exp_model_list[[1]]

# then, we can use the fitted object `rf_model` for futher analysis.
# the rest of codes are ommited since they are unconcerned for this package.
# ...

## ----plot rf, echo=FALSE, fig.align='center', fig.cap='Variable Importance Extracted from the User-defined Random Forest Model', fig.height=6, fig.width=6, dpi=200, fig.show='hold', out.width="80%"----
# the code below uses randomForest and ggplot2 to plot out variable importance
library('ggplot2') # visualization
library('ggthemes') # visualization
library('dplyr') # data manipulation
importance    <- randomForest::importance(rf_model)
# varImportance <- data.frame(Variables = row.names(importance), Importance = round(importance[ ,'MeanDecreaseGini'],2))  # for classification problem
varImportance <- data.frame(Variables = row.names(importance), Importance = round(importance[ ,'IncNodePurity'],2))  # for regression problem

rankImportance <- varImportance %>%
  mutate(Rank = paste0('#',dense_rank(desc(Importance))))

ggplot(rankImportance, aes(x = reorder(Variables, Importance), 
    y = Importance, fill = Importance)) +
  geom_bar(stat='identity') + 
  # geom_text(aes(x = Variables, y = 0.5, label = Rank),
  #   hjust=0, vjust=0.55, size = 4, colour = 'red') +
  labs(x = 'Variables', y= 'Importance (measured in squared OOB residuals)') +
  coord_flip() + 
  theme_few()

## ----user-defined inf---------------------------------------------------------
HDCATE.inference(model, sig_level=0.05, n_rep_boot=3000)

## ----bw-----------------------------------------------------------------------
HDCATE.set_bw(model, 0.15)  # for example, set bandwidth=0.15

## ----recover bw---------------------------------------------------------------
HDCATE.set_bw(model)

## ----help, eval=FALSE---------------------------------------------------------
#  help(HDCATE)
#  help(HDCATE.set_condition_var)
#  help(HDCATE.fit)
#  help(HDCATE.inference)
#  help(HDCATE.plot)
#  help(HDCATE.use_cross_fitting)
#  help(HDCATE.use_full_sample)
#  help(HDCATE.set_first_stage)
#  help(HDCATE.unset_first_stage)
#  help(HDCATE.set_bw)

