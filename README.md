# hdcate
Estimation of Conditional Average Treatment Effects With High-Dimensional Data

This package uses a two-step procedure to estimate the conditional average treatment effects (CATE) with potentially high-dimensional covariate(s). In the first stage, the nuisance functions necessary for identifying CATE can be estimated by machine learning methods, allowing the number of covariates to be comparable to or larger than the sample size. The second stage consists of a low-dimensional local linear regression, reducing CATE to a function of the covariate(s) of interest.

The CATE estimator implemented in this package not only allows for high-dimensional data, but also has the “double robustness” property: either the model for the propensity score or the models for the conditional means of the potential outcomes are allowed to be misspecified (but not both).

The algorithm used in this package is described in the paper by Fan et al., “Estimation of Conditional Average Treatment Effects With High-Dimensional Data” (2022), Journal of Business & Economic Statistics. (doi:10.1080/07350015.2020.1811102)

This package can be downloaded/installed from CRAN, "hdcate", check CRAN for updates also.

## User Manual

User manual is provided at <https://microfan1.github.io/hdcate/>

## Installation (Github Version)

```R
devtools::install_github('microfan1/hdcate')
```
