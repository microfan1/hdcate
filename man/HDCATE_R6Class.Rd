% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hdcate_r6_class.R
\name{HDCATE_R6Class}
\alias{HDCATE_R6Class}
\title{High-Dimensional Conditional Average Treatment Effects (HDCATE) Estimator}
\format{
\link[R6:R6Class]{R6::R6Class} object.
}
\description{
Use a two-step procedure to estimate the conditional average treatment effects (CATE) for all possible values of the covariate(s).
}
\concept{HDCATE}
\keyword{internal}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{HDCATE_R6Class$new()}}
\item \href{#method-propensity_hd_estimate}{\code{HDCATE_R6Class$propensity_hd_estimate()}}
\item \href{#method-conditional_expectations_hd_estimate}{\code{HDCATE_R6Class$conditional_expectations_hd_estimate()}}
\item \href{#method-first_stage}{\code{HDCATE_R6Class$first_stage()}}
\item \href{#method-second_stage}{\code{HDCATE_R6Class$second_stage()}}
\item \href{#method-get_bw}{\code{HDCATE_R6Class$get_bw()}}
\item \href{#method-fit}{\code{HDCATE_R6Class$fit()}}
\item \href{#method-inference}{\code{HDCATE_R6Class$inference()}}
\item \href{#method-plot}{\code{HDCATE_R6Class$plot()}}
\item \href{#method-get_confidence_bands}{\code{HDCATE_R6Class$get_confidence_bands()}}
\item \href{#method-draw_weights}{\code{HDCATE_R6Class$draw_weights()}}
\item \href{#method-set_condition_var}{\code{HDCATE_R6Class$set_condition_var()}}
\item \href{#method-clone}{\code{HDCATE_R6Class$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{HDCATE_R6Class$new(data, y_name, d_name, x_formula)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-propensity_hd_estimate"></a>}}
\if{latex}{\out{\hypertarget{method-propensity_hd_estimate}{}}}
\subsection{Method \code{propensity_hd_estimate()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{HDCATE_R6Class$propensity_hd_estimate(data = NA)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-conditional_expectations_hd_estimate"></a>}}
\if{latex}{\out{\hypertarget{method-conditional_expectations_hd_estimate}{}}}
\subsection{Method \code{conditional_expectations_hd_estimate()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{HDCATE_R6Class$conditional_expectations_hd_estimate(data = NA)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-first_stage"></a>}}
\if{latex}{\out{\hypertarget{method-first_stage}{}}}
\subsection{Method \code{first_stage()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{HDCATE_R6Class$first_stage(data = NA)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-second_stage"></a>}}
\if{latex}{\out{\hypertarget{method-second_stage}{}}}
\subsection{Method \code{second_stage()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{HDCATE_R6Class$second_stage(
  predictor_eta_hat = NA,
  eta_hat = NA,
  subsample_idx = NULL,
  local_weight = NULL,
  estimate_std = TRUE,
  quietly = FALSE,
  save_model = TRUE
)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-get_bw"></a>}}
\if{latex}{\out{\hypertarget{method-get_bw}{}}}
\subsection{Method \code{get_bw()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{HDCATE_R6Class$get_bw(phi, use_sample_idx)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-fit"></a>}}
\if{latex}{\out{\hypertarget{method-fit}{}}}
\subsection{Method \code{fit()}}{
Fit the HDCATE function
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{HDCATE_R6Class$fit()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
estimated HDCATE
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-inference"></a>}}
\if{latex}{\out{\hypertarget{method-inference}{}}}
\subsection{Method \code{inference()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{HDCATE_R6Class$inference(
  sig_level = 0.01,
  boot_method = "normal",
  n_rep_boot = 1000
)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-plot"></a>}}
\if{latex}{\out{\hypertarget{method-plot}{}}}
\subsection{Method \code{plot()}}{
Plot the results.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{HDCATE_R6Class$plot(
  output_pdf = FALSE,
  pdf_name = "hdcate_plot.pdf",
  include_band = TRUE,
  test_side = "both",
  y_axis_min = "auto",
  y_axis_max = "auto",
  display.hdcate = "HDCATEF",
  display.ate = "ATE",
  display.siglevel = "sig_level"
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{output_pdf}}{if \code{TRUE}, save image to a pdf file named as \code{pdf_name}}

\item{\code{pdf_name}}{the name of the output PDF file}

\item{\code{include_band}}{if \code{TRUE}, plot uniform confidence bands as well.}

\item{\code{test_side}}{\code{'both'} for a 2-sided test, \code{'left'} for a left-sided test or \code{'right'} for a right-sided test}

\item{\code{y_axis_min}}{the lowest value plotted in Y axis, the default is \code{'auto'}}

\item{\code{y_axis_max}}{the largest value plotted in Y axis, the default is \code{'auto'}}

\item{\code{display.hdcate}}{the name of HDCATE function in the legend, the default is 'HDCATEF'}

\item{\code{display.ate}}{the name of average treatment effect in the legend, the default is 'ATE'}

\item{\code{display.siglevel}}{the name of the significant level for confidence bands in the legend, the default is 'sig_level'}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-get_confidence_bands"></a>}}
\if{latex}{\out{\hypertarget{method-get_confidence_bands}{}}}
\subsection{Method \code{get_confidence_bands()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{HDCATE_R6Class$get_confidence_bands(test_side = "both")}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-draw_weights"></a>}}
\if{latex}{\out{\hypertarget{method-draw_weights}{}}}
\subsection{Method \code{draw_weights()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{HDCATE_R6Class$draw_weights(method, n_rep_boot, n_obs)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-set_condition_var"></a>}}
\if{latex}{\out{\hypertarget{method-set_condition_var}{}}}
\subsection{Method \code{set_condition_var()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{HDCATE_R6Class$set_condition_var(name = NA, min = NA, max = NA, step = NA)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{HDCATE_R6Class$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}