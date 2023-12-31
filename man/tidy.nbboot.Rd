% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tidy_boots.R
\name{tidy.nbboot}
\alias{tidy.nbboot}
\title{Tidy function for nbboot}
\usage{
\method{tidy}{nbboot}(
  x,
  coef_type = c("original", "bootstrap_mean", "bootstrap_median"),
  standard_error = TRUE,
  p_value = c("bootstrapped", "approx", "none"),
  confint = c("none", "bootstrapped", "approx"),
  conf_level = 0.95,
  approx_t_value = TRUE,
  symmetric_bootstrap_p = TRUE,
  include_ylev = FALSE,
  ...
)
}
\arguments{
\item{x}{A fitted bootstrapped zero-inflated model}

\item{coef_type}{What type of coefficient should be reported, original, bootstrapped mean, or bootstrapped median}

\item{standard_error}{Should bootstrapped standard errors be reported?}

\item{p_value}{What type of p-value should be reported? Bootstrapped p_values, approximate p-values, or none?}

\item{confint}{What type of confidence intervals should be reported? Bootstrapped p_values, approximate p-values, or none?}

\item{conf_level}{Confidence level for confidence intervals}

\item{approx_t_value}{Should approximate t_values be reported}

\item{symmetric_bootstrap_p}{Should bootstrap p-values be computed as symmetric (leaving alpha/2 percent in each tail)? FALSE gives non-symmetric, but narrower, intervals. TRUE corresponds most closely to conventional p-values.}

\item{include_ylev}{Logical. Should y.lev be included in the tidy output? Makes for nicer tables when using modelsummary}

\item{...}{Other arguments to be passed to tidy}
}
\value{
A tidy function for a bootstrapped nb model
}
\description{
Tidy function for nbboot
}
\examples{
data(genevzinb2)
model <- evzinb(y~x1+x2+x3,data=genevzinb2, n_bootstraps = 10, multicore = TRUE, ncores = 2)
zinb_comp <- compare_models(model)
tidy(zinb_comp$nb) 

}
