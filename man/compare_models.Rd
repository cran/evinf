% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/zinb_comparison.R
\name{compare_models}
\alias{compare_models}
\title{Function to compare evzinb or evinb models with zinb and nb models}
\usage{
compare_models(
  object,
  nb_comparison = TRUE,
  zinb_comparison = TRUE,
  winsorize = FALSE,
  razorize = FALSE,
  cutoff_value = 10,
  init_theta = NULL,
  multicore = FALSE,
  ncores = NULL
)
}
\arguments{
\item{object}{A fitted evzinb or evinb model object}

\item{nb_comparison}{Should comparison be made with a negative binomial model?}

\item{zinb_comparison}{Should comparions be made with the zinb model?}

\item{winsorize}{Should winsorizing be done in the comparisons?}

\item{razorize}{Should razorizing (trimming) be done in the comparisons?}

\item{cutoff_value}{Integer: Which observation should be used as a basis for winsorizing/razorising. E.g. 10 means that everything larger than the 10th observation will be winsorized/razorised}

\item{init_theta}{Optional initial value for theta in the NB specification}

\item{multicore}{Logical: should multiple cores be used}

\item{ncores}{Number of cores if multicore is used}
}
\value{
A list with the original model as the first object and compared models as the following objects
}
\description{
Function to compare evzinb or evinb models with zinb and nb models
}
\examples{
data(genevzinb2)
model <- evzinb(y~x1+x2+x3,data=genevzinb2, n_bootstraps = 10, multicore = TRUE, ncores = 2)
compare_models(model)

}
