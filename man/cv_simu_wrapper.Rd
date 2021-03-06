% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simu.R
\name{cv_simu_wrapper}
\alias{cv_simu_wrapper}
\title{Run the interaction model with cross-validation on simulated data}
\usage{
cv_simu_wrapper(ind_wrapper, design, sample_df, pen_vect,
  scaled_ground_haz, nfold)
}
\arguments{
\item{nfold}{constant \code{K} in K-fold cross-validation}
}
\description{
Run the interaction model with cross-validation on simulated data
}
\seealso{
Other simu_wrapper: \code{\link{ac_simu_wrapper}},
  \code{\link{criterion_simu_wrapper}},
  \code{\link{design_simu}},
  \code{\link{ridge_simu_wrapper}}
}
\concept{simu_wrapper}
