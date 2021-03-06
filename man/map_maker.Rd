% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utility.R
\name{map_maker}
\alias{map_maker}
\title{Create two dimensional hazard rate from Age-Cohort-Interaction model}
\usage{
map_maker(mu, age_coef, cohort_coef, age_seq, cohort_seq,
  islet_ls = NULL)
}
\arguments{
\item{mu}{Constant part of the log-hazard}

\item{age_coef}{Age-effect parameter sequence}

\item{cohort_coef}{Cohort-effect parameter sequence}

\item{age_seq}{Age discretizatioin vector}

\item{cohort_seq}{Cohort discretizatioin vector}

\item{islet_ls}{List of parameters of the interaction terms.}
}
\value{
The interactions terms are densities of two dimensional normal distributions.
\code{islet_ls} gives the mean, variance-covariance matrix and multiplicative coefficients of each
interaction term.
}
\description{
Create two dimensional hazard rate from Age-Cohort-Interaction model
}
\examples{
\dontrun{
library(lattice)
mu <- log(1e-2)
age_seq <- seq(0, 100, length.out = 100)
cohort_seq <- seq(1900, 1900 + age_seq[101], length.out = age_length + 1)
age_coef <- seq(0, 2, length.out = age_length - 1)
cohort_coef <- seq(0, 0.8, length.out = age_length - 1)
islet_ls = list(list(mean = c(1945, 45), sigma = 100 * diag(2), coef = 40))
hazard <- map_maker(mu, age_coef, cohort_coef, age_seq, cohort_seq, islet_ls = islet_ls)
persp(hazard$sheet)
}
}
