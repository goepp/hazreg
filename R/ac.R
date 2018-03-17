par2haz_ac <- function(par, J, K) {
  if (length(par) != J + K - 1) {
    stop("Error: length of param not equal to J + K - 1")
  }
  mu <- par[1]
  ext_alpha <- c(0, par[2:J])
  ext_beta <- c(0, par[(J + 1):(J + K - 1)])
  exp(outer(ext_beta, mu + ext_alpha, FUN = "+"))
}
loglik_ac <- function(par, O, R) {
  K <- nrow(O); J <- ncol(O)
  if (length(par) != J + K - 1) {
    stop("Error: length of param not equal to J + K - 1")
  }
  mu <- par[1]
  alpha <- c(0, par[2:J])
  beta <- c(0, par[(J + 1):(J + K - 1)])
  eta <- outer(beta, mu + alpha, FUN = "+")
  sum(exp(eta) * R - eta * O)
}
mle_ac <- function(O, R) {
  K <- nrow(O)
  J <- ncol(O)
  par_init <- rep(0, J + K - 1)
  par_res <- optim(par_init, fn = loglik_ac, O = O, R = R,
                   method = "BFGS", control = list(maxit = 5000))$par
  haz_res <- par2haz_ac(par_res, J, K)
  list("par" = par_res, "haz" = haz_res)
}
par2grid_ac <- function(par, cuts_age, cuts_cohort) {
  J <- length(cuts_age) + 1
  K <- length(cuts_cohort) + 1
  par2haz_ac(par, J, K) %>% haz2grid(cuts_age, cuts_cohort)
}
par2haz_ac <- function(par, J, K) {
  if (length(par) != J + K - 1) stop("Error: length of param not equal to J + K - 1")
  mu <- par[1]
  ext_alpha <- c(0, par[2:J])
  ext_beta <- c(0, par[(J + 1):(J + K - 1)])
  exp(outer(ext_beta, mu + ext_alpha, FUN = "+"))
}
#' Negative lok-likelihood in the age-cohort model
#'
#' @family ac_likelihood
#' @param par parameters of the age-cohort model in the form \code{c(mu, alpha and beta)}
#' @param O Observed events as returned by \code{\link{exhaustive_stat_2d}}
#' @param R Time at risk as returned by \code{\link{exhaustive_stat_2d}}
#' @return The negative log-likelihood as described in TODO NAME OF REFERENCE SHEET
#' @examples
#' par <- rnorm(ncol(O) + nrow(O))
#' loglik_ac(par, O, R)
loglik_ac <- function(par, O, R) {
  K <- nrow(O)
  J <- ncol(O)
  if (length(par) != J + K - 1) {
    stop("Error: length of param not equal to J + K - 1")
  }
  mu <- par[1]
  alpha <- c(0, par[2:J])
  beta <- c(0, par[(J + 1):(J + K - 1)])
  eta <- outer(beta, mu + alpha, FUN = "+")
  sum(exp(eta) * R - eta * O)
}
#' First order derivate of the negative lok-likelihood in the age-cohort model
#'
#' @family ac_likelihood
#' @param par parameters of the age-cohort model in the form \code{c(mu, alpha and beta)}
#' @param O Observed events as returned by \code{\link{exhaustive_stat_2d}}
#' @param R Time at risk as returned by \code{\link{exhaustive_stat_2d}}
#' @return The vector of derivatives of the negative log-likelihood as
#' described in TODO NAME OF REFERENCE SHEET
#' @examples
#' par <- rnorm(ncol(O) + nrow(O))
#' score_ac(par, O, R)
score_ac <- function(par, O, R) {
  K <- nrow(O)
  J <- ncol(O)
  if (length(par) != J + K - 1) {
    stop("error: length of param not equal to J + K - 1")
  }
  mu <- par[1]
  alpha <- par[2:J]
  ext_alpha <- c(0, alpha)
  beta <- par[(J + 1):(J + K - 1)]
  ext_beta <- c(0, beta)
  deriv_mu <- sum(exp(outer(ext_beta, mu + ext_alpha, FUN = "+")) * R - O)
  deriv_alpha <- sapply(2:J, function(ind_j) sum( exp(outer(ext_beta, mu + ext_alpha[ind_j], FUN = "+")) *
                                                    R[, ind_j] - O[, ind_j]))
  deriv_beta <- sapply(2:K, function(ind_k) sum( exp(outer(ext_beta[ind_k], mu + ext_alpha, FUN = "+")) *
                                                   R[ind_k, ] - O[ind_k, ]))
  c(deriv_mu, deriv_alpha, deriv_beta)
}
#' Second order derivate of the negative lok-likelihood in the age-cohort model
#'
#' @family ac_likelihood
#' @param par parameters of the age-cohort model in the form \code{c(mu, alpha and beta)}
#' @param O Observed events as returned by \code{\link{exhaustive_stat_2d}}
#' @param R Time at risk as returned by \code{\link{exhaustive_stat_2d}}
#' @return The matrix of second order derivatives of the negative log-likelihood as
#' described in TODO NAME OF REFERENCE SHEET
#' @examples
#' par <- rnorm(ncol(O) + nrow(O))
#' hessian_ac(par, O, R)
hessian_ac <- function(par, O, R) {
  K <- nrow(O)
  J <- ncol(O)
  if (length(par) != J + K - 1) {
    stop("Error: length of param not equal to J + K - 1")
  }
  mu <- par[1]
  alpha <- par[2:J]
  ext_alpha <- c(0, alpha)
  beta <- par[(J + 1):(J + K - 1)]
  ext_beta <- c(0, beta)
  deriv_diag_mu <- sum(exp(outer(ext_beta, mu + ext_alpha, FUN = "+")) * R)
  deriv_diag_alpha <- diag(sapply(2:J, function(ind_j) sum( exp(outer(ext_beta, mu + ext_alpha[ind_j], FUN = "+")) *
                                                              R[, ind_j])), J - 1, J - 1)
  deriv_diag_beta <- diag(sapply(2:K, function(ind_k) sum( exp(outer(ext_beta[ind_k], mu + ext_alpha, FUN = "+")) *
                                                             R[ind_k, ])), K - 1, K - 1)
  deriv_alpha_mu <- matrix(sapply(2:J, function(ind_j) sum( exp(outer(ext_beta, mu + ext_alpha[ind_j], FUN = "+")) *
                                                              R[, ind_j])), J - 1, 1)
  deriv_beta_mu <- matrix(sapply(2:K, function(ind_k) sum( exp(outer(ext_beta[ind_k], mu + ext_alpha, FUN = "+")) *
                                                             R[ind_k,])), K - 1, 1)
  deriv_beta_alpha <- (exp(outer(ext_beta, mu + ext_alpha, FUN = "+")) * R)[-1, -1]

  cbind(rbind(deriv_diag_mu, deriv_alpha_mu, deriv_beta_mu),
        rbind(t(deriv_alpha_mu), deriv_diag_alpha, deriv_beta_alpha),
        rbind(t(deriv_beta_mu), t(deriv_beta_alpha), deriv_diag_beta)) %>%
    "dimnames<-"(list(c("mu", rep("alpha", J - 1), rep("beta", K - 1)),
                      c("mu", rep("alpha", J - 1), rep("beta", K - 1))))
}
#' Solver of the Age-cohort model
#'
#' @family ac_likelihood
#' @param O Observed events as returned by \code{\link{exhaustive_stat_2d}}
#' @param R Time at risk as returned by \code{\link{exhaustive_stat_2d}}
#' @return The maximum likelihood estimate of the Age-Cohort Model. Minimization is made
#' using the Newton-Raphson algorithm (\url{https://en.wikipedia.org/wiki/Newton%27s_method})
#' @seealso \code{\link{exhaustive_stat_2d}}, \code{\link{exhaustive_stat_1d}}
#' @examples
#' par <- rnorm(ncol(O) + nrow(O))
#' hessian_ac(par, O, R)
solver_ac <- function(O, R, maxiter = 1000, verbose = FALSE) {
  K <- nrow(O)
  J <- ncol(O)
  old_par <- rep(0, J + K - 1)
  for (iter in 1:maxiter) {
    if (verbose) old_par %>% par2grid_ac(cuts_age, cuts_cohort) %>%
      make_plot() %>% print()
    par <- old_par - Solve(hessian_ac(old_par, O, R),
                           score_ac(old_par, O, R))
    if (sum(is.na(abs(par - old_par) / abs(old_par)))) break
    if (max(abs(par - old_par) / abs(old_par)) <= sqrt(.Machine$double.eps)) break
    old_par <- par
  }
  if (iter == maxiter) {
    warning("Warning: Newton-Raphson procedure did not converge")
  }
  list("par" = par, "convergence" = iter != maxiter, "niter" = iter)
}
