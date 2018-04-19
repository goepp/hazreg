#' @export
par2haz_interaction <- function(par, J, K) {
  if (length(par) != J * K) stop("error: length of param not equal to J * K")
  exp(matrix(par, K, J))
}
#' @export
par2grid_interaction <- function(par, cuts_age, cuts_cohort) {
  J <- length(cuts_age) + 1
  K <- length(cuts_cohort) + 1
  par_df <- par2haz_interaction(unname(par), J, K) %>%
    melt(varnames = c("cohort", "age")) %>%
    mutate(age = levels(cut(-1, breaks = c(0, cuts_age, Inf), right = FALSE, dig.lab = 3))[matrix(age)]) %>%
    mutate(age = factor(age, levels = unique(age))) %>%
    mutate(cohort = levels(cut(-1, breaks = c(0, cuts_cohort, Inf), right = FALSE, dig.lab = 4))[matrix(cohort)]) %>%
    mutate(cohort = factor(cohort, levels = unique(cohort)))
  par_df
}
#' @export
par_sel2par_interaction <- function(par_sel, sel) {
  J <- ncol(sel)
  K <- nrow(sel)
  interaction <- mapvalues(sel,
                           from = levels(as.factor(sel)),
                           to = par_sel) %>%
    as.vector()
  stopifnot(length(interaction) == J * K)
  interaction
}
#' @export
par2haz_sel_interaction <- function(par, sel, J, K, haz.log = FALSE) {
  L <- nlevels(as.factor(sel))
  if (length(par) != L) {
    stop("Error: selection nlevels not consistent with J * K")
  }
  eta <- sel[] * NA
  for (level_ind in 1:L) {
    eta[sel == level_ind] <- par[level_ind]
  }
  if (haz.log) {
    return(eta)
  } else {
    return(exp(eta))
  }
}
#' @export
loglik_sel_interaction <- function(par, O, R) {
  if (any(O < 0 & R < 0)) stop('Error: O and R must have non-negative values')
  sum(exp(par) * R - "[<-"(par * O, which(is.nan(par * O), arr.ind = TRUE), 0))
}
#' Negative lok-likelihood in the interaction model
#'
#' @family interaction_likelihood
#' @param par parameters of the interaction model in the form \code{c(mu, alpha and beta)}
#' @param O Observed events as returned by \code{\link{exhaustive_stat_2d}}
#' @param R Time at risk as returned by \code{\link{exhaustive_stat_2d}}
#' @return The negative log-likelihood as described in TODO NAME OF REFERENCE SHEET
#' @examples
#' \dontrun{
#' J <- 10
#' K <- 15
#' set.seed(0)
#' O <- matrix(rpois(K * J, 2), K, J)
#' R <- matrix(rpois(K * J, 10), K, J)
#' par <- rnorm(J * K)
#' loglik_interaction(par, O, R)
#' }
#' @export
loglik_interaction <- function(par, O, R, pen, weights_age = NULL,
                               weights_cohort = NULL) {
  if (any(dim(O) != dim(R)) | (length(par) != length(R))) stop('Error: dimensions of O, R, and pen must agree')
  if (O < 0 || R < 0 || pen < 0) stop('Error: O, R, and pen must have non-negative values')
  K <- nrow(O)
  J <- ncol(O)
  if (is.null(weights_age)) weights_age <- matrix(1, K , J - 1)
  if (is.null(weights_cohort)) weights_cohort <- matrix(1, K - 1, J)
  eta <- matrix(par, K, J)
  pen_term <- pen / 2 * (
    sum(weights_age * t(apply(eta, 1, diff)) ^ 2) +
      sum(weights_cohort * apply(eta, 2, diff) ^ 2)
  )
  sum(exp(eta) * R - eta * O) + pen_term
}
#' First order derivate of the negative lok-likelihood in the interaction model
#'
#' @family likelihood_likelihood
#' @param par parameters of the interaction model in the form \code{as.vector(eta)}
#' @param O Observed events as returned by \code{\link{exhaustive_stat_2d}}
#' @param R Time at risk as returned by \code{\link{exhaustive_stat_2d}}
#' @return The vector of derivatives of the negative log-likelihood as
#' described in TODO NAME OF REFERENCE SHEET
#' @examples
#' \dontrun{
#' J <- 10
#' K <- 15
#' set.seed(0)
#' O <- matrix(rpois(K * J, 2), K, J)
#' R <- matrix(rpois(K * J, 10), K, J)
#' par <- rnorm(J * K)
#' score_interaction(par, O, R)
#' }
#' @export
score_interaction <- function(par, O, R, pen, weights_age = NULL,
                              weights_cohort = NULL) {
  if (any(dim(O) != dim(R)) | (length(par) != length(R))) stop('Error: dimensions of O, R, and pen must agree')
  if (O < 0 || R < 0 || pen < 0) stop('Error: O, R, and pen must have non-negative values')
  K <- nrow(O)
  J <- ncol(O)
  if (is.null(weights_age)) weights_age <- matrix(1, K, J - 1)
  if (is.null(weights_cohort)) weights_cohort <- matrix(1, K - 1, J)
  eta <- matrix(par, K, J)
  unpenalized_score <- as.vector(exp(eta) * R - O)
  pen_term <- pen * as.vector(
    cbind(0, weights_age * t(apply(eta, 1, diff))) -
      cbind(weights_age * t(apply(eta, 1, diff)), 0)) +
    pen * as.vector(rbind(0, weights_cohort * apply(eta, 2, diff)) -
                      rbind(weights_cohort * apply(eta, 2, diff), 0))
  unpenalized_score + pen_term
}
#' Second order derivate of the negative lok-likelihood in the interaction model
#'
#' @family ac_likelihood
#' @param par parameters of the interaction model in the form \code{c(mu, alpha and beta)}
#' @param O Observed events as returned by \code{\link{exhaustive_stat_2d}}
#' @param R Time at risk as returned by \code{\link{exhaustive_stat_2d}}
#' @return The matrix of second order derivatives of the negative log-likelihood as
#' described in TODO NAME OF REFERENCE SHEET
#' @examples
#' \dontrun{
#' J <- 10
#' K <- 15
#' set.seed(0)
#' O <- matrix(rpois(K * J, 2), K, J)
#' R <- matrix(rpois(K * J, 10), K, J)
#' par <- rnorm(J * K)
#' hessian_interaction(par, O, R)
#' }
#' @export
hessian_interaction <- function(par, O, R, pen, weights_age = NULL,
                                weights_cohort = NULL, use_band = FALSE) {
  if (any(dim(O) != dim(R)) | (length(par) != length(R))) stop('Error: dimensions of O, R, and par must agree')
  if (O < 0 || R < 0 || pen < 0) stop('Error: O, R, and pen must have non-negative values')
  K <- nrow(O)
  J <- ncol(O)
  if (is.null(weights_age))    weights_age <- matrix(1, K, J - 1)
  if (is.null(weights_cohort)) weights_cohort <- matrix(1, K - 1, J)
  if (nrow(weights_age) != K | ncol(weights_age) != J - 1 | nrow(weights_cohort) != K - 1 |
      ncol(weights_cohort) != J) {
    stop('Error: dimensions of weights and exhaustive statistics must agree')
  }
  eta <- matrix(par, K, J)
  if (!use_band) {
    unpenalized_hessian <- diag(as.vector(exp(eta) * R))
    pen_term_diag <- diag(as.vector(
      cbind(0, weights_age) + cbind(weights_age, 0) +
        rbind(0, weights_cohort) + rbind(weights_cohort, 0)
    ))
    index_weights_age <- cbind(
      1:((J - 1) * K) + K,
      1:((J - 1) * K)
    )
    index_weights_cohort <- cbind(
      (1:(J * K))[-seq(1, K * J, by = K)],
      (1:(J * K))[-K * seq(1, J)]
    )
    pen_term_subdiag <- "[<-"("[<-"(matrix(0, K * J, K * J),
                                    index_weights_age, -as.vector(weights_age)),
                              index_weights_cohort, -as.vector(weights_cohort))
    matrix_hessian <- unpenalized_hessian + pen * (pen_term_diag + pen_term_subdiag +
                                                     t(pen_term_subdiag))
    return(matrix_hessian)
  } else {
    unpenalized_diag <- as.vector(exp(eta) * R)
    pen_diag <- as.vector(cbind(0, weights_age) + cbind(weights_age, 0) +
                            rbind(0, weights_cohort) + rbind(weights_cohort, 0))
    temp1 <- rep(0, K * J)
    temp1[(1:(J * K))[-K * seq(1, J)]] <- pen * -as.vector(weights_cohort)
    upper_diag_1 <- temp1
    upper_diag_Kplus1 <- c(pen * -as.vector(weights_age), rep(0, K))
    band_hessian <- cbind(unpenalized_diag + pen * pen_diag,
                          upper_diag_1,
                          matrix(rep(0, K * J * (K - 2)), K * J, K - 2),
                          upper_diag_Kplus1)
    return(band_hessian)
  }
}
ridge_solver_interaction_old <- function(O, R, pen, maxiter = 1000, verbose = FALSE,
                                         use_band = TRUE) {
  old_par <- rep(0, ncol(O) * nrow(O))
  for (iter in 1:maxiter) {
    if (verbose) {
      old_par %>% par2grid_interaction(cuts_age, cuts_cohort) %>%
        make_plot() %>% print()
    }
    if (use_band) {
      par <- old_par - bandsolve(hessian_interaction(old_par, O, R, pen, use_band = TRUE),
                                 score_interaction(old_par, O, R, pen))
    } else {
      par <- old_par - Solve(hessian_interaction(old_par, O, R, pen, use_band = FALSE),
                             score_interaction(old_par, O, R, pen))
    }
    if (sum(is.na(abs(par - old_par) / abs(old_par)))) break
    if (max(abs(par - old_par) / abs(old_par)) <= sqrt(1e-8)) break
    old_par <- par
  }
  if (iter == maxiter) {
    warning("Warning: Newton-Raphson procedure did not converge")
  }
  list("par" = par, "convergence" = iter != maxiter, "niter" = iter)
}
#' Computes the ridge regularized estimates of the interaction model
#'
#' @param O Observed events as returned by \code{\link{exhaustive_stat_2d}}
#' @param R Time at risk as returned by \code{\link{exhaustive_stat_2d}}
#' @param pen The penalty constant.
#' @param maxiter The maximal number of iteractions of the Newton-Raphson procedure.
#' @param verbose Whether to display the progress of the Newton-Raphson procedure.
#' @param use_band Whether to use faster inversion of the hessian function. \code{TRUE}
#' @param weights_age weights of the age differences in the log-hazard. See ADD REFERENCE SHEET
#' @param weights_cohort weights of the cohort differences in the log-hazard. See ADD REFERENCE SHEET
#' @param old_par initial point of the iteration. The default choice is \code{0}
#' @param tol Absolute error tolerance used as stopping criteria of the Newton-Raphson minimization.
#' is recommended
#' @return The vector estimate of the ridge regularized interaction model.
#' @examples
#' \dontrun{
#' J <- 10
#' K <- 15
#' set.seed(0)
#' O <- matrix(rpois(K * J, 2), K, J)
#' R <- matrix(rpois(K * J, 10), K, J)
#' pen <-  50
#' ridge <- ridge_solver_interaction(O, R, pen)
#' ridge$par
#' ridge$iter
#' }
#' @export
#' @family ridge
ridge_solver_interaction <- function(O, R, pen, weights_age = NULL,
                                     weights_cohort = NULL, use_band = TRUE,
                                     old_par = rep(0, ncol(O) * nrow(O)), maxiter = 1000,
                                     verbose = FALSE, tol = 1e-8) {
  for (iter in 1:maxiter) {
    if (verbose) old_par %>% raster('numeric')
    if (use_band) {
      par <- old_par - bandsolve(
        hessian_interaction(old_par, O, R, pen, weights_age, weights_cohort, use_band = TRUE),
        score_interaction(old_par, O, R, pen, weights_age, weights_cohort))
    } else {
      par <- old_par - Solve(
        hessian_interaction(old_par, O, R, pen, weights_age, weights_cohort, use_band = FALSE),
        score_interaction(old_par, O, R, pen, weights_age, weights_cohort))
    }
    if (sum(is.na(abs(par - old_par) / abs(old_par)))) break
    if (max(abs(par - old_par) / abs(old_par)) <= 1e-8) break
    old_par <- par
  }
  if (iter == maxiter) {
    warning("Warning: Newton-Raphson procedure did not converge")
  }
  list("par" = par, "convergence" = iter != maxiter, "niter" = iter)
}
#' @rdname ridge_solver_interaction
#' @export
#' @family adaptive_ridge
aridge_solver_interaction <- function(O, R, pen, sample_size,
                                      use_band = TRUE,
                                      maxiter = 1000,
                                      verbose = FALSE,
                                      progress = TRUE) {
  if (progress) {
    pb <- progress_bar$new(
      format = "  adaptive ridge [:bar] :percent in :elapsed",
      total = length(pen), clear = FALSE, width = 100)
    pb$tick(0)
  }
  sel <- par_sel <- haz_sel <- vector('list', length(pen))
  bic <- aic <- ebic <- NA * pen
  epsilon_age <- 1e-6
  epsilon_cohort <- 1e-6
  K <- nrow(O)
  J <- ncol(O)
  weights_age <- matrix(1, K, J - 1)
  weights_cohort <- matrix(1, K - 1, J)
  old_par <- rep(0, K * J)
  valve_age <- weights_age %>%
    "colnames<-"(diag(sapply(colnames(O)[-1], paste0, "-",
                             colnames(O)[-length(cuts_age) - 1]))) %>%
    "rownames<-"(rownames(O))
  valve_cohort <- weights_cohort %>%
    "rownames<-"(diag(sapply(rownames(O)[-1], paste0, "-",
                             rownames(O)[-length(cuts_age) - 1]))) %>%
    "colnames<-"(colnames(O))
  ind_pen <- 1
  for (iter in 1:maxiter) {
    old_valve_age <- valve_age
    old_valve_cohort <- valve_cohort
    par <- ridge_solver_interaction(O, R, pen[ind_pen],
                                    weights_age, weights_cohort,
                                    use_band,
                                    old_par)$par
    eta <- matrix(par, K, J)
    weights_age[, ]    <- 1 / (t(apply(eta, 1, diff)) ^ 2 + epsilon_age ^ 2)
    weights_cohort[, ] <- 1 / (apply(eta, 2, diff) ^ 2 + epsilon_cohort ^ 2)
    valve_age[, ] <- (weights_age * t(apply(eta, 1, diff)) ^ 2)
    valve_cohort[, ] <- (weights_cohort * apply(eta, 2, diff) ^ 2)
    converge <- max(abs(old_valve_age - valve_age),
                      abs(old_valve_cohort - valve_cohort)) <= 1e-6
    if (converge) {
      sel[[ind_pen]] <- valve2sel(valve_age, valve_cohort)
      if (verbose) raster(sel[[ind_pen]], 'factor')
      exhaust_sel <- exhaustive_stat_sel(list("O" = O, "R" = R), sel[[ind_pen]])
      par_sel[[ind_pen]] <- log(exhaust_sel$O / exhaust_sel$R) %>%
        '[<-'(which(log(exhaust_sel$O / exhaust_sel$R) == -Inf), -37) %>%
        '[<-'(which(is.nan(log(exhaust_sel$O / exhaust_sel$R))), 0)
      haz_sel[[ind_pen]] <- par2haz_sel_interaction(
        par_sel[[ind_pen]], sel[[ind_pen]], J, K, haz.log = FALSE)
      bic[ind_pen] <- log(sample_size) * length(par_sel[[ind_pen]]) +
        2 * loglik_sel_interaction(par_sel[[ind_pen]], exhaust_sel$O, exhaust_sel$R)
      ebic[ind_pen] <- bic[ind_pen] +
        2 * log(choose(K * J, length(par_sel[[ind_pen]])))
      aic[ind_pen] <- 2 * length(par_sel[[ind_pen]]) +
        2 * loglik_sel_interaction(par_sel[[ind_pen]], exhaust_sel$O, exhaust_sel$R)
      if (progress) pb$tick()
      ind_pen <-  ind_pen + 1
    }
    old_par <- par
    if (ind_pen == length(pen) + 1) break
  }
  if (iter == maxiter) {
    warning("Warning: aridge did not converge")
  }
  return(list("sel" = sel, "par" = par_sel, "haz" = haz_sel,
              "bic" = bic, "aic" = aic, 'ebic' = ebic))
}
aridge_solver_interaction_q <- function(O, R, pen, sample_size, q = 2,
                                        use_band = TRUE,
                                        maxiter = 2000) {
  sel <- par_sel <- haz_sel <- vector('list', length(pen))
  bic <- aic <- ebic <- NA * pen
  epsilon_age <- 1e-6
  epsilon_cohort <- 1e-6
  K <- nrow(O)
  J <- ncol(O)
  weights_age <- matrix(1, K, J - 1)
  weights_cohort <- matrix(1, K - 1, J)
  old_par <- rep(0, K * J)
  ind_pen <- 1
  for (iter in 1:maxiter) {
    par <- ridge_solver_interaction(O, R, pen[ind_pen],
                                    weights_age, weights_cohort,
                                    old_par = old_par,
                                    use_band = use_band)$par
    eta <- matrix(par, K, J)
    weights_age[, ]    <- 1 / (t(apply(eta, 1, diff)) ^ 2 + epsilon_age ^ 2) ^ (1 - q / 2)
    weights_cohort[, ] <- 1 / (apply(eta, 2, diff) ^ 2 + epsilon_cohort ^ 2) ^ (1 - q / 2)
    if (sum(is.na(abs(par - old_par) / abs(old_par))) ||
        (max(abs(par - old_par) / abs(old_par)) <= 1e-8)) {
      valve_age <- (weights_age * t(apply(eta, 1, diff)) ^ 2) %>% round(digits = 15) %>%
        "colnames<-"(diag(sapply(colnames(O)[-1], paste0, "-",
                                 colnames(O)[-length(cuts_age) - 1]))) %>%
        "rownames<-"(rownames(O))
      valve_cohort <- (weights_cohort * apply(eta, 2, diff) ^ 2) %>% round(digits = 15) %>%
        "rownames<-"(diag(sapply(rownames(O)[-1], paste0, "-",
                                 rownames(O)[-length(cuts_age) - 1]))) %>%
        "colnames<-"(colnames(O))
      sel[[ind_pen]] <- valve2sel(valve_age, valve_cohort)
      exhaust_sel <- exhaustive_stat_sel(list("O" = O, "R" = R), sel[[ind_pen]])
      par_sel[[ind_pen]] <- log(exhaust_sel$O / exhaust_sel$R) %>%
        '[<-'(which(log(exhaust_sel$O / exhaust_sel$R) == -Inf), -37) %>%
        '[<-'(which(is.nan(log(exhaust_sel$O / exhaust_sel$R))), 0)
      haz_sel[[ind_pen]] <- par2haz_sel_interaction(par_sel[[ind_pen]], sel[[ind_pen]],
                                                       J, K, haz.log = FALSE)
      bic[ind_pen] <- log(sample_size) * length(par_sel[[ind_pen]]) +
        2 * loglik_sel_interaction(par_sel[[ind_pen]], exhaust_sel$O, exhaust_sel$R)
      ebic[ind_pen] <- bic[ind_pen] +
        2 * log(choose(nrow(O) * ncol(O), length(par_sel[[ind_pen]])))
      aic[ind_pen] <- 2 * length(par_sel[[ind_pen]]) +
        2 * loglik_sel_interaction(par_sel[[ind_pen]], exhaust_sel$O, exhaust_sel$R)
      ind_pen <-  ind_pen + 1
    }
    old_par <- par
    if (ind_pen == length(pen) + 1) break
  }
  if (iter == maxiter) {
    warning("Warning: aridge did not converge")
  }
  return(list("sel" = sel, "par" = par_sel, "haz" = haz_sel,
              "bic" = bic, "aic" = aic, 'ebic' = ebic))
}
#' @export
cv_aridge_interaction <- function(data, cuts_age, cuts_cohort,
                                  pen = 10 ^ seq(-4, 4, length = 50), nfold = 10) {
  score_matrix <- matrix(NA, nfold, length(pen))
  for (ind in 1:nfold) {
    sample_test <- seq(floor(nrow(data)/nfold) * (ind - 1) + 1,
                       floor(nrow(data)/nfold) * ind)
    sample_train <- setdiff(1:nrow(data), sample_test)
    exhaust_train <- exhaustive_stat(dplyr::slice(data, sample_train), cuts_age, cuts_cohort)
    train <- aridge_solver_interaction(exhaust_train$O, exhaust_train$R, pen = pen,
                                       sample_size = length(sample_train))
    train_par_ls <- lapply(train$par, function(par) '[<-'(par, which(is.nan(par)), 0))
    train_sel_ls <- train$sel
    exhaust_test_ls <-  lapply(train_sel_ls, function(sel) exhaustive_stat_sel(
      exhaustive_stat(dplyr::slice(data, sample_test),  cuts_age, cuts_cohort), sel))
    score_matrix[ind, ] <- mapply(
      FUN = function(par, exhaust_test) loglik_sel_interaction(par, exhaust_test$O, exhaust_test$R),
      par = train_par_ls, exhaust_test = exhaust_test_ls)
    if (any(is.null(score_matrix[ind, ]))) {
      stop('Error in call to aridge_solver_interaction')
    }
  }
  score <- colSums(score_matrix)
  exhaust <- exhaustive_stat_2d(data, cuts_age, cuts_cohort)
  fit <- aridge_solver_interaction(exhaust$O, exhaust$R, pen[which.min(score)], nrow(data))
  list('haz' = fit$haz[[1]],
       'par' = fit$par[[1]],
       'sel' = fit$sel[[1]],
       'score' = colSums(score_matrix))
}
#' @export
cv_ridge_interaction <- function(data, cuts_age, cuts_cohort,
                                 pen = 10 ^ seq(-4, 4, length = 50), nfold = 10) {
  score_matrix <- matrix(NA, nfold, length(pen))
  for (ind in 1:nfold) {
    sample_test <- seq(floor(nrow(data) / nfold) * (ind - 1) + 1,
                       floor(nrow(data) / nfold) * ind)
    sample_train <- setdiff(1:nrow(data), sample_test)
    exhaust_train <- exhaustive_stat_2d(dplyr::slice(data, sample_train), cuts_age, cuts_cohort)
    exhaust_test <- exhaustive_stat_2d(dplyr::slice(data, sample_test), cuts_age, cuts_cohort)
    train_par_ls <- lapply(pen, ridge_solver_interaction, O = exhaust_train$O, R = exhaust_train$R) %>%
      lapply(., function(element) element$par)
    score_matrix[ind, ] <- lapply(train_par_ls, loglik_interaction, exhaust_test$O, exhaust_test$R, pen = 0) %>%
      unlist()
    if (any(is.null(score_matrix[ind, ]))) {
      stop('Error in call to ridge_solver_interaction')
    }
  }
  score <- colSums(score_matrix)
  exhaust <- exhaustive_stat_2d(data, cuts_age, cuts_cohort)
  par <- ridge_solver_interaction(exhaust$O, exhaust$R, pen[which.min(score)])$par
  haz <- par2haz_interaction(par, length(cuts_age) + 1, length(cuts_cohort) + 1)
  list('haz' = haz,
       'par' = par,
       'score' = colSums(score_matrix))
}
#' @export
bootstrap <- function(data, cuts_age, cuts_cohort, times = 100,
                      pen = 10 ^ seq(-3, 4, length = 50),
                      crit = c("ebic")) {
  # fit <- vector("list", length = times)
  model_wrapper <- function(index, data, cuts_age, cuts_cohort) {
    sample <- data[sample(1:nrow(data), replace = TRUE), ]
    exhaust <- exhaustive_stat_2d(sample, cuts_age, cuts_cohort)
    fit <- aridge_solver_interaction(exhaust$O, exhaust$R, pen, nrow(sample), progress = TRUE)
    fit$haz[[which.min(fit[[crit]])]]
  }
  fit <- mclapply(1:times, model_wrapper, data, cuts_age, cuts_cohort)
  list('col' = fit,
       'median' = apply(simplify2array(fit), 1:2, median),
       'mean' = apply(simplify2array(fit), 1:2, mean),
       'sd' = apply(simplify2array(fit), 1:2, sd))
}


