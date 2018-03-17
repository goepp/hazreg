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
par2haz_interaction <- function(par, J, K) {
  if (length(par) != J * K) stop("error: length of param not equal to J * K")
  exp(matrix(par, K, J))
}
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
loglik_sel_interaction <- function(par, O, R) {
  sum(exp(par) * R - "[<-"(par * O, which(is.nan(par * O), arr.ind = TRUE), 0))
}
loglik_interaction <- function(par, O, R, pen, weights_age = NULL,
                               weights_cohort = NULL) {
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
score_interaction <- function(par, O, R, pen, weights_age = NULL,
                              weights_cohort = NULL) {
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
hessian_interaction <- function(par, O, R, pen, weights_age = NULL,
                                weights_cohort = NULL, band = FALSE) {
  K <- nrow(O)
  J <- ncol(O)
  if (any(dim(O) != dim(R))) {
    stop('Error: dimensions of O and R must agree')
  }
  if (length(par) != J * K) {
    stop('Error: dimensions of parameter and exhaustive statistics must agree')
  }
  if (is.null(weights_age)) weights_age <- matrix(1, K, J - 1)
  if (is.null(weights_cohort)) weights_cohort <- matrix(1, K - 1, J)
  if (nrow(weights_age) != K || ncol(weights_age) != J - 1 || nrow(weights_cohort) != K - 1 ||
      ncol(weights_cohort) != J) {
    stop('Error: dimensions of weights and exhaustive statistics must agree')
  }
  eta <- matrix(par, K, J)
  if (!band) {
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
ridge_solver_interaction <- function(O, R, pen, maxiter = 1000, verbose = FALSE,
                                     use_band = FALSE) {
  old_par <- rep(0, ncol(O) * nrow(O))
  for (iter in 1:maxiter) {
    if (verbose) {
      old_par %>% par2grid_interaction(cuts_age, cuts_cohort) %>%
        make_plot() %>% print()
    }
    if (use_band) {
      par <- old_par - bandsolve(hessian_interaction(old_par, O, R, pen, band = TRUE),
                                 score_interaction(old_par, O, R, pen))
    } else {
      par <- old_par - Solve(hessian_interaction(old_par, O, R, pen, band = FALSE),
                             score_interaction(old_par, O, R, pen))
    }
    if (sum(is.na(abs(par - old_par) / abs(old_par)))) break
    if (max(abs(par - old_par) / abs(old_par)) <= sqrt(.Machine$double.eps)) break
    old_par <- par
  }
  if (iter == maxiter) {
    warning("Warning: Newton-Raphson procedure did not converge")
  }
  list("par" = par, "convergence" = iter != maxiter, "niter" = iter)
}
wridge_solver_interaction <- function(O, R, pen, weights_age, weights_cohort,
                                      maxiter = 1000, verbose = FALSE, old_par = NULL) {
  if (is.null(old_par)) old_par <- rep(0, ncol(O) * nrow(O))
  for (iter in 1:maxiter) {
    if (verbose) old_par %>% par2grid_interaction(cuts_age, cuts_cohort) %>%
      make_plot() %>% print()
    par <- old_par - bandsolve(
      hessian_interaction(old_par, O, R, pen, weights_age, weights_cohort, band = TRUE),
      score_interaction(old_par, O, R, pen, weights_age, weights_cohort))
    if (sum(is.na(abs(par - old_par) / abs(old_par)))) break
    if (max(abs(par - old_par) / abs(old_par)) <= sqrt(.Machine$double.eps)) break
    old_par <- par
  }
  if (iter == maxiter) {
    warning("Warning: Newton-Raphson procedure did not converge")
  }
  list("par" = par, "convergence" = iter != maxiter, "niter" = iter)
}
aridge_solver_interaction <- function(O, R, pen_vect, sample_size,
                                      maxiter = 1000 * length(pen_vect)) {
  sel_ls <- par_sel_ls <- haz_sel_ls <- bic_ls <- aic_ls <- ebic_ls <- vector('list', length(pen_vect))
  epsilon_age <- 1e-6
  epsilon_cohort <- 1e-6
  K <- nrow(O)
  J <- ncol(O)
  weights_age <- matrix(1, K, J - 1)
  weights_cohort <- matrix(1, K - 1, J)
  old_par <- rep(0, J * K)
  ind_pen <- 1
  for (iter in 1:maxiter) {
    par <- wridge_solver_interaction(O = O, R = R, pen = pen_vect[ind_pen], weights_age = weights_age,
                                     weights_cohort = weights_cohort,
                                     maxiter = 1000, old_par = old_par)$par
    eta <- matrix(par, K, J)
    weights_age[, ]    <- 1 / (t(apply(eta, 1, diff)) ^ 2 + epsilon_age ^ 2)
    weights_cohort[, ] <- 1 / (apply(eta, 2, diff) ^ 2 + epsilon_cohort ^ 2)
    if (sum(is.na(abs(par - old_par) / abs(old_par))) ||
        (max(abs(par - old_par) / abs(old_par)) <= sqrt(.Machine$double.eps))) {
      valve_age <- (weights_age * t(apply(eta, 1, diff)) ^ 2) %>% round(digits = 15) %>%
        "colnames<-"(diag(sapply(colnames(O)[-1], paste0, "-",
                                 colnames(O)[-length(cuts_age) - 1]))) %>%
        "rownames<-"(rownames(O))
      valve_cohort <- (weights_cohort * apply(eta, 2, diff) ^ 2) %>% round(digits = 15) %>%
        "rownames<-"(diag(sapply(rownames(O)[-1], paste0, "-",
                                 rownames(O)[-length(cuts_age) - 1]))) %>%
        "colnames<-"(colnames(O))
      sel_ls[[ind_pen]] <- valve2sel(valve_age, valve_cohort)
      exhaust_sel <- exhaustive_stat_sel(list("O" = O, "R" = R), sel_ls[[ind_pen]])
      par_sel_ls[[ind_pen]] <- log(exhaust_sel$O / exhaust_sel$R) %>%
        '[<-'(which(log(exhaust_sel$O / exhaust_sel$R) == -Inf), -37) %>%
        '[<-'(which(is.nan(log(exhaust_sel$O / exhaust_sel$R))), 0)
      haz_sel_ls[[ind_pen]] <- par2haz_sel_interaction(par_sel_ls[[ind_pen]], sel_ls[[ind_pen]],
                                                       J, K, haz.log = FALSE)
      bic_ls[[ind_pen]] <- log(sample_size) * length(par_sel_ls[[ind_pen]]) +
        2 * loglik_sel_interaction(par_sel_ls[[ind_pen]], exhaust_sel$O, exhaust_sel$R)
      ebic_ls[[ind_pen]] <- bic_ls[[ind_pen]] +
        2 * log(choose(nrow(O) * ncol(O), length(par_sel_ls[[ind_pen]])))
      aic_ls[[ind_pen]] <- 2 * length(par_sel_ls[[ind_pen]]) +
        2 * loglik_sel_interaction(par_sel_ls[[ind_pen]], exhaust_sel$O, exhaust_sel$R)
      ind_pen <-  ind_pen + 1
    }
    old_par <- par
    if (ind_pen == length(pen_vect) + 1) break
  }
  if (iter == maxiter) {
    warning("Warning: aridge did not converge")
  }
  return(list("sel" = sel_ls, "par" = par_sel_ls,
              "haz" = haz_sel_ls, "bic" = bic_ls, "aic" = aic_ls, 'ebic' = ebic_ls))
}
