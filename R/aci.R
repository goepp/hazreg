#' Transform parameter to hazard in the age-cohort-interaction model
#'
#' @param par parameter of the age-cohort-interaction model in the
#' form \code{c(mu, alpha and beta)}
#' @param J length of \code{cuts_age}
#' @param K length of \code{cuts_cohort}
#' @return a matrix of dimension \code{K * J}
#' @examples
#' \dontrun{
#' J <-  10
#' K <-  15
#' set.seed(0)
#' persp(par2haz_aci(rnorm(J * K), K, J))
#' }
#' @export
par2haz_aci <- function(par, J, K) {
  if (length(par) != J * K) {
    stop("Error: length of param not equal to J * K")
  }
  mu <- par[1]
  ext_alpha <- c(0, par[2:J])
  ext_beta <- c(0, par[(J + 1):(J + K - 1)])
  ext_delta <- "[<-"(matrix(0, K, J), -1, -1, matrix(par[-(1:(J + K - 1))], K - 1, J - 1))
  exp(outer(ext_beta, mu + ext_alpha, FUN = "+") + ext_delta)
}
#' @rdname par2haz_aci
#' @export
par2grid_aci <- function(par, cuts_age, cuts_cohort) {
  J <- length(cuts_age) + 1
  K <- length(cuts_cohort) + 1
  par_df <- par2haz_aci(unname(par), J, K) %>%
    melt(varnames = c("cohort", "age")) %>%
    mutate(age = levels(cut(
      -1,
      breaks = c(0, cuts_age, Inf),
      right = FALSE,
      dig.lab = 3))[matrix(age)]) %>%
    mutate(age = factor(age, levels = unique(age))) %>%
    mutate(cohort = levels(cut(
      -1,
      breaks = c(0, cuts_cohort, Inf),
      right = FALSE,
      dig.lab = 4))[matrix(cohort)]) %>%
    mutate(cohort = factor(cohort, levels = unique(cohort)))
  par_df
}
#' @rdname par2haz_aci
#' @export
par_sel2par_aci <- function(par_sel, sel) {
  J <- ncol(sel) + 1
  K <- nrow(sel) + 1
  interaction_sel <- par_sel[-(1:(J + K - 1))]
  interaction <- as.vector(mapvalues(sel,
                                     from = levels(as.factor(sel)),
                                     to = interaction_sel))
  stopifnot(length(interaction) == (J - 1) * (K - 1))
  c(par_sel[1:(J + K - 1)], interaction)
}
#' @export
haz_aci2i <- function(par, sel, J, K,
                      cohort_tab_center = NULL,
                      age_tab_center = NULL) {
  par_i <- par_sel2par_aci(par, sel) *
    (1:(K * J) >= (J + K))
  par2haz_aci(par_i, J, K) %>%
    '['(-1, -1) %>%
    'dimnames<-'(list(cohort_tab_center[-1],
                      age_tab_center[-1]))
}
#' @export
valve2sel_aci <- function(valve_age, valve_cohort, epsilon = 1e-8) {
  if (any(dim(valve_age) != dim(valve_cohort))) {
    stop("Error: dimensions of valve matrices must agree")
  }
  K <- nrow(valve_age) + 1
  J <- ncol(valve_age) + 1;
  adjacency <- diag((J - 1) * (K  - 1));
  node_names <- matrix(1:((J - 1) * (K - 1)), K - 1, J - 1)
  for (j in 1:(J - 1)) {
    for (k in 1:(K - 1)) {
      if (k > 1) {
        if (valve_cohort[k - 1, j] < epsilon) {
          adjacency[node_names[k, j], node_names[k - 1, j]] <-
            adjacency[node_names[k - 1, j], node_names[k, j]] <- 1
        }
      }
    }
  }
  for (j in 1:(J - 1)) {
    for (k in 1:(K - 1)) {
      if (j > 1) {
        if (valve_age[k, j - 1] < epsilon) {
          adjacency[node_names[k, j], node_names[k, j - 1]] <-
            adjacency[node_names[k, j - 1], node_names[k, j]] <- 1
        }
      }
    }
  }
  graph <- graph_from_adjacency_matrix(adjacency, mode = "undirected")
  sel_fct <- (matrix(clusters(graph)$membership, K - 1, J - 1) %>%
                'rownames<-'(rownames(valve_age)) %>%
                'colnames<-'(colnames(valve_cohort)))
  L <- nlevels(as.factor(sel_fct))
  sel_array <- lapply(1:L, function(ind) sel_fct == ind) %>%
    unlist() %>%
    array(., dim = c(K - 1, J - 1, L))
  list('fct' = sel_fct, 'array' = sel_array)
}
#' @export
par2haz_sel_aci_old <- function(par, sel, J, K, haz.log = FALSE) {
  L <- nlevels(as.factor(sel))
  if (length(par) != J + K - 1 + L) {
    stop("Error: selection nlevels not consistent with J * K")
  }
  mu <- par[1]
  ext_alpha <- c(0, par[2:J])
  ext_beta <- c(0, par[(J + 1):(J + K - 1)])
  z <- par[-(1:(J + K - 1))]
  delta <- matrix(NA, K - 1, J - 1)
  for (level_ind in 1:nlevels(as.factor(sel))) {
    delta[sel == level_ind] <- z[level_ind]
  }
  ext_delta <- "[<-"(matrix(0, K, J), -1, -1, delta)
  if (haz.log) outer(ext_beta, mu + ext_alpha, FUN = "+") + ext_delta
  else exp(outer(ext_beta, mu + ext_alpha, FUN = "+") + ext_delta)
}
#' @export
loglik_sel_aci_old <- function(par, O, R, sel) {
  K <- nrow(O)
  J <- ncol(O)
  eta <- par2haz_sel(par, sel, J, K, haz.log = TRUE)
  sum(exp(eta) * R - eta * O)
}
#' Negative lok-likelihood in the age-cohort-interaction model
#'
## @family aci_likelihood
#' @param par parameter of the age-cohort-interaction model in the
#' form \code{c(mu, alpha and beta)}
#' @param O Observed events as returned by \code{\link{exhaustive_stat_2d}}
#' @param R Time at risk as returned by \code{\link{exhaustive_stat_2d}}
#' @return The vector of derivatives of the negative log-likelihood as
#' described in TODO NAME OF REFERENCE SHEET
#' @examples
#' \dontrun{
#' par <- rnorm(ncol(O) * nrow(O))
#' loglik_aci(par, O, R) # length 1 vector
#' score_aci(par, O, R) # vector of length (ncol(O) + nrow(O))
#' hessian_aci(par, O, R) # square matrix of dimension (ncol(O) + nrow(O))
#' }
#' @export
loglik_aci  <- function(par, O, R, pen, weights_age = NULL,
                        weights_cohort = NULL) {
  if (any(dim(O) != dim(R)) & (length(pen) != length(R))) stop('Error: dimensions of O, R, and pen must agree')
  if (O < 0 || R < 0 || pen < 0) stop('Error: O, R, and pen must have non-negative values')
  K <- nrow(O)
  J <- ncol(O)
  if (is.null(weights_age)) weights_age <- matrix(1, K - 1, J - 1)
  if (is.null(weights_cohort)) weights_cohort <- matrix(1, K - 1, J - 1)
  mu <- par[1]
  alpha <- c(0, par[2:J])
  beta <- c(0, par[(J + 1):(J + K - 1)])
  delta <- matrix(par[-(1:(J + K - 1))], K - 1, J - 1)
  ext_delta <- "[<-"(matrix(0, K, J), -1, -1, delta)
  eta <- outer(beta, mu + alpha, FUN = "+") + ext_delta
  pen_term <- pen / 2 * (
    sum(weights_age * t(apply(ext_delta[-1, , drop = FALSE], 1, diff)) ^ 2) +
      sum(weights_cohort * apply(ext_delta[, -1, drop = FALSE], 2, diff) ^ 2)
  )
  sum(exp(eta) * R - eta * O) + pen_term
}
#' @rdname loglik_aci
#' @export
score_aci <- function(par, O, R, pen, weights_age = NULL,
                      weights_cohort = NULL) {
  if (any(dim(O) != dim(R)) & (length(pen) != length(R))) stop('Error: dimensions of O, R, and pen must agree')
  if (O < 0 || R < 0 || pen < 0) stop('Error: O, R, and pen must have non-negative values')
  K <- nrow(O)
  J <- ncol(O)
  if (is.null(weights_age)) weights_age <- matrix(1, K - 1, J - 1)
  if (is.null(weights_cohort)) weights_cohort <- matrix(1, K - 1, J - 1)
  mu <- par[1]
  alpha <- par[2:J]
  ext_alpha <- c(0, alpha)
  beta <- par[(J + 1):(J + K - 1)]
  ext_beta <- c(0, beta)
  delta <- matrix(par[-(1:(J + K - 1))], K - 1, J - 1)
  ext_delta <- "[<-"(matrix(0, K, J), -1, -1, delta)
  deriv_mu <- sum(exp(outer(ext_beta, mu + ext_alpha, FUN = "+") + ext_delta) * R - O)
  deriv_alpha <- sapply(2:J, function(ind_j) sum( exp(outer(ext_beta, mu + ext_alpha[ind_j], FUN = "+")
                                                      + ext_delta[, ind_j]) * R[, ind_j] - O[, ind_j]))
  deriv_beta <- sapply(2:K, function(ind_k) sum( exp(outer(ext_beta[ind_k], mu + ext_alpha, FUN = "+")
                                                     + ext_delta[ind_k, ]) * R[ind_k, ] - O[ind_k, ]))
  pen_term <- pen *  (
    weights_age * t(apply(ext_delta[-1, , drop = FALSE], 1, diff)) -
      (weights_age * t(apply(ext_delta[-1, , drop = FALSE], 1, diff))) %>% cbind(0) %>% "["(1:(K - 1), -1) +
      weights_cohort * apply(ext_delta[, -1, drop = FALSE], 2, diff) -
      (weights_cohort * apply(ext_delta[, -1, drop = FALSE], 2, diff)) %>% rbind(0) %>% "["(-1, 1:(J - 1))
  )
  deriv_delta_mat <- (exp(outer(ext_beta, mu + ext_alpha, FUN = "+") + ext_delta) * R - O)[-1, -1]
  deriv_delta <- as.vector(deriv_delta_mat + pen_term)
  c(deriv_mu, deriv_alpha, deriv_beta, deriv_delta) %>%
    setNames(c('mu', rep('alpha', J - 1), rep('beta', K - 1), rep('delta', (J - 1) * (K - 1))))
}
#' @rdname loglik_aci
#' @export
hessian_aci <- function(par, O, R, pen, weights_age = NULL,
                        weights_cohort = NULL, use_band = FALSE) {
  if (any(dim(O) != dim(R)) | (length(par) != length(R))) stop('Error: dimensions of O, R, and pen must agree')
  if (O < 0 || R < 0 || pen < 0) stop('Error: O, R, and pen must have non-negative values')
  K <- nrow(O)
  J <- ncol(O)
  if (K <= 2 || J <= 2) stop("Error: J and K must both be equal or greater than 3")
  if (is.null(weights_age)) {
    weights_age <- matrix(1, K - 1, J - 1)
  }
  if (is.null(weights_cohort)) {
    weights_cohort <- matrix(1, K - 1, J - 1)
  }
  mu <- par[1]
  alpha <- par[2:J]
  ext_alpha <- c(0, alpha)
  beta <- par[(J + 1):(J + K - 1)]
  ext_beta <- c(0, beta)
  delta <- matrix(par[-(1:(J + K - 1))], K - 1, J - 1)
  ext_delta <- "[<-"(matrix(0, K, J), -1, -1, delta)
  deriv_diag_mu <- sum(exp(outer(ext_beta, mu + ext_alpha, FUN = "+") + ext_delta) * R)
  deriv_diag_alpha <- diag(sapply(2:J,
                                  function(ind_j) sum(exp(outer(ext_beta, mu + ext_alpha[ind_j], FUN = "+") +
                                                            ext_delta[, ind_j, drop = FALSE]) * R[, ind_j])), J - 1, J - 1)
  deriv_diag_beta <- diag(sapply(2:K,
                                 function(ind_k) sum(exp(outer(ext_beta[ind_k], mu + ext_alpha, FUN = "+") +
                                                           ext_delta[ind_k, ]) * R[ind_k, ])), K - 1, K - 1)
  deriv_alpha_mu <- matrix(sapply(2:J,
                                  function(ind_j) sum(exp(outer(ext_beta, mu + ext_alpha[ind_j], FUN = "+") +
                                                            ext_delta[, ind_j]) * R[, ind_j])), J - 1, 1)
  deriv_beta_mu <- matrix(sapply(2:K,
                                 function(ind_k) sum(exp(outer(ext_beta[ind_k], mu + ext_alpha, FUN = "+") +
                                                           ext_delta[ind_k,]) * R[ind_k,])), K - 1, 1)
  deriv_delta_mu <- matrix(as.vector((exp(outer(ext_beta, mu + ext_alpha, FUN = "+") + ext_delta) * R)
                                     [-1, -1, drop = FALSE]), (J - 1) * (K - 1), 1)
  deriv_beta_alpha <- (exp(outer(ext_beta, mu + ext_alpha, FUN = "+") + ext_delta) * R)[-1, -1, drop = FALSE]
  deriv_delta_alpha <- "[<-"(matrix(0, (J - 1) * (K - 1), J - 1),
                             cbind(1:((K - 1) * (J - 1)), rep(1:(J - 1), each = K - 1)),
                             as.vector((exp(outer(ext_beta, mu + ext_alpha, FUN = "+") + ext_delta) * R)
                                       [-1, -1, drop = FALSE]))
  deriv_delta_beta <- "[<-"(matrix(0, (J - 1) * (K - 1), K - 1),
                            cbind(seq(1, (J - 1) * (K - 1)), rep(1:(K - 1), J - 1)),
                            as.vector((exp(outer(ext_beta, mu + ext_alpha, FUN = "+") + ext_delta) * R)
                                      [-1, -1, drop = FALSE]))
  if (!use_band) {
    index_weights_age <- cbind(
      1:((K - 1) * (J - 2)) + K - 1,
      1:((K - 1) * (J - 2))
    )
    index_weights_cohort <- cbind(
      (1:((J - 1) * (K - 1)))[-seq(1, (K - 1) * (J - 1), by = (K - 1))],
      (1:((J - 1) * (K - 1)))[-(K - 1)  * seq(1, J - 1)]
    )
    deriv_pen_diag <- diag(as.vector(
      weights_age + "["(cbind(weights_age, 0), 1:(K - 1), -1) +
        weights_cohort + weights_cohort %>% rbind(0) %>% "["(-1, 1:(J - 1))
    ))
    deriv_pen_weights <- "[<-"("[<-"(matrix(0, (K - 1) * (J - 1), (K - 1) * (J - 1)),
                                     index_weights_age, -as.vector(weights_age[, -1, drop = FALSE])),
                               index_weights_cohort, -as.vector(weights_cohort[-1, , drop = FALSE]))
    deriv_diag_delta <- diag(as.vector(
      (exp(outer(ext_beta, mu + ext_alpha, FUN = "+") + ext_delta) * R)[-1, -1, drop = FALSE]),
      (J - 1) * (K - 1), (J - 1) * (K - 1)) +
      pen * (deriv_pen_weights + t(deriv_pen_weights) + deriv_pen_diag)

    cbind(rbind(deriv_diag_mu, deriv_alpha_mu, deriv_beta_mu, deriv_delta_mu),
          rbind(t(deriv_alpha_mu), deriv_diag_alpha, deriv_beta_alpha, deriv_delta_alpha),
          rbind(t(deriv_beta_mu), t(deriv_beta_alpha), deriv_diag_beta, deriv_delta_beta),
          rbind(t(deriv_delta_mu), t(deriv_delta_alpha), t(deriv_delta_beta), deriv_diag_delta)) %>%
      "dimnames<-"(list(c("mu", rep("alpha", J - 1), rep("beta", K - 1), rep("delta", (J - 1) * (K - 1))),
                        c("mu", rep("alpha", J - 1), rep("beta", K - 1), rep("delta", (J - 1) * (K - 1)))))
  } else {
    index_cohort <- (1:((J - 1) * (K - 1)))[-(K - 1)  * seq(1, J - 1)]
    deriv_diag_delta <- cbind(
      pen * as.vector(weights_age + "["(cbind(weights_age, 0), 1:(K - 1), -1)) +
        pen * as.vector(weights_cohort + weights_cohort %>% rbind(0) %>% "["(-1, 1:(J - 1))) +
        as.vector((exp(outer(ext_beta, mu + ext_alpha, FUN = "+") + ext_delta) * R)[-1, -1, drop = FALSE]),
      -pen * "[<-"(rep(0, (K - 1) * (J - 1)), index_cohort, as.vector(weights_cohort[-1, , drop = FALSE])),
      matrix(rep(0, (K - 1) * (J - 1) * (K - 3)), (K - 1) * (J - 1), K - 3),
      -pen * "[<-"(rep(0, (K - 1) * (J - 1)), 1:((J - 2) * (K - 1)), as.vector(weights_age[, -1, drop = FALSE]))
    )
    list('A' = cbind(rbind(deriv_diag_mu, deriv_alpha_mu, deriv_beta_mu),
                     rbind(t(deriv_alpha_mu), deriv_diag_alpha, deriv_beta_alpha),
                     rbind(t(deriv_beta_mu), t(deriv_beta_alpha), deriv_diag_beta)) %>%
           'dimnames<-'(list(c("mu", rep("alpha", J - 1), rep("beta", K - 1)),
                             c("mu", rep("alpha", J - 1), rep("beta", K - 1)))),
         'B' = t(cbind(deriv_delta_mu, deriv_delta_alpha, deriv_delta_beta)) %>%
           'dimnames<-'(list(c("mu", rep("alpha", J - 1), rep("beta", K - 1)),
                             rep("delta", (J - 1) * (K - 1)))),
         'D' = deriv_diag_delta)
  }
}
hessian_aci_no_band <- function(par, O, R, pen, weights_age = NULL,
                                weights_cohort = NULL) {
  K <- nrow(O)
  J <- ncol(O)
  if (is.null(weights_age)) weights_age <- matrix(1, K - 1, J - 1)
  if (is.null(weights_cohort)) weights_cohort <- matrix(1, K - 1, J - 1)
  mu <- par[1]
  alpha <- par[2:J]
  ext_alpha <- c(0, alpha)
  beta <- par[(J + 1):(J + K - 1)]
  ext_beta <- c(0, beta)
  delta <- matrix(par[-(1:(J + K - 1))], K - 1, J - 1)
  ext_delta <- "[<-"(matrix(0, K, J), -1, -1, delta)
  deriv_diag_mu <- sum(exp(outer(ext_beta, mu + ext_alpha, FUN = "+") + ext_delta) * R)
  deriv_diag_alpha <- diag(sapply(2:J, function(ind_j) sum( exp(outer(ext_beta, mu + ext_alpha[ind_j], FUN = "+")
                                                                + ext_delta[, ind_j]) * R[, ind_j])), J - 1, J - 1)
  deriv_diag_beta <- diag(sapply(2:K, function(ind_k) sum( exp(outer(ext_beta[ind_k], mu + ext_alpha, FUN = "+")
                                                               + ext_delta[ind_k, ]) * R[ind_k, ])), K - 1, K - 1)
  deriv_alpha_mu <- matrix(sapply(2:J, function(ind_j) sum( exp(outer(ext_beta, mu + ext_alpha[ind_j], FUN = "+")
                                                                + ext_delta[, ind_j]) * R[, ind_j])), J - 1, 1)
  deriv_beta_mu <- matrix(sapply(2:K, function(ind_k) sum( exp(outer(ext_beta[ind_k], mu + ext_alpha, FUN = "+")
                                                               + ext_delta[ind_k,]) * R[ind_k,])), K - 1, 1)
  deriv_delta_mu <- matrix(as.vector((exp(outer(ext_beta, mu + ext_alpha, FUN = "+") + ext_delta) * R)[-1, -1]), (J - 1) * (K - 1), 1)
  deriv_beta_alpha <- (exp(outer(ext_beta, mu + ext_alpha, FUN = "+") + ext_delta) * R)[-1, -1]
  deriv_delta_alpha <- "[<-"(matrix(0, (J - 1) * (K - 1), J - 1),
                             cbind(1:((K - 1) * (J - 1)), rep(1:(J - 1), each = K - 1)),
                             as.vector((exp(outer(ext_beta, mu + ext_alpha, FUN = "+") + ext_delta) * R)[-1, -1]))
  deriv_delta_beta <- "[<-"(matrix(0, (J - 1) * (K - 1), K - 1),
                            cbind(seq(1, (J - 1) * (K - 1)), rep(1:(K - 1), J - 1)),
                            as.vector((exp(outer(ext_beta, mu + ext_alpha, FUN = "+") + ext_delta) * R)[-1, -1]))
  index_weights_age <- cbind(
    1:((K - 1) * (J - 2)) + K - 1,
    1:((K - 1) * (J - 2))
  )
  index_weights_cohort <- cbind(
    (1:((J - 1) * (K - 1)))[-seq(1, (K - 1) * (J - 1), by = (K - 1))],
    (1:((J - 1) * (K - 1)))[-(K - 1)  * seq(1, J - 1)]
  )
  deriv_pen_diag <- diag(as.vector(
    weights_age + "["(cbind(weights_age, 0), 1:(K - 1), -1) +
      weights_cohort + weights_cohort %>% rbind(0) %>% "["(-1, 1:(J - 1))
  ))
  deriv_pen_weights <- "[<-"("[<-"(matrix(0, (K - 1) * (J - 1), (K - 1) * (J - 1)),
                                   index_weights_age, -as.vector(weights_age[, -1])),
                             index_weights_cohort, -as.vector(weights_cohort[-1, ]))
  deriv_diag_delta <- diag(
    as.vector((exp(outer(ext_beta, mu + ext_alpha, FUN = "+") + ext_delta) * R)[-1, -1]),
    (J - 1) * (K - 1), (J - 1) * (K - 1)) + pen * (deriv_pen_weights + t(deriv_pen_weights) + deriv_pen_diag)

  cbind(rbind(deriv_diag_mu, deriv_alpha_mu, deriv_beta_mu, deriv_delta_mu),
        rbind(t(deriv_alpha_mu), deriv_diag_alpha, deriv_beta_alpha, deriv_delta_alpha),
        rbind(t(deriv_beta_mu), t(deriv_beta_alpha), deriv_diag_beta, deriv_delta_beta),
        rbind(t(deriv_delta_mu), t(deriv_delta_alpha), t(deriv_delta_beta), deriv_diag_delta)) %>%
    "dimnames<-"(list(c("mu", rep("alpha", J - 1), rep("beta", K - 1), rep("delta", (J - 1) * (K - 1))),
                      c("mu", rep("alpha", J - 1), rep("beta", K - 1), rep("delta", (J - 1) * (K - 1)))))
}
loglik_aci_sel_old <- function(par, O, R, sel, L) {
  K <- nrow(O)
  J <- ncol(O)
  if (length(par) != J + K - 1 + L) {
    stop("Error: selection nlevels not consistent with J * K")
  }
  mu <- par[1]
  ext_alpha <- c(0, par[2:J])
  ext_beta <- c(0, par[(J + 1):(J + K - 1)])
  z <- par[-(1:(J + K - 1))]
  delta <- matrix(NA, K - 1, J - 1)
  for (level_ind in 1:L) {
    delta[sel == level_ind] <- z[level_ind]
  }
  ext_delta <- "[<-"(matrix(0, K, J), -1, -1, delta)
  eta <- outer(ext_beta, mu + ext_alpha, FUN = "+") + ext_delta
  sum(exp(eta) * R - eta * O)
}
#' Foobar
#' @export
par2haz_aci_sel <- function(par, O, R, sel_array) {
  K <- nrow(O)
  J <- ncol(O)
  if (any(dim(sel_array)[1:2] != c(K - 1, J - 1))) {
    stop("Error : sel_array dimensions do not agree with exhaustive statistics")
  }
  # L <- dim(sel_array)[3]
  mu <- par[1]
  ext_alpha <- c(0, par[2:J])
  ext_beta <- c(0, par[(J + 1):(J + K - 1)])
  delta <- sweep(sel_array, MARGIN = 3, par[-(1:(J + K - 1))], '*') %>%
    apply(., MARGIN = c(1, 2), sum)
  ext_delta <- "[<-"(matrix(0, K, J), -1, -1, delta)
  exp(outer(ext_beta, mu + ext_alpha, FUN = "+") + ext_delta) %>%
    'dimnames<-'(dimnames(O))
}
#' Negative lok-likelihood in the age-cohort-interaction model with constaints
#'
## @family aci_sel_likelihood
#' @param par parameter of the age-cohort-interaction model in the
#' form \code{c(mu, alpha and beta)}
#' @param O Observed events as returned by \code{\link{exhaustive_stat_2d}}
#' @param R Time at risk as returned by \code{\link{exhaustive_stat_2d}}
#' @param sel_array an array defining the selection found by the adaptive ridge
#' procedure. (see \code{aridge_solver_aci})
#' @return The vector of derivatives of the negative log-likelihood as
#' described in TODO NAME OF REFERENCE SHEET
#' @examples
#' \dontrun{
#' par <- rnorm(ncol(O) * nrow(O))
#' loglik_aci_sel(par, O, R) # length 1 vector
#' score_aci_sel(par, O, R) # vector of length (ncol(O) + nrow(O))
#' hessian_aci_sel(par, O, R) # square matrix of dimension (ncol(O) + nrow(O))
#' }
#' @export
loglik_aci_sel <- function(par, O, R, sel_array) {
  if (any(dim(O) != dim(R))) stop('Error: dimensions of O, R, and par must agree')
  if (O < 0 || R < 0) stop('Error: O and R must have non-negative values')
  K <- nrow(O)
  J <- ncol(O)
  if (any(dim(sel_array)[1:2] != c(K - 1, J - 1))) {
    stop("Error : dimensions of O, R, and sel_array must agree")
  }
  if (length(par) != dim(sel_array)[3] + K + J - 1) {
    stop("Error: dimensions of par and sel_array must agree")
  }
  # L <- dim(sel_array)[3]
  mu <- par[1]
  ext_alpha <- c(0, par[2:J])
  ext_beta <- c(0, par[(J + 1):(J + K - 1)])
  delta <- sweep(sel_array, MARGIN = 3, par[-(1:(J + K - 1))], '*') %>%
    apply(., MARGIN = c(1, 2), sum)
  ext_delta <- "[<-"(matrix(0, K, J), -1, -1, delta)
  eta <- outer(ext_beta, mu + ext_alpha, FUN = "+") + ext_delta
  sum(exp(eta) * R - eta * O)
}
#' @rdname loglik_aci_sel
#' @export
score_aci_sel <- function(par, O, R, sel_array) {
  K <- nrow(O)
  J <- ncol(O)
  if (any(dim(sel_array)[1:2] != c(K - 1, J - 1))) {
    stop("Error : sel_array dimensions do not agree with exhaustive statistics")
  }
  L <- dim(sel_array)[3]
  mu <- par[1]
  ext_alpha <- c(0, par[2:J])
  ext_beta <- c(0, par[(J + 1):(J + K - 1)])
  delta <- sweep(sel_array, MARGIN = 3, par[-(1:(J + K - 1))], '*') %>%
    apply(., MARGIN = c(1, 2), sum)
  ext_delta <- "[<-"(matrix(0, K, J), -1, -1, delta)
  eta <- outer(ext_beta, mu + ext_alpha, FUN = "+") + ext_delta
  deriv_mu <- sum(exp(outer(ext_beta, mu + ext_alpha, FUN = "+") + ext_delta) * R - O)
  deriv_alpha <- sapply(2:J, function(ind_j) sum( exp(outer(ext_beta, mu + ext_alpha[ind_j], FUN = "+")
                                                      + ext_delta[, ind_j]) * R[, ind_j] - O[, ind_j]))
  deriv_beta <- sapply(2:K, function(ind_k) sum( exp(outer(ext_beta[ind_k], mu + ext_alpha, FUN = "+")
                                                     + ext_delta[ind_k, ]) * R[ind_k, ] - O[ind_k, ]))
  temp_array <- (array(replicate(L, exp(eta) * R - O)[-1, -1, ], dim = c(K - 1, J - 1, L)) * sel_array)
  deriv_delta <- apply(temp_array, MARGIN = 3, sum)
  c(deriv_mu, deriv_alpha, deriv_beta, deriv_delta) %>%
    setNames(c('mu', rep('alpha', J - 1), rep('beta', K - 1), rep('delta', L)))
}
#' @rdname loglik_aci_sel
#' @export
hessian_aci_sel <- function(par, O, R, sel_array) {
  K <- nrow(O)
  J <- ncol(O)
  if (any(dim(sel_array)[1:2] != c(K - 1, J - 1))) {
    stop("Error : sel_array dimensions do not agree with exhaustive statistics")
  }
  L <- dim(sel_array)[3]
  mu <- par[1]
  ext_alpha <- c(0, par[2:J])
  ext_beta <- c(0, par[(J + 1):(J + K - 1)])
  delta <- sweep(sel_array, MARGIN = 3, par[-(1:(J + K - 1))], '*') %>%
    apply(., MARGIN = c(1, 2), sum)
  ext_delta <- "[<-"(matrix(0, K, J), -1, -1, delta)
  eta <- outer(ext_beta, mu + ext_alpha, FUN = "+") + ext_delta
  deriv_diag_mu <- sum(exp(outer(ext_beta, mu + ext_alpha, FUN = "+") + ext_delta) * R)
  deriv_diag_alpha <- diag(sapply(2:J, function(ind_j) sum( exp(outer(ext_beta, mu + ext_alpha[ind_j], FUN = "+")
                                                                + ext_delta[, ind_j]) * R[, ind_j])), J - 1, J - 1)
  deriv_diag_beta <- diag(sapply(2:K, function(ind_k) sum( exp(outer(ext_beta[ind_k], mu + ext_alpha, FUN = "+")
                                                               + ext_delta[ind_k, ]) * R[ind_k, ])), K - 1, K - 1)
  deriv_alpha_mu <- matrix(sapply(2:J, function(ind_j) sum( exp(outer(ext_beta, mu + ext_alpha[ind_j], FUN = "+")
                                                                + ext_delta[, ind_j]) * R[, ind_j])), J - 1, 1)
  deriv_beta_mu <- matrix(sapply(2:K, function(ind_k) sum( exp(outer(ext_beta[ind_k], mu + ext_alpha, FUN = "+")
                                                               + ext_delta[ind_k,]) * R[ind_k,])), K - 1, 1)
  temp_array <- (array(replicate(L, exp(eta) * R)[-1, -1, ], dim = c(K - 1, J - 1, L)) * sel_array)
  deriv_delta_mu <- temp_array %>%
    apply(., MARGIN = 3, sum) %>%
    matrix(., L, 1)
  deriv_beta_alpha <- (exp(outer(ext_beta, mu + ext_alpha, FUN = "+") + ext_delta) * R)[-1, -1]
  deriv_delta_alpha <- temp_array %>%
    apply(., MARGIN = c(2, 3), sum) %>%
    t()
  deriv_delta_beta <- temp_array %>%
    apply(., MARGIN = c(1, 3), sum) %>%
    t()
  deriv_diag_delta <- temp_array %>%
    apply(., MARGIN = 3, sum) %>%
    diag(., nrow = length(.))
  cbind(rbind(deriv_diag_mu, deriv_alpha_mu, deriv_beta_mu, deriv_delta_mu),
        rbind(t(deriv_alpha_mu), deriv_diag_alpha, deriv_beta_alpha, deriv_delta_alpha),
        rbind(t(deriv_beta_mu), t(deriv_beta_alpha), deriv_diag_beta, deriv_delta_beta),
        rbind(t(deriv_delta_mu), t(deriv_delta_alpha), t(deriv_delta_beta), deriv_diag_delta)) %>%
    "dimnames<-"(list(c("mu", rep("alpha", J - 1), rep("beta", K - 1), rep("delta", L)),
                      c("mu", rep("alpha", J - 1), rep("beta", K - 1), rep("delta", L))))
}
#' Computes the estimate of the constrained aci model
#'
#' @param O Observed events as returned by \code{\link{exhaustive_stat_2d}}
#' @param R Time at risk as returned by \code{\link{exhaustive_stat_2d}}
#' @param sel_array an array defining the selection found by the adaptive ridge
#' procedure. (see \code{aridge_solver_aci})
#' @return The vector estimate of the constrained age-cohort-interaction model.
#' The \code{sel_array} is an array of dimension \code{K  J  L}
#' @export
ridge_solver_aci_sel <- function(O, R, sel_array, maxiter = 1000) {
  K <- nrow(O)
  J <- ncol(O)
  L <- dim(sel_array)[3]
  old_par <- rep(0, 1 + J - 1 + K - 1 + L)
  for (iter in 1:maxiter) {
    par <- old_par - Solve(hessian_aci_sel(old_par, O, R, sel_array),
                           score_aci_sel(old_par, O, R, sel_array))
    if (sum(is.na(abs(par - old_par) / abs(old_par)))) break
    if (max(abs(par - old_par)) <= 1e-8) break
    old_par <- par
  }
  if (iter == maxiter) {
    warning("Warning: Newton-Raphson procedure did not converge")
  }
  list("par" = par, "convergence" = iter != maxiter, "niter" = iter)
}
ridge_solver_aci_old <- function(O, R, pen, use_band = FALSE, maxiter = 1000) {
  old_par <- rep(0, ncol(O) * nrow(O))
  for (iter in 1:maxiter) {
    if (use_band) {
      par <- old_par - block_bandsolve_rotated(hessian_aci(old_par, O, R, pen, use_band = TRUE),
                                               score_aci(old_par, O, R, pen))
    } else {
      par <- old_par - Solve(hessian_aci(old_par, O, R, pen),
                             score_aci(old_par, O, R, pen))
    }
    if (sum(is.na(abs(par - old_par) / abs(old_par)))) break
    if (max(abs(par - old_par)/abs(old_par)) <= sqrt(.Machine$double.eps)) break
    old_par <- par
  }
  if (iter == maxiter) {
    warning("Warning: Newton-Raphson procedure did not converge")
  }
  list("par" = par, "convergence" = iter != maxiter, "niter" = iter)
}
#' Computes the ridge regularized estimates of the age-cohort-interaction
#'  model
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
#' is recommended
#' @return The vector estimate of the ridge regularized age-cohort-interaction model.
#' @examples
#' \dontrun{
#' J <- 10
#' K <- 15
#' set.seed(0)
#' O <- matrix(rpois(K * J, 2), K, J)
#' R <- matrix(rpois(K * J, 10), K, J)
#' pen <-  50
#' ridge <- ridge_solver_aci(O, R, pen)
#' ridge$par
#' ridge$iter
#' }
#' @export
#' @family ridge
ridge_solver_aci <- function(O, R, pen, weights_age = NULL,
                             weights_cohort = NULL, old_par = NULL,
                             use_band = FALSE, maxiter = 100) {
  if (is.null(old_par)) old_par <- rep(0, ncol(O) * nrow(O))
  for (iter in 1:maxiter) {
    if (use_band) {
      par <- old_par - block_bandsolve(
        hessian_aci(old_par, O, R, pen, weights_age, weights_cohort, use_band = TRUE),
        score_aci(old_par, O, R, pen, weights_age, weights_cohort),
        mat_as_rotated = TRUE,
        lim = NULL)
    } else {
      par <- old_par - Solve(
        hessian_aci(old_par, O, R, pen, weights_age, weights_cohort),
        score_aci(old_par, O, R, pen, weights_age, weights_cohort))
    }
    if (sum(is.na(abs(par - old_par) / abs(old_par)))) break
    if (max(abs(par - old_par)) <= 1e-4) break
    old_par <- par
  }
  if (iter == maxiter) {
    warning("Warning: Newton-Raphson procedure did not converge")
  }
  list("par" = par, "convergence" = iter != maxiter, "niter" = iter)
}
#' @rdname ridge_solver_aci
#' @export
#' @family adaptive_ridge
aridge_solver_aci <- function(O, R, pen, sample_size,
                              use_band = FALSE,
                              maxiter = 1000,
                              verbose = FALSE) {
  pb <- progress_bar$new(
    format = "  adaptive ridge [:bar] :percent in :elapsed",
    total = length(pen), clear = FALSE, width = 100)
  pb$tick(0)
  sel <- par_sel <- haz_sel <- vector('list', length(pen))
  bic <- aic <- ebic <- rep(NA, length(pen))
  epsilon_age <- 1e-5
  epsilon_cohort <- 1e-5
  K <- nrow(O)
  J <- ncol(O)
  weights_age <- matrix(1, K - 1, J - 1)
  weights_cohort <- matrix(1, K - 1, J - 1)
  valve_age <- matrix(1, K - 1, J - 1) %>%
    "rownames<-"(diag(sapply(rownames(O)[-1], paste0, "-", rownames(O)[-J]))) %>%
    "colnames<-"(diag(sapply(colnames(O)[-1], paste0, "-", colnames(O)[-J])))
  valve_cohort <- matrix(1, K - 1, J - 1) %>%
    "rownames<-"(diag(sapply(rownames(O)[-1], paste0, "-", rownames(O)[-J]))) %>%
    "colnames<-"(diag(sapply(colnames(O)[-1], paste0, "-", colnames(O)[-J])))
  old_par <- rep(0, J * K)
  ind_pen <- 1
  for (iter in 1:maxiter) {
    if (verbose) cat("iter =", iter, '\n')
    old_valve_age <- valve_age
    old_valve_cohort <- valve_cohort
    par <- ridge_solver_aci(O, R, pen[ind_pen],
                            weights_age, weights_cohort,
                            old_par, use_band)$par
    delta <- matrix(par[-(1:(J + K - 1))], K - 1, J - 1)
    weights_age[, ] <- 1 / (t(apply(cbind(0, delta), 1, diff)) ^ 2 + epsilon_age ^ 2)
    weights_cohort[, ] <- 1 / (apply(rbind(0, delta), 2, diff) ^ 2 + epsilon_cohort ^ 2)
    valve_age[, ] <- (weights_age * t(apply(cbind(0, delta), 1, diff)) ^ 2)
    valve_cohort[, ] <- (weights_cohort * apply(rbind(0, delta), 2, diff) ^ 2)
    converge <- max(abs(old_valve_age - valve_age),
                      abs(old_valve_cohort - valve_cohort)) <= 1e-6
    if (converge) {
      selection <- valve2sel_aci(valve_age, valve_cohort)
      sel[[ind_pen]] <- selection$fct
      if (verbose) sel[[ind_pen]] %>% raster('factor')
      sel_array <- selection$array
      L <- dim(sel_array)[3]
      par_sel[[ind_pen]] <- ridge_solver_aci_sel(O, R, sel_array)$par
      haz_sel[[ind_pen]] <- par2haz_aci_sel(par_sel[[ind_pen]], O, R, sel_array)
      aic[ind_pen] <- 2 * L + 2 * loglik_aci_sel(par_sel[[ind_pen]], O, R, sel_array)
      bic[ind_pen] <- log(sample_size) * L +
        2 * loglik_aci_sel(par_sel[[ind_pen]], O, R, sel_array)
      ebic[ind_pen] <- bic[ind_pen] +  2 * log(choose(J * K, L))
      pb$tick()
      ind_pen <-  ind_pen + 1
    }
    old_par <- par
    if (ind_pen == length(pen) + 1) break
  }
  if (iter == maxiter) {
    warning("Warning: aridge did not converge")
  }
  return(list("par" = par_sel, "haz" = haz_sel, "sel" = sel,
              "bic" = bic, "aic" = aic, 'ebic' = ebic))
}
#' @export
aridge_solver_aci_q <- function(O, R, pen, sample_size,
                              use_band = FALSE, q = 2,
                              maxiter = 1000,
                              verbose = FALSE) {
  pb <- progress_bar$new(
    format = "  adaptive ridge [:bar] :percent in :elapsed",
    total = length(pen), clear = FALSE, width = 100)
  pb$tick(0)
  sel <- par_sel <- haz_sel <- vector('list', length(pen))
  bic <- aic <- ebic <- rep(NA, length(pen))
  epsilon_age <- 1e-5
  epsilon_cohort <- 1e-5
  K <- nrow(O)
  J <- ncol(O)
  weights_age <- matrix(1, K - 1, J - 1)
  weights_cohort <- matrix(1, K - 1, J - 1)
  valve_age <- matrix(1, K - 1, J - 1) %>%
    "colnames<-"(diag(sapply(colnames(O)[-1], paste0, "-", colnames(O)[-J]))) %>%
    "rownames<-"(diag(sapply(rownames(O)[-1], paste0, "-", rownames(O)[-J])))
  valve_cohort <- matrix(1, K - 1, J - 1) %>%
    "rownames<-"(diag(sapply(rownames(O)[-1], paste0, "-", rownames(O)[-J]))) %>%
    "colnames<-"(diag(sapply(colnames(O)[-1], paste0, "-", colnames(O)[-J])))
  old_par <- rep(0, J * K)
  ind_pen <- 1
  for (iter in 1:maxiter) {
    if (verbose) cat("iter =", iter, '\n')
    old_valve_age <- valve_age
    old_valve_cohort <- valve_cohort
    par <- ridge_solver_aci(O, R, pen[ind_pen],
                            weights_age, weights_cohort,
                            old_par, use_band)$par
    delta <- matrix(par[-(1:(J + K - 1))], K - 1, J - 1)
    weights_age[, ] <- 1 / (t(apply(cbind(0, delta), 1, diff)) ^ 2 + epsilon_age ^ 2) ^ (1 - q / 2)
    weights_cohort[, ] <- 1 / (apply(rbind(0, delta), 2, diff) ^ 2 + epsilon_cohort ^ 2) ^ (1 - q / 2)
    valve_age[, ] <- (weights_age * t(apply(cbind(0, delta), 1, diff)) ^ 2)
    valve_cohort[, ] <- (weights_cohort * apply(rbind(0, delta), 2, diff) ^ 2)
    converged2 <- max(abs(old_valve_age - valve_age),
                      abs(old_valve_cohort - valve_cohort)) <= 1e-6
    # old_par <- par
    if (converged2) {
      selection <- valve2sel_aci(valve_age, valve_cohort)
      sel[[ind_pen]] <- selection$fct
      if (verbose) sel[[ind_pen]] %>% raster('factor')
      sel_array <- selection$array
      L <- dim(sel_array)[3]
      par_sel[[ind_pen]] <- ridge_solver_aci_sel(O, R, sel_array)$par
      haz_sel[[ind_pen]] <- par2haz_aci_sel(par_sel[[ind_pen]], O, R, sel_array)
      aic[ind_pen] <- 2 * L + 2 * loglik_aci_sel(par_sel[[ind_pen]], O, R, sel_array)
      bic[ind_pen] <- log(sample_size) * L +
        2 * loglik_aci_sel(par_sel[[ind_pen]], O, R, sel_array)
      ebic[ind_pen] <- bic[ind_pen] +  2 * log(choose(J * K, L))
      # cat('progress:', ind_pen / length(pen) * 100, '% \n')
      pb$tick()
      ind_pen <-  ind_pen + 1
      # par <- par * 0
      weights_age[] <- 1
      weights_cohort[] <- 1
      valve_age[] <- 1
      valve_cohort[] <- 1
    }
    old_par <- par
    if (ind_pen == length(pen) + 1) break
  }
  if (iter == maxiter) {
    warning("Warning: aridge did not converge")
  }
  return(list("par" = par_sel, "haz" = haz_sel, "sel" = sel,
              "bic" = bic, "aic" = aic, 'ebic' = ebic))
}

