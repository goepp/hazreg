#' Generate geometric sequences
#'
#' @param x A strictly positive number
#' @param y A strictly positive number
#' @param length.out The length of the sequence. It must be a strictly positive integer
#' @return A geometric sequence
#' @examples
#' seq_geom(1, 10, 10)
#' seq_geom(1, 10, 10)
seq_geom <- function(a, b, length.out) {
  return(a * (b / a) ^ ((0:(length.out - 1))/(length.out - 1)))
}
reverselog_trans <- function(base = exp(1)) {
  library(scales)
  trans <- function(x) -log(x, base)
  inv <- function(x) base ^ (-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv,
            log_breaks(base = base),
            domain = c(1e-100, Inf))
}
#' Create two dimensional hazard rate from Age-Cohort-Interaction model
#'
#' @param mu Constant part of the log-hazard
#' @param age_coef Age-effect parameter sequence
#' @param cohort_coef Cohort-effect parameter sequence
#' @param age_seq
#' @param cohort_seq
#' @param islet_ls List of parameters of the interaction terms.
#' @return
#' The interactions terms are densities of two dimensional normal distributions.
#' \code{islet_ls} gives the mean, variance-covariance matrix and multiplicative coefficients of each
#' interaction term.
#' @examples
#' library(lattice)
#' mu <- log(1e-2)
#' age_seq <- seq(0, 100, length.out = 100)
#' cohort_seq <- seq(1900, 1900 + age_seq[101], length.out = age_length + 1)
#' age_coef <- seq(0, 2, length.out = age_length - 1)
#' cohort_coef <- seq(0, 0.8, length.out = age_length - 1)
#' islet_ls = list(list(mean = c(1945, 45), sigma = 100 * diag(2), coef = 40))
#' hazard <- map_maker(mu, age_coef, cohort_coef, age_seq, cohort_seq, islet_ls = islet_ls)
#' persp(hazard$sheet)
map_maker <- function(mu, age_coef, cohort_coef, age_seq, cohort_seq, islet_ls = NULL) {
  library(pryr)
  library(mvtnorm)
  ac.map <- exp(outer(c(0, cohort_coef), c(0, age_coef) + mu, FUN = "+"))
  colnames(ac.map) <- levels(cut(-1, age_seq, right = FALSE, dig.lab = 3))
  rownames(ac.map) <- levels(cut(-1, cohort_seq, right = FALSE, dig.lab = 4))
  islet.map <- 0 * ac.map
  if (!is.null(islet_ls)) {
    for (ind in 1:length(islet_ls)) {
      myFun <- function(ind1, ind2) dmvnorm(c(ind1, ind2), mean = islet_ls[[ind]]$mean, sigma = islet_ls[[ind]]$sigma)
      myVecFun <- Vectorize(myFun, vectorize.args = c('ind1', 'ind2'))
      islet.map <- islet.map +
        outer(cohort_seq[-length(cohort_seq)], age_seq[-length(age_seq)], FUN = myVecFun) *
        islet_ls[[ind]]$coef
    }
  }
  sheet <- ac.map + islet.map
  list("sheet" = sheet, "cohort" = cohort_seq, "age" = age_seq)
}
#' Generate data from piecewise constant hazard
#'
#' @param n number of observations. A positive integer
#' @param cuts cuts of the piecewise constant hazard
#' @param alpha values of the piecewise constant hazard
#' @return A sequence of observations of length \code{n} generated from the
#' distribution given by \code{cuts} and \code{alpha}
#' @examples
#' map_maker(100, c(20, 40, 60), c(0.01, 0.03, 0.05, 0.08))
rsurv <- function(n, cuts, alpha) {
  u <- runif(n)
  k <- length(alpha)
  if (length(cuts) != (k - 1)) stop("error: length(cuts) must be equal to length(alpha)-1 !")
  cuts0 <- c(0, cuts)
  if (k > 1) {
    thresh <- exp(-cumsum(alpha[-k] * diff(cuts0)))
    if (n <= 200) {
      seg <- apply(matrix(rep(thresh, n), byrow = T, nrow = n) > u, 1, sum) + 1
    } else {
      seg <- rep(NA, n)
      for (i in 1:n)
        seg[i] <- sum(thresh > u[i]) + 1
    }
  } else {
    seg <- rep(1, n)
  }
  cuts0[seg] - (log(u) + cumsum(c(0, alpha[-k] * diff(cuts0)))[seg]) / alpha[seg]
}
#' Generate age and cohort observations
#'
#' @param map A matrix representing the hazard, as provided by \code{\link{map_maker}}.
#' @param sample_size Number of observations
#' @param dob_dist Marginal distribution of the cohort. Defaults to the uniform distribution
#' with bounds given by \code{dob_dist_par}.
#' @param dob_dist_par Parameter of the distribution of the cohort
#' @return \code{surv_data}, a data frame of observations with age and cohort.
#' @examples
#' mu <- log(1e-2)
#' age_seq <- seq(0, 100, length.out = 100)
#' cohort_seq <- seq(1900, 1900 + age_seq[101], length.out = age_length + 1)
#' age_coef <- seq(0, 2, length.out = age_length - 1)
#' cohort_coef <- seq(0, 0.8, length.out = age_length - 1)
#' islet_ls = list(list(mean = c(1945, 45), sigma = 100 * diag(2), coef = 40))
#' hazard <- map_maker(mu, age_coef, cohort_coef, age_seq, cohort_seq, islet_ls = islet_ls)
#' (map_surv(hazard, 10))
map_surv <- function(map, sample_size, dob_dist, dob_dist_par = NULL) {
  require(pryr)
  if (is.character(dob_dist)) {
    if (is.null(dob_dist_par)) error("error: no parameter supplied for dob distribution")
    distf = paste0("r", dob_dist)
    if (length(formals(fget(distf))) - 1 != length(dob_dist_par)) error("error: too few parameters supplied for this distribution")
    cohort_sample <- do.call(fget(distf), # simulated dates of birth
                             as.list(c(sample_size, dob_dist_par[1:(length(formals(fget(distf))) - 1)])))
  } else {
    cohort_sample <- sample(dob_dist, sample_size, replace = T)
  }
  cohort_splitting <- split(cohort_sample, cut(cohort_sample, map$cohort, dig.lab = 4, right = FALSE))
  age_splitting <- vector("list", length(map$cohort) - 1)
  for (ind in seq_along(age_splitting)) {
    age_splitting[[ind]] <- rsurv(length(cohort_splitting[[ind]]), map$age[2:(length(map$age) - 1)], unname(map$sheet[ind, ]))
  }
  age_sample <- unsplit(age_splitting, cut(cohort_sample, map$cohort, dig.lab = 4, right = FALSE))
  data.frame("age" = age_sample, "cohort" = cohort_sample)
}
exhaustive_stat <- function(surv.data, cuts_age, cuts_cohort) {
  ext_cuts_age <- c(0, cuts_age)
  ext_cuts_cohort <- c(min(surv.data$cohort), cuts_cohort)
  temp1 <- model.matrix(~0 + cut(surv.data$age, breaks = c(ext_cuts_age, Inf)))
  temp2 <- t(1 - apply(temp1, 1, cummax)) + temp1
  temp3 <- t(apply(apply(
    cbind(t(matrix(rep(ext_cuts_age, length(surv.data$age)), ncol = length(surv.data$age))),
          surv.data$age), 1, sort),
    2, diff))
  R_old <- as.data.frame(unname(temp2 * temp3))
  rearrange_by_cohort <- function(vect) split(vect, cut(
    surv.data$cohort, breaks = c(ext_cuts_cohort, Inf),right = FALSE, dig.lab = 4))
  sum_by_cohort <- function(vect) unlist(sapply(rearrange_by_cohort(vect), sum))
  R <- sapply(R_old, sum_by_cohort)
  O_old <- as.data.frame(unname(temp1 * matrix(rep(surv.data$delta, length(ext_cuts_age)), ncol = length(ext_cuts_age))))
  O <- sapply(O_old, sum_by_cohort)
  colnames(O) <- colnames(R) <- levels(cut(0, breaks = c(ext_cuts_age, Inf), right = FALSE, dig.lab = 3))
  list("O" = O, "R" = R)
}
exhaustive_stat_sel <- function(exhaust, sel) {
  if (any(dim(exhaust$O) != dim(sel))) stop("error: dimensions of exhaust and selection must agree")
  O <- lapply(levels(as.factor(sel)),
              function(sel_val) sum(exhaust$O[which(sel == sel_val, arr.ind = TRUE)])) %>%
    unlist()
  R <- lapply(levels(as.factor(sel)),
              function(sel_val) sum(exhaust$R[which(sel == sel_val, arr.ind = TRUE)])) %>%
    unlist()
  return(list("O" = O, "R" = R))
}
#' Solve linear system when the matrix is block-wise with a band diagonal block
#'
#' @param hessian_mat Square matrix
#' @param score_vect Vector of same length as \code{hessian_mat}
#' @param lim Positive integer
#' @return The vector \code{hessian_mat} ^ {-1} \code{score_vect} is returned.
#' If \code{hessian_mat[lim:end, lim:end]} is not a band matrix or if \code{hessian_mat} is
#' symmetric, an error is returned.
#' If these two conditions are met, the returned value is the same as
#' \code{solve(hessian_mat, score_vect)}. There might be nummerical impercision errors.
#' @examples
#' set.seed(0)
#' A <- matrix(rnorm(36), 6, 6)
#' B <- matrix(rnorm(24), 6, 4)
#' C <- diag(rnorm(4)) + rbind(cbind(0, diag(rnorm(3))), 0)
#' hessian_mat <- rbind(cbind(A + t(A), B),
#'                      cbind(t(B), C + t(C)))
#' score_vect <- rnorm(10)
#' lim <- 6
#' max(band_prod(hessian_mat, score_vect, lim) - solve(hessian_mat, score_vect))
band_prod <- function(hessian_mat, score_vect, lim) {
  if (!isSymmetric(hessian_mat) || ncol(hessian_mat) != length(as.vector(score_vect))) {
    error("error: dimensions of hessian and score do not match")
  }
  S1 <- score_vect[1:lim]
  S2 <- score_vect[(lim + 1):(ncol(hessian_mat))]
  A <- hessian_mat[1:lim, 1:lim]
  B <- hessian_mat[1:lim, (lim + 1):ncol(hessian_mat)]
  D <- hessian_mat[(lim + 1):nrow(hessian_mat), (lim + 1):ncol(hessian_mat)]
  rotD <- mat2rot(D)
  schur <- A - B %*% bandsolve(rotD, t(B))
  temp1 <- Solve(schur, S1)
  temp2 <- Solve(schur, B %*% bandsolve(rotD, S2))
  c(as.vector(temp1 - temp2),
    as.vector(-bandsolve(rotD, t(B) %*% temp1) + bandsolve(rotD, S2) +
                bandsolve(rotD, t(B) %*% temp2))
  )
}
band_prod_rotated <- function(hessian_mat, score_vect) {
  if (nrow(hessian_mat$A) != ncol(hessian_mat$A)) {
    stop("Error: upper-left matrix not square")
  }
  if (max(hessian_mat$A - t(hessian_mat$A)) >= 1e-14) {
    stop("Error: upper-left matrix not symmetric")
  }
  if ((nrow(hessian_mat$A) + nrow(hessian_mat$D)) != length(score_vect)) {
    stop('Error: dimensions of matrix and vector must agree')
  }
  S1 <- score_vect[1:ncol(hessian_mat$A)]
  S2 <- score_vect[-seq_along(S1)]
  # A <- hessian_mat$A
  # B <- hessian_mat$B
  # rotD <- hessian_mat$D
  schur <- hessian_mat$A - hessian_mat$B %*% bandsolve(hessian_mat$D, t(hessian_mat$B))
  temp1 <- Solve(schur, S1)
  temp2 <- Solve(schur, hessian_mat$B %*% bandsolve(hessian_mat$D, S2))
  c(as.vector(temp1 - temp2),
    as.vector(-bandsolve(hessian_mat$D, t(hessian_mat$B) %*% temp1) +
                bandsolve(hessian_mat$D, S2) +
                bandsolve(hessian_mat$D, t(hessian_mat$B) %*% temp2))
  )
}
valve2sel <- function(valve_age, valve_cohort, epsilon = 1e-8) {
  library(igraph)
  if (any(dim(valve_age) + c(0, 1) != dim(valve_cohort) + c(1, 0))) {
    stop("Error: dimensions of valve matrices must agree")
  }
  K <- nrow(valve_age)
  J <- ncol(valve_age) + 1
  adjacency <- diag(J * K);
  node_names <- matrix(1:(J * K), K, J)
  for (j in 1:J) {
    for (k in 1:K) {
      if (k > 1) {
        if (valve_cohort[k - 1, j] < epsilon) {
          adjacency[node_names[k, j], node_names[k - 1, j]] <-
            adjacency[node_names[k - 1, j], node_names[k, j]] <- 1
        }
      }
      if (j > 1) {
        if (valve_age[k, j - 1] < epsilon) {
          adjacency[node_names[k, j], node_names[k, j - 1]] <-
            adjacency[node_names[k, j - 1], node_names[k, j]] <- 1
        }
      }
    }
  }
  graph <- graph_from_adjacency_matrix(adjacency, mode = "undirected")
  (matrix(clusters(graph)$membership, K, J) %>%
      'rownames<-'(rownames(valve_age)) %>%
      'colnames<-'(colnames(valve_cohort)))
}
grid2raster <- function(grid_df, title = NULL) {
  ggplot(grid_df, aes(cohort, age, fill = value)) + geom_raster() +
    ggtitle(title) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
}
sel2segment <- function(sel, cuts_age, cuts_cohort) {
  library(dplyr)
  if (sel %>% as.factor() %>% nlevels() == 1) {
    return(tibble())
  } else {
    diff_age <- apply(sel, 1, diff) %>%
      t() %>%
      "rownames<-"(rownames(sel)) %>%
      "colnames<-"(cuts_age)
    diff_cohort <- apply(sel, 2, diff) %>%
      "rownames<-"(cuts_cohort) %>%
      "colnames<-"(colnames(sel))
    diff_age_df <- diff_age %>%
      melt(varnames = c("cohort", "age")) %>%
      filter(value != 0) %>%
      mutate(x = mapvalues(cohort, rownames(sel), seq(0, length(cuts_cohort)) + 0.5, warn_missing = FALSE)) %>%
      mutate(y = mapvalues(age, cuts_age, seq(1, length(cuts_age)) + 0.5, warn_missing = FALSE)) %>%
      select(-c(cohort, age, value)) %>%
      mutate_all(funs(as.character)) %>%
      mutate_all(funs(as.numeric)) %>%
      mutate(xend = x + 1) %>%
      mutate(yend = y)
    diff_cohort_df <- diff_cohort %>%
      melt(varnames = c('cohort', 'age')) %>%
      filter(value != 0) %>%
      mutate(x = mapvalues(cohort, cuts_cohort, seq(1, length(cuts_cohort)) + 0.5, warn_missing = FALSE)) %>%
      mutate(y = mapvalues(age, colnames(sel), seq(0, length(cuts_age)) + 0.5, warn_missing = FALSE)) %>%
      select(-c(cohort, age, value)) %>%
      mutate_all(funs(as.character)) %>%
      mutate_all(funs(as.numeric)) %>%
      mutate(xend = x) %>%
      mutate(yend = y + 1)
    return(bind_rows(diff_age_df, diff_cohort_df) %>% dplyr::as_data_frame())
  }
}
#' Compute exhaustive statistics for the piecewise constant hazard model
#'
#' @param surv_data data frame containing individuals' age and delta
#' @param cuts_age Cuts defining the age intervals of the piecewise constant hazard model
#' @return A list containing O, the observed events in each interval, and R, the time
#' at risk in each interval.
#' @examples
#' mu <- log(1e-2)
#' age_seq <- seq(0, 100, length.out = 100)
#' cohort_seq <- seq(1900, 1900 + age_seq[101], length.out = age_length + 1)
#' age_coef <- seq(0, 2, length.out = age_length - 1)
#' cohort_coef <- seq(0, 0.8, length.out = age_length - 1)
#' islet_ls = list(list(mean = c(1945, 45), sigma = 100 * diag(2), coef = 40))
#' hazard <- map_maker(mu, age_coef, cohort_coef, age_seq, cohort_seq, islet_ls = islet_ls)
#' set.seed(0)
#' sample_size <- 1000
#' censoring <- runif(sample_size, 75, 100)
#' surv_data <- map_surv(hazard, sample_size = sample_size,
#'                       dob_dist = "unif", dob_dist_par = c(1900, 2000)) %>%
#'                         mutate(delta = as.numeric(censoring >= age)) %>%
#'                           mutate(age = pmin(age, censoring))
#' surv_data
exhaustive_stat_1d <- function(surv_data, cuts_age) {
  ext_cuts_age <- c(0, cuts_age)
  temp1 <- model.matrix(~0 + cut(surv_data$age,
                                 breaks = c(ext_cuts_age, Inf),
                                 labels = 1:length(ext_cuts_age),
                                 include.lowest = TRUE))
  attr(temp1, "contrasts") <- NULL
  attr(temp1, "assign") <- NULL
  temp2 <- t(1 - apply(temp1, 1, cummax)) + temp1
  temp3 <- t(apply(apply(cbind(t(matrix(rep(ext_cuts_age, length(surv_data$age)),
                                        ncol = length(surv_data$age))),
                               surv_data$age), 1, sort), 2, diff))
  O <- temp1 * matrix(rep(surv_data$delta, length(ext_cuts_age)), ncol = length(ext_cuts_age))
  R <- temp3 * temp2
  colnames(O) <- colnames(R) <- levels(cut(0, breaks = c(ext_cuts_age, Inf), right = FALSE, dig.lab = 3))
  list(O = O %>% colSums(), R = R %>% colSums())
}
#' Compute exhaustive statistics for the piecewise constant hazard model
#'
#' @param surv_data data frame containing individuals' age, cohort, and delta
#' @param cuts_age Cuts defining the age intervals
#' @param cuts_cohort Cuts defining the cohort intervals
#' @return A list containing O, the observed events in each interval, and R, the time
#' at risk in each age-cohort rectangle
#' @examples
#' mu <- log(1e-2)
#' age_seq <- seq(0, 100, length.out = 100)
#' cohort_seq <- seq(1900, 1900 + age_seq[101], length.out = age_length + 1)
#' age_coef <- seq(0, 2, length.out = age_length - 1)
#' cohort_coef <- seq(0, 0.8, length.out = age_length - 1)
#' islet_ls = list(list(mean = c(1945, 45), sigma = 100 * diag(2), coef = 40))
#' hazard <- map_maker(mu, age_coef, cohort_coef, age_seq, cohort_seq, islet_ls = islet_ls)
#' set.seed(0)
#' sample_size <- 1000
#' censoring <- runif(sample_size, 75, 100)
#' surv_data <- map_surv(hazard, sample_size = sample_size,
#'                       dob_dist = "unif", dob_dist_par = c(1900, 2000)) %>%
#'                         mutate(delta = as.numeric(censoring >= age)) %>%
#'                           mutate(age = pmin(age, censoring))
#' surv_data
exhaustive_stat_2d <- function(surv_data, cuts_age, cuts_cohort) {
  surv_data <- surv_data %>%
    mutate(cohort_lvl = cut(surv_data$cohort,
                            breaks = c(-Inf, cuts_cohort, Inf),
                            include.lowest = TRUE)) # compute to which cohort interval an individual belongs
  surv_data_ls <- split(surv_data, surv_data$cohort_lvl) # split the data following the cohort intervals
  res <- lapply(surv_data_ls,
                exhaustive_stat_1d,
                cuts_age = cuts_age) # compute exhaustive stats for each cohort interval
  O <- lapply(res, function(x) x$O) %>% do.call(rbind, .) # retrive O and merge it into one matrix
  R <- lapply(res, function(x) x$R) %>% do.call(rbind, .) # retrive R and merge it into one matrix
  list(O = O, R = R)
}
###
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
###
par2grid_aci <- function(par, cuts_age, cuts_cohort) {
  J <- length(cuts_age) + 1
  K <- length(cuts_cohort) + 1
  par_df <- par2haz_aci(unname(par), J, K) %>%
    melt(varnames = c("cohort", "age")) %>%
    mutate(age = levels(cut(-1, breaks = c(0, cuts_age, Inf), right = FALSE, dig.lab = 3))[matrix(age)]) %>%
    mutate(age = factor(age, levels = unique(age))) %>%
    mutate(cohort = levels(cut(-1, breaks = c(0, cuts_cohort, Inf), right = FALSE, dig.lab = 4))[matrix(cohort)]) %>%
    mutate(cohort = factor(cohort, levels = unique(cohort)))
  par_df
}
valve2sel_aci <- function(valve_age, valve_cohort, epsilon = 1e-8) {
  library(igraph)
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
  (matrix(clusters(graph)$membership, K - 1, J - 1) %>%
      'rownames<-'(rownames(valve_age)) %>%
      'colnames<-'(colnames(valve_cohort)))
}
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
haz2grid <- function(haz, cuts_age, cuts_cohort) {
  haz %>%
    unname() %>%
    melt(varnames = c("cohort", "age")) %>%
    mutate(age = levels(cut(-1, breaks = c(0, cuts_age, Inf),
                            right = FALSE, dig.lab = 3))[matrix(age)]) %>%
    mutate(age = factor(age, levels = unique(age))) %>%
    mutate(cohort = levels(cut(-1, breaks = c(0, cuts_cohort, Inf),
                               right = FALSE, dig.lab = 4))[matrix(cohort)]) %>%
    mutate(cohort = factor(cohort, levels = unique(cohort)))
}
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
loglik_sel_aci_old <- function(par, O, R, sel) {
  K <- nrow(O)
  J <- ncol(O)
  eta <- par2haz_sel(par, sel, J, K, haz.log = TRUE)
  sum(exp(eta) * R - eta * O)
}
#' Negative lok-likelihood in the age-cohort-interaction model
#'
#' @param par parameters of the age-cohort model in the form \code{c(mu, alpha and beta)}
#' @param O Observed events as returned by \code{\link{exhaustive_stat_2d}}
#' @param R Time at risk as returned by \code{\link{exhaustive_stat_2d}}
#' @return The vector of derivatives of the negative log-likelihood as
#' described in TODO NAME OF REFERENCE SHEET
#' @examples
#' par <- rnorm(ncol(O) + nrow(O))
#' loglik_ac(par, O, R) # length 1 vector
#' score_ac(par, O, R) # vector of length (ncol(O) + nrow(O))
#' hessian_ac(par, O, R) # square matrix of dimension (ncol(O) + nrow(O))
loglik_aci  <- function(par, O, R, pen, weights_age = NULL, weights_cohort = NULL) {
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
    sum(weights_age * t(apply(ext_delta[-1, ], 1, diff)) ^ 2) +
      sum(weights_cohort * apply(ext_delta[, -1], 2, diff) ^ 2)
  )
  sum(exp(eta) * R - eta * O) + pen_term
}
score_aci <- function(par, O, R, pen, weights_age = NULL, weights_cohort = NULL) {
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
    weights_age * t(apply(ext_delta[-1, ], 1, diff)) -
      (weights_age * t(apply(ext_delta[-1, ], 1, diff))) %>% cbind(0) %>% "["(1:(K - 1), -1) +
      weights_cohort * apply(ext_delta[, -1], 2, diff) -
      (weights_cohort * apply(ext_delta[, -1], 2, diff)) %>% rbind(0) %>% "["(-1, 1:(J - 1))
  )
  deriv_delta_mat <- (exp(outer(ext_beta, mu + ext_alpha, FUN = "+") + ext_delta) * R - O)[-1, -1]
  deriv_delta <- as.vector(deriv_delta_mat + pen_term)
  c(deriv_mu, deriv_alpha, deriv_beta, deriv_delta)
}
hessian_aci_no_band <- function(par, O, R, pen, weights_age = NULL, weights_cohort = NULL) {
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
hessian_aci <- function(par, O, R, pen, weights_age = NULL, weights_cohort = NULL,
                            use_band = FALSE) {
  K <- nrow(O)
  J <- ncol(O)
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
                                   index_weights_age, -as.vector(weights_age[, -1])),
                             index_weights_cohort, -as.vector(weights_cohort[-1, ]))
  deriv_diag_delta <- diag(as.vector(
    (exp(outer(ext_beta, mu + ext_alpha, FUN = "+") + ext_delta) * R)[-1, -1]),
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
      as.vector((exp(outer(ext_beta, mu + ext_alpha, FUN = "+") + ext_delta) * R)[-1, -1]),
      -pen * "[<-"(rep(0, (K - 1) * (J - 1)), index_cohort, as.vector(weights_cohort[-1, ])),
      matrix(rep(0, (K - 1) * (J - 1) * (K - 3)), (K - 1) * (J - 1), K - 3),
      -pen * "[<-"(rep(0, (K - 1) * (J - 1)), 1:((K - 2) * (J - 1)), as.vector(weights_age[-1, ]))
    )
    list('A' = cbind(rbind(deriv_diag_mu, deriv_alpha_mu, deriv_beta_mu),
                     rbind(t(deriv_alpha_mu), deriv_diag_alpha, deriv_beta_alpha),
                     rbind(t(deriv_beta_mu), t(deriv_beta_alpha), deriv_diag_beta)),
         'B' = t(cbind(deriv_delta_mu, deriv_delta_alpha, deriv_delta_beta)),
         'D' = deriv_diag_delta)
  }
}
loglik_aci_sel <- function(par, O, R, sel_array) {
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
  eta <- outer(ext_beta, mu + ext_alpha, FUN = "+") + ext_delta
  sum(exp(eta) * R - eta * O)
}
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
  exp(outer(ext_beta, mu + ext_alpha, FUN = "+") + ext_delta)
}
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
  deriv_delta <- temp_array %>%
    apply(., MARGIN = 3, sum)
  c(deriv_mu, deriv_alpha, deriv_beta, deriv_delta)
}
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
ridge_solver_aci <- function(O, R, pen, use_band = FALSE, maxiter = 1000) {
  old_par <- rep(0, ncol(O) * nrow(O))
  for (iter in 1:maxiter) {
    if (use_band) {
      par <- old_par - band_prod_rotated(hessian_aci(old_par, O, R, pen, use_band = TRUE),
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
wridge_solver_aci <- function(O, R, pen, weights_age, weights_cohort,
                              maxiter = 1000, old_par = NULL, use_band = FALSE) {
  if (is.null(old_par)) old_par <- rep(0, ncol(O) * nrow(O))
  for (iter in 1:maxiter) {
    if (use_band) {
      par <- old_par - band_prod_rotated(
        hessian_aci(old_par, O, R, pen, weights_age, weights_cohort, use_band = TRUE),
        score_aci(old_par, O, R, pen, weights_age, weights_cohort))
    } else {
      par <- old_par - Solve(
        hessian_aci(old_par, O, R, pen, weights_age, weights_cohort),
        score_aci(old_par, O, R, pen, weights_age, weights_cohort))
    }
    if (sum(is.na(abs(par - old_par) / abs(old_par)))) break
    if (max(abs(par - old_par) / abs(old_par)) <= 1e-7) break
    old_par <- par
  }
  if (iter == maxiter) {
    warning("Warning: Newton-Raphson procedure did not converge")
  }
  list("par" = par, "convergence" = iter != maxiter, "niter" = iter)
}
aridge_solver_aci <- function(O, R, pen_vect, sample_size,
                              use_band = FALSE,
                              maxiter = 1000 * length(pen_vect)) {
  sel_ls <- par_sel_ls <- haz_sel_ls <- vector('list', length(pen_vect))
  bic <- aic <- ebic <- rep(NA, length(pen_vect))
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
    cat("iter =", iter, '\n')
    old_valve_age <- valve_age
    old_valve_cohort <- valve_cohort
    par <- wridge_solver_aci(O, R, pen = pen_vect[ind_pen],
                             weights_age, weights_cohort,
                             old_par = old_par, use_band = use_band)$par
    delta <- matrix(par[-(1:(J + K - 1))], K - 1, J - 1)
    weights_age[, ] <- 1 / (t(apply(cbind(0, delta), 1, diff)) ^ 2 + epsilon_age ^ 2)
    weights_cohort[, ] <- 1 / (apply(rbind(0, delta), 2, diff) ^ 2 + epsilon_cohort ^ 2)
    valve_age[, ] <- (weights_age * t(apply(cbind(0, delta), 1, diff)) ^ 2)
    valve_cohort[, ] <- (weights_cohort * apply(rbind(0, delta), 2, diff) ^ 2)
    # persp(cuts_cohort, cuts_age, par2haz_aci(par, J, K)[-1, -1],
    #       xlab = 'cohort', ylab = 'age', zlab = 'hazard', phi = 30, expand = 0.75,
    #       theta = -30, shade = 0.25, ltheta = -15, lphi = 120, ticktype = 'detailed')
    # grid.arrange(levelplot(valve_age), levelplot(valve_cohort))
    # max(abs(par - old_par))
    # converged <- sum(is.na(abs(par - old_par) / abs(old_par))) ||
    #   (max(abs(par - old_par) / abs(old_par)) <= 1e-3)
    max(abs(old_valve_age - valve_age),
        abs(old_valve_cohort - valve_cohort)) %>% print()
    converged2 <- max(abs(old_valve_age - valve_age),
                      abs(old_valve_cohort - valve_cohort)) <= 1e-6
    # old_par <- par
    if (converged2) {
      sel_ls[[ind_pen]] <- valve2sel_aci(valve_age, valve_cohort)
      sel_ls[[ind_pen]] %>%
        melt(varnames = c('cohort', 'age')) %>%
        dplyr::as_data_frame() %>%
        mutate(value = as.factor(value)) %>%
        ggplot(., aes(cohort, age, fill = value)) +
        geom_raster() +
        theme(legend.position = 'none')
      L <- max(sel_ls[[ind_pen]])
      sel_array <- lapply(1:L, function(ind) sel_ls[[ind_pen]] == ind) %>%
        unlist() %>%
        array(., dim = c(K - 1, J - 1, L))
      par_sel_ls[[ind_pen]] <- ridge_solver_aci_sel(O, R, sel_array)$par
      haz_sel_ls[[ind_pen]] <- par2haz_aci_sel(par_sel_ls[[ind_pen]], O, R, sel_array) %>%
        'dimnames<-'(dimnames(O))
      bic[ind_pen] <- log(sample_size) * L +
        2 * loglik_aci_sel(par_sel_ls[[ind_pen]], O, R, sel_array)
      ebic[ind_pen] <- bic[ind_pen] +  2 * log(choose(J * K, L) )
      aic[ind_pen] <- 2 * L + 2 * loglik_aci_sel(par_sel_ls[[ind_pen]], O, R, sel_array)
      cat('progress:', ind_pen / length(pen_vect) * 100, '% \n')
      ind_pen <-  ind_pen + 1
      # par <- par * 0
      # weights_age <- weights_age * 0 + 1
      # weights_cohort <- weights_cohort * 0 + 1
    }
    old_par <- par
    if (ind_pen == length(pen_vect) + 1) break
  }
  if (iter == maxiter) {
    warning("Warning: aridge did not converge")
  }
  return(list("par" = par_sel_ls, "haz" = haz_sel_ls, "sel" = sel_ls,
              "bic" = bic, "aic" = aic, 'ebic' = ebic))
}
###
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
loglik_sel_interaction <- function(par, O, R) {
  sum(exp(par) * R - "[<-"(par * O, which(is.nan(par * O), arr.ind = TRUE), 0))
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
###
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
solver_mle_sel <- function(O, R, sel) {
  L <- nlevels(as.factor(sel))
  K <- nrow(O)
  J <- ncol(O)
  par_init <- rep(0, 1 + J - 1 + K - 1 + L)
  par_res <- optim(par_init, fn = loglik_sel, O = O, R = R,
                   sel = sel, method = "BFGS", control = list(maxit = 5000))$par
  haz_res <- par2haz_sel(par_res, sel, J, K, haz.log = FALSE)
  list("par" = par_res, "haz" = haz_res)
}
cv_aridge <- function(pen_vect, nfold, data, cuts_age, cuts_cohort) {
  score_mat <- matrix(NA, nfold, length(pen_vect))
  for (ind in 1:nfold) {
    sample_test <- seq(floor(nrow(data)/nfold) * (ind - 1) + 1,
                       floor(nrow(data)/nfold) * ind)
    sample_train <- setdiff(1:nrow(data), sample_test)
    exhaust_train <- exhaustive_stat(dplyr::slice(data, sample_train), cuts_age, cuts_cohort)
    train <- aridge_solver_interaction(exhaust_train$O, exhaust_train$R, pen_vect = pen_vect, sample_size = length(sample_train))
    train_par_ls <- lapply(train$par, function(par) '[<-'(par, which(is.nan(par)), 0))
    train_sel_ls <- train$sel
    exhaust_test_ls <-  lapply(train_sel_ls, function(sel) exhaustive_stat_sel(
      exhaustive_stat(dplyr::slice(data, sample_test),  cuts_age, cuts_cohort), sel))
    score_mat[ind, ] <- mapply(FUN = function(par, exhaust_test) loglik_sel_interaction(par, exhaust_test$O,
                                                                                        exhaust_test$R),
                               par = train_par_ls, exhaust_test = exhaust_test_ls)
    if (any(is.null(score_mat[ind, ]))) {
      stop('Error in call to aridge_solver_interaction')
    }
  }
  colSums(score_mat)
}
cv_ridge <- function(pen_vect, nfold, data, cuts_age, cuts_cohort) {
  score_mat <- matrix(NA, nfold, length(pen_vect))
  for (ind in 1:nfold) {
    sample_test <- seq(floor(nrow(data)/nfold) * (ind - 1) + 1,
                       floor(nrow(data)/nfold) * ind)
    sample_train <- setdiff(1:nrow(data), sample_test)
    exhaust_train <- exhaustive_stat(dplyr::slice(data, sample_train), cuts_age, cuts_cohort)
    exhaust_test <- exhaustive_stat(dplyr::slice(data, sample_test), cuts_age, cuts_cohort)
    train_par_ls <- lapply(pen_vect, ridge_solver_interaction, O = exhaust_train$O, R = exhaust_train$R) %>%
      lapply(., function(element) element$par)
    score_mat[ind, ] <- lapply(train_par_ls, loglik_interaction, exhaust_test$O,
                               exhaust_test$R, pen = 0) %>% unlist()
    if (any(is.null(score_mat[ind, ]))) {
      stop('Error in call to ridge_solver_interaction')
    }
  }
  colSums(score_mat)
}
###
par2haz_sel <- function(par, sel, J, K, haz.log = FALSE) {
  L <- nlevels(as.factor(sel))
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
  if (haz.log) outer(ext_beta, mu + ext_alpha, FUN = "+") + ext_delta
  else exp(outer(ext_beta, mu + ext_alpha, FUN = "+") + ext_delta)
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
mle_sel_old <- function(O, R, sel) {
  L <- nlevels(as.factor(sel))
  K <- nrow(O)
  J <- ncol(O)
  par_init <- rep(0, 1 + J - 1 + K - 1 + L)
  par_res <- optim(par_init, fn = loglik_aci_sel, O = O, R = R,
                   sel = sel, L = L, method = "BFGS")$par
  haz_res <- par2haz_sel(par_res, sel, J, K, haz.log = FALSE)
  list("par" = par_res, "haz" = haz_res)
}
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
##
aci_cv_aridge <- function(pen_vect, nfold, data, cuts_age, cuts_cohort) {
  score_mat <- matrix(NA, nfold, length(pen_vect))
  for (ind in 1:nfold) {
    sample_test <- seq(floor(nrow(data)/nfold) * (ind - 1) + 1,
                       floor(nrow(data)/nfold) * ind)
    sample_train <- setdiff(1:nrow(data), sample_test)
    exhaust_train <- exhaustive_stat(dplyr::slice(data, sample_train), cuts_age, cuts_cohort)
    train <- aridge_solver_interaction(exhaust_train$O, exhaust_train$R, pen_vect = pen_vect, sample_size = length(sample_train))
    train_par_ls <- lapply(train$par, function(par) '[<-'(par, which(is.nan(par)), 0))
    train_sel_ls <- train$sel
    exhaust_test_ls <-  lapply(train_sel_ls, function(sel) exhaustive_stat_sel(
      exhaustive_stat(dplyr::slice(data, sample_test),  cuts_age, cuts_cohort), sel))
    score_mat[ind, ] <- mapply(FUN = function(par, exhaust_test) loglik_sel_interaction(par, exhaust_test$O,
                                                                                        exhaust_test$R),
                               par = train_par_ls, exhaust_test = exhaust_test_ls)
    if (any(is.null(score_mat[ind, ]))) {
      stop('Error in call to aridge_solver_interaction')
    }
  }
  colSums(score_mat)
}
aci_cv_ridge <- function(pen_vect, nfold, data, cuts_age, cuts_cohort) {
  score_mat <- matrix(NA, nfold, length(pen_vect))
  for (ind in 1:nfold) {
    sample_test <- seq(floor(nrow(data)/nfold) * (ind - 1) + 1,
                       floor(nrow(data)/nfold) * ind)
    sample_train <- setdiff(1:nrow(data), sample_test)
    exhaust_train <- exhaustive_stat(dplyr::slice(data, sample_train), cuts_age, cuts_cohort)
    exhaust_test <- exhaustive_stat(dplyr::slice(data, sample_test), cuts_age, cuts_cohort)
    train_par_ls <- lapply(pen_vect, ridge_solver, O = exhaust_train$O, R = exhaust_train$R) %>%
      lapply(., function(element) element$par)
    score_mat[ind, ] <- lapply(train_par_ls, loglik, exhaust_test$O, exhaust_test$R, pen = 0) %>%
      unlist()
    if (any(is.null(score_mat[ind, ]))) {
      stop('Error in call to ridge_solver_interaction')
    }
  }
  colSums(score_mat)
}
check_derivate_score <- function(point, func, deriv, ...) {
  epsilon_vect <- 10 ^ -seq(1, 6, by = 1)
  p <- length(point)
  l <- length(epsilon_vect)
  score_df <- data.frame(
    score_ = c(rep(deriv(point, ...), l), rep(NA, p * l)),
    # score_ = c(rep(deriv(point, exhaust, adj, pen), l), rep(NA, p * l)),
    type = c(rep("functional", p * l), rep("numerical", p * l)),
    index = rep(rep(seq_along(point), l), 2),
    epsilon = factor(rep(rep(epsilon_vect, each = p), 2), levels = epsilon_vect)
  )
  score_num <- rep(NA, p)
  for (epsilon_ind in epsilon_vect) {
    for (pos in seq_along(score_num)) {
      epsilon_vect_pos <- "[<-"(rep(0, p), pos, epsilon_ind)
      score_num[pos] <- (func(point + epsilon_vect_pos, ...) - func(point, ...)) / epsilon_ind
      # score_num[pos] <- (loglik_graph(point + epsilon_vect_pos, exhaust, adj, pen) -
      #                      loglik_graph(point, exhaust, adj, pen)) / epsilon_ind
    }
    score_df$score_[score_df$epsilon == epsilon_ind & score_df$type == "numerical"] <- score_num
  }
  functional_df <- score_df %>% filter(type == "functional")
  numerical_df <- score_df %>% filter(type == "numerical")
  mse_df <- functional_df %>%
    mutate(squared_diff = (functional_df$score_ - numerical_df$score_) ^ 2) %>%
    group_by(epsilon) %>%
    summarise(mse = sum(squared_diff)) %>%
    ungroup() %>%
    as.data.frame() %>%
    mutate(epsilon = as.numeric(levels(epsilon)))
  list('vect' = score_df, 'mse' = mse_df)
}
check_derivate_hessian <- function(point, func, deriv, ...) {
  epsilon_vect <- 10 ^ -seq(1, 6, by = 1)
  p <- length(point)
  l <- length(epsilon_vect)
  hessian_df <- data.frame(
    hessian_ = c(rep(as.vector(deriv(point, ...)), l),
                 rep(NA, p ^ 2 * l)),
    type = c(rep("functional", p ^ 2 * l),
             rep("numerical", p ^ 2 * l)),
    index = rep(rep(1:(p ^ 2), l), 2),
    epsilon = factor(rep(rep(epsilon_vect, each = p ^ 2), 2), levels = epsilon_vect)
  )
  hessian_num <- rep(NA, p ^ 2)
  for (epsilon_ind in epsilon_vect) {
    for (pos in seq(1, p)) {
      epsilon_vect_pos <- "[<-"(rep(0, p), pos, epsilon_ind)
      hessian_num[(pos - 1) * p + seq(p)] <-
        (func(point + epsilon_vect_pos, ...) - func(point, ...)) / epsilon_ind
    }
    hessian_df$hessian_[hessian_df$epsilon == epsilon_ind & hessian_df$type == "numerical"] <- hessian_num
  }
  functional_df <- hessian_df %>% filter(type == "functional")
  numerical_df <- hessian_df %>% filter(type == "numerical")
  mse_df <- functional_df %>% mutate(squared_diff = (functional_df$hessian_ - numerical_df$hessian_) ^ 2) %>%
    group_by(epsilon) %>%
    summarise(mse = sum(squared_diff)) %>%
    ungroup() %>%
    as.data.frame() %>%
    mutate(epsilon = as.numeric(levels(epsilon)))
  list('vect' = hessian_df, 'mse' = mse_df)
}
