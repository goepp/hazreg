#' hazreg: A package for computing regularized hazard in survival analysis models
#'
#' The foo package provides functions for the estimation of hazard
#' of different models in time-to-event data.
#'
#' @section Utility functions
#' @section The Age-cohort model
#' @section The Age-cohort interaction model
#' @section The interaction model
#'
#' @docType package
#' @name hazreg
NULL
#' Generate geometric sequences
#'
#' @param x A strictly positive number
#' @param y A strictly positive number
#' @param length.out The length of the sequence. It must be a strictly positive integer
#' @return A geometric sequence
#' @examples
#' seq_geom(1, 10, 10)
#' seq_geom(1, 10, 10)
#' @export
seq_geom <- function(a, b, length.out) {
  return(a * (b / a) ^ ((0:(length.out - 1))/(length.out - 1)))
}
#' @export
reverselog_trans <- function(base = exp(1)) {
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
#' @param age_seq Age discretizatioin vector
#' @param cohort_seq Cohort discretizatioin vector
#' @param islet_ls List of parameters of the interaction terms.
#' @return
#' The interactions terms are densities of two dimensional normal distributions.
#' \code{islet_ls} gives the mean, variance-covariance matrix and multiplicative coefficients of each
#' interaction term.
#' @examples
#' \dontrun{
#' library(lattice)
#' mu <- log(1e-2)
#' age_seq <- seq(0, 100, length.out = 100)
#' cohort_seq <- seq(1900, 1900 + age_seq[101], length.out = age_length + 1)
#' age_coef <- seq(0, 2, length.out = age_length - 1)
#' cohort_coef <- seq(0, 0.8, length.out = age_length - 1)
#' islet_ls = list(list(mean = c(1945, 45), sigma = 100 * diag(2), coef = 40))
#' hazard <- map_maker(mu, age_coef, cohort_coef, age_seq, cohort_seq, islet_ls = islet_ls)
#' persp(hazard$sheet)
#' }
#' @export
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
#' \dontrun{
#' map_maker(100, c(20, 40, 60), c(0.01, 0.03, 0.05, 0.08))
#' }
#' @export
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
#' \dontrun{
#' library(tidyverse)
#' mu <- log(1e-2)
#' age_seq <- seq(0, 100, length.out = 100)
#' cohort_seq <- seq(1900, 1900 + age_seq[101], length.out = age_length + 1)
#' age_coef <- seq(0, 2, length.out = age_length - 1)
#' cohort_coef <- seq(0, 0.8, length.out = age_length - 1)
#' islet_ls = list(list(mean = c(1945, 45), sigma = 100 * diag(2), coef = 40))
#' hazard <- map_maker(mu, age_coef, cohort_coef, age_seq, cohort_seq, islet_ls = islet_ls)
#' (map_surv(hazard, 10))
#' }
#' @export
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
exhaustive_stat <- function(surv_data, cuts_age, cuts_cohort) {
  ext_cuts_age <- c(0, cuts_age)
  ext_cuts_cohort <- c(min(surv_data$cohort), cuts_cohort)
  temp1 <- model.matrix(~0 + cut(surv_data$age, breaks = c(ext_cuts_age, Inf)))
  temp2 <- t(1 - apply(temp1, 1, cummax)) + temp1
  temp3 <- t(apply(apply(
    cbind(t(matrix(rep(ext_cuts_age, length(surv_data$age)), ncol = length(surv_data$age))),
          surv_data$age), 1, sort),
    2, diff))
  R_old <- as.data.frame(unname(temp2 * temp3))
  rearrange_by_cohort <- function(vect) split(vect, cut(
    surv_data$cohort, breaks = c(ext_cuts_cohort, Inf),right = FALSE, dig.lab = 4))
  sum_by_cohort <- function(vect) unlist(sapply(rearrange_by_cohort(vect), sum))
  R <- sapply(R_old, sum_by_cohort)
  O_old <- as.data.frame(unname(temp1 * matrix(rep(surv_data$delta,
                                                   length(ext_cuts_age)), ncol = length(ext_cuts_age))))
  O <- sapply(O_old, sum_by_cohort)
  colnames(O) <- colnames(R) <- levels(cut(0, breaks = c(ext_cuts_age, Inf), right = FALSE, dig.lab = 3))
  list("O" = O, "R" = R)
}
#' @export
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
#' @param mat Square matrix
#' @param vect Vector of same length as \code{mat}
#' @param lim Positive integer
#' @return The vector \code{mat} ^ {-1} \code{vect} is returned.
#' If \code{mat[lim:end, lim:end]} is not a band matrix or if \code{mat} is
#' symmetric, an error is returned.
#' If these two conditions are met, the returned value is the same as
#' \code{solve(mat, vect)}. There might be nummerical impercision errors.
#' @examples
#' \dontrun{
#' ## With mat_as_rotated = FALSE
#' set.seed(0)
#' A <- matrix(rnorm(36), 6, 6)
#' B <- matrix(rnorm(24), 6, 4)
#' D <- diag(rnorm(4)) + rbind(cbind(0, diag(rnorm(3))), 0)
#' mat <- rbind(cbind(A + t(A), B),
#'                      cbind(t(B), D + t(D)))
#' vect <- rnorm(10)
#' lim <- 6
#' max(block_bandsolve(mat, vect, lim) - solve(mat, vect))
#'
#' ## With mat_as_rotated = TRUE
#' set.seed(0)
#' A <- matrix(rnorm(36), 6, 6)
#' B <- matrix(rnorm(24), 6, 4)
#' D <- cbind(rnorm(4), c(rnorm(3), 0))
#' mat_rotated <- list('A' = A + t(A),
#'             'B' = B,
#'             'D' = D)
#' mat <- rbind(cbind(mat_rotated$A, mat_rotated$B),
#'                      cbind(t(mat_rotated$B), rot2mat(mat_rotated$D)))
#' vect <- rnorm(10)
#' max(block_bandsolve(mat_rotated, vect, mat_as_rotated = TRUE) - solve(mat, vect))
#' }
#' @export
block_bandsolve <- function(mat, vect, mat_as_rotated = FALSE, lim = NULL) {
  if (mat_as_rotated) {
    if (nrow(mat$A) != ncol(mat$A)) {
      stop("Error: upper-left matrix not square")
    }
    # This code is sensible to extreme values of mat$A
    # if (max(mat$A - t(mat$A)) >= 1e-10) {
    #   stop("Error: upper-left matrix not symmetric")
    # }
    if ((nrow(mat$A) + nrow(mat$D)) != length(vect)) {
      stop('Error: dimensions of matrix and vector must agree')
    }
    vect1 <- vect[1:ncol(mat$A)]
    vect2 <- vect[-seq_along(vect1)]
    A <- mat$A
    B <- mat$B
    rotD <- mat$D
  } else {
    if (is.null(lim)) {
      stop('Error: lim must be specified when mat_as_rotated is FALSE')
    }
    if (!isSymmetric(mat) || ncol(mat) != length(as.vector(vect))) {
      error("Error: dimensions of hessian and score do not match")
    }
    vect1 <- vect[1:lim]
    vect2 <- vect[(lim + 1):(ncol(mat))]
    A <- mat[1:lim, 1:lim]
    B <- mat[1:lim, (lim + 1):ncol(mat)]
    D <- mat[(lim + 1):nrow(mat), (lim + 1):ncol(mat)]
    # rotD <- rbind(
    #   mat2rot(D[, apply(D, 1, function(x) any(x != 0))]),
    #   rep(0, sum(apply(D, 1, function(x) all(x == 0))))
    #   )
    rotD <- mat2rot(D)
  }
  schur <- A - B %*% bandsolve(rotD, t(B))
  temp1 <- Solve(schur, vect1)
  temp2 <- Solve(schur, B %*% bandsolve(rotD, vect2))
  c(as.vector(temp1 - temp2),
    as.vector(-bandsolve(rotD, t(B) %*% temp1) +
                bandsolve(rotD, vect2) +
                bandsolve(rotD, t(B) %*% temp2))
  )
}
#' @export
valve2sel <- function(valve_age, valve_cohort, epsilon = 1e-8) {
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
#' @export
grid2raster <- function(grid_df, title = NULL) {
  ggplot(grid_df, aes(cohort, age, fill = value)) + geom_raster() +
    ggtitle(title) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
}
#' @export
sel2segment <- function(sel, cuts_age, cuts_cohort) {
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
#' @export
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
#' @export
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
#' @export
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
#' @export
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
#' @export
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
#' @export
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
#' @export
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
#' @export
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
#' @export
check_derivate_score <- function(par, func, deriv, ...) {
  epsilon_vect <- 10 ^ -seq(1, 6, by = 1)
  p <- length(par)
  l <- length(epsilon_vect)
  score_df <- data.frame(
    score_ = c(rep(deriv(par, ...), l), rep(NA, p * l)),
    # score_ = c(rep(deriv(par, exhaust, adj, pen), l), rep(NA, p * l)),
    type = c(rep("functional", p * l), rep("numerical", p * l)),
    index = rep(rep(seq_along(par), l), 2),
    epsilon = factor(rep(rep(epsilon_vect, each = p), 2), levels = epsilon_vect)
  )
  score_num <- rep(NA, p)
  for (epsilon_ind in epsilon_vect) {
    for (pos in seq_along(score_num)) {
      epsilon_vect_pos <- "[<-"(rep(0, p), pos, epsilon_ind)
      score_num[pos] <- (func(par + epsilon_vect_pos, ...) - func(par, ...)) / epsilon_ind
      # score_num[pos] <- (loglik_graph(par + epsilon_vect_pos, exhaust, adj, pen) -
      #                      loglik_graph(par, exhaust, adj, pen)) / epsilon_ind
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
#' @export
check_derivate_hessian <- function(par, func, deriv, ...) {
  epsilon_vect <- 10 ^ -seq(1, 6, by = 1)
  p <- length(par)
  l <- length(epsilon_vect)
  hessian_df <- data.frame(
    hessian_ = c(rep(as.vector(deriv(par, ...)), l),
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
        (func(par + epsilon_vect_pos, ...) - func(par, ...)) / epsilon_ind
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
