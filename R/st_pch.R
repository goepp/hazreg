rsurv_graph <- function(sample_size_ls, cuts_ls, alpha_ls, censoring_interval = c(65, 85)) {
  real_time <- mapply(rsurv, sample_size_ls, cuts_ls, alpha_ls)
  censoring <- lapply(sample_size_ls, function(a)
    runif(a, censoring_interval[1], censoring_interval[2]))
  time <- mapply(pmin, real_time, censoring)
  delta <- mapply(function(a, b) a <= b, real_time, censoring)
  list('time' = time, 'delta' = delta)
}
par2df <- function(par, adj, cuts, time_seq = seq(0, 80, 0.1), plot = TRUE) {
  K <- nrow(adj)
  J <- length(par) / K
  alpha <- par %>%
    matrix(., K, J) %>%
    split(., rep(1:K, J))
  hazard <- mapply(stepfun, rep(list(cuts), K), alpha)
  df <- data_frame(time = rep(time_seq, K),
                   hazard = lapply(
                     hazard, function(a) exp(a(time_seq))
                   ) %>% unlist(),
                   region = rep(colnames(adj),
                                each = length(time_seq)))
  if (plot) {
    ggplot(df, aes(time, hazard, color = region)) + geom_step()
  }
  # df
}
exhaustive_stat <- function(time, delta, cuts) {
  k <- length(cuts) + 1
  J <- cut(time, breaks = c(0, cuts, Inf), labels = 1:k, right = FALSE)
  cuts0 <- c(0, cuts)
  O <- sapply(split(delta, J), sum)
  aux <- sapply(split(rep(1, length(delta)), J), sum)[-1]
  R <- sapply(split(time - cuts0[J], J), sum)
  R[-k] <- R[-k] + rev(cumsum(rev(aux))) * diff(cuts0)
  list('O' = O, 'R' = R)
}
exhaustive_stat_graph <- function(sample, cuts) {
  temp <- mapply(exhaustive_stat, time = sample$time, delta = sample$delta, MoreArgs = list(cuts))
  O <- temp[1, ] %>% do.call(rbind, .)
  R <- temp[2, ] %>% do.call(rbind, .)
  list('O' = O, 'R' = R)
}
loglik_graph <- function(par, exhaust, adj, pen = NULL, weights = NULL) {
  K <- nrow(exhaust$O)
  J <- ncol(exhaust$O)
  eta <- matrix(par, K, J)
  if (is.null(weights)) {
    weights_spatial_support <- outer(1:K, 1:K, function(a, b) a > b)
    weights$spatial <- replicate(J, adj * weights_spatial_support)
    weights$temporal <- matrix(1, K, J - 1)
  }
  if (is.null(pen)) {
    pen$spatial <- 0
    pen$temporal <- 0
  }
  # Unpenalized likelihood
  unpenalized <- sum(exp(eta) * exhaust$R - eta * exhaust$O)
  # Term of spatial regularization
  temp1 <- replicate(K, eta) %>% aperm(c(3, 1, 2))
  temp2 <- replicate(K, eta) %>% aperm(c(1, 3, 2))
  spatial_pen <- pen$spatial / 2 *
    sum(weights$spatial * (temp1 - temp2) ^ 2)
  # Term of temporal regularization
  temporal_pen <- pen$temporal / 2 *
    sum(weights$temporal * t(apply(eta, 1, diff)) ^ 2)
  # Return the penalized negative log-likelihood
  unpenalized + spatial_pen + temporal_pen
}
score_graph <- function(par, exhaust, adj, pen = NULL, weights = NULL) {
  K <- nrow(exhaust$O)
  J <- ncol(exhaust$O)
  eta <- matrix(par, K, J)
  if (is.null(weights)) {
    weights_spatial_support <- outer(1:K, 1:K, function(a, b) a > b)
    weights$spatial <- replicate(J, adj * weights_spatial_support)
    weights$temporal <- matrix(1, K, J - 1)
  }
  if (is.null(pen)) {
    pen$spatial <- 0
    pen$temporal <- 0
  }
  # Unpenalized likelihood
  score_unpenalized <- exp(eta) * exhaust$R - exhaust$O
  # Term of spatial regularization
  temp1 <- replicate(K, eta) %>% aperm(c(3, 1, 2))
  temp2 <- replicate(K, eta) %>% aperm(c(1, 3, 2))
  score_spatial_pen <- pen$spatial * (
    apply(weights$spatial * (temp1 - temp2), c(2, 3), sum) -
      apply(weights$spatial * (temp1 - temp2), c(1, 3), sum)
  )
  # Term of temporal regularization
  temp3 <- weights$temporal * t(apply(eta, 1, diff))
  score_temporal_pen <- pen$temporal * (
    cbind(0, temp3) - cbind(temp3, 0)
  )
  # Return the penalized negative log-likelihood
  as.vector(score_unpenalized + score_spatial_pen + score_temporal_pen)
}
hessian_graph <- function(par, exhaust, adj, pen = NULL, weights = NULL) {
  K <- nrow(exhaust$O)
  J <- ncol(exhaust$O)
  eta <- matrix(par, K, J)
  if (is.null(weights)) {
    weights_spatial_support <- outer(1:K, 1:K, function(a, b) a > b)
    weights$spatial <- replicate(J, adj * weights_spatial_support)
    weights$temporal <- matrix(1, K, J - 1)
  }
  if (is.null(pen)) {
    pen$spatial <- 0
    pen$temporal <- 0
  }
  # Unpenalized likelihood
  hessian_unpenalized <- as.vector((exp(eta) * exhaust$R)) %>% diag()
  # Term of spatial regularization
  hessian_spatial_pen1 <- as.vector(
    apply(weights$spatial, c(1, 3), sum) +
      apply(weights$spatial, c(2, 3), sum)
  ) %>% diag()
  # New spatial_pen2
  hessian_spatial_pen2 <- (weights$spatial + aperm(weights$spatial, c(2, 1, 3))) %>%
    lapply(seq(dim(.)[3]), asub, x = ., dims = 3) %>%
    bdiag() %>%
    as.matrix()
  # Old spatial_pen2
  # temp1 <- (weights$spatial + aperm(weights$spatial, c(2, 1, 3))) %>%
  #   aperm(., dim = c(2, 3, 1)) %>%
  #   matrix(., J, K * K) %>%
  #   as.vector()
  # index_small_diag <- cbind(
  #   rep(1:(J * K), K),
  #   rep(1:J, K * K) + J * rep(0:(K - 1), each = K * J)
  # )
  # hessian_spatial_pen2 <- '[<-'(matrix(0, J * K, J * K),
  #                               index_small_diag,
  #                               temp1)
  hessian_spatial_pen <- pen$spatial * (
    hessian_spatial_pen1 - hessian_spatial_pen2
  )
  # Term of temporal regularization
  hessian_temporal_pen1 <- as.vector(
    cbind(0, weights$temporal) + cbind(weights$temporal, 0)
  ) %>% diag()
  index_subdiag <- cbind(
    1:(K * (J - 1)) + K,
    1:(K * (J - 1))
  )
  hessian_temporal_pen2 <- '[<-'(matrix(0, K * J, K * J),
                                 index_subdiag,
                                 as.vector((weights$temporal)))
  hessian_temporal_pen <- pen$temporal * (
    hessian_temporal_pen1 -
      hessian_temporal_pen2 -
      t(hessian_temporal_pen2)
  )
  # Return the penalized negative log-likelihood
  hessian_unpenalized + hessian_spatial_pen + hessian_temporal_pen
}
ridge_solver_graph <- function(exhaust, adj, pen = NULL, maxiter = 1000) {
  K <- nrow(exhaust$O)
  J <- ncol(exhaust$O)
  if (is.null(pen)) {
    pen$spatial <- 0
    pen$temporal <- 0
  }
  old_par <- rep(0, K * J)
  for (iter in 1:maxiter) {
    par <- old_par - Solve(hessian_graph(old_par, exhaust, adj, pen),
                           score_graph(old_par, exhaust, adj, pen))
    max(abs(par - old_par) / abs(old_par))
    # par2df(par, adj, cuts, time_seq = seq(0, 80, 0.1), plot = TRUE)
    if (sum(is.na(abs(par - old_par) / abs(old_par)))) {
      break
    }
    if (max(abs(par - old_par) / abs(old_par)) <= 0.000001) {
      break
    }
    old_par <- par
  }
  if (iter == maxiter) warning("Newton-Raphson procedure did not converge")
  list("par" = par, "convergence" = iter != maxiter, "niter" = iter)
}
wridge_solver_graph <- function(exhaust, adj, pen, weights, maxiter = 1000,
                                old_par = NULL) {
  if (is.null(old_par)) {
    old_par <- rep(0, ncol(exhaust$O) * nrow(exhaust$O))
  }
  for (iter in 1:maxiter) {
    par <- old_par - Solve(hessian_graph(old_par, exhaust, adj, pen, weights),
                           score_graph(old_par, exhaust, adj, pen, weights))
    if (sum(is.na(abs(par - old_par) / abs(old_par)))) {
      break
    }
    max(abs(par - old_par) / abs(old_par))
    if (max(abs(par - old_par) / abs(old_par)) <= 1e-8) {
      break
    }
    old_par <- par
  }
  if (iter == maxiter) {
    warning("Newton-Raphson procedure did not converge")
  }
  list("par" = par, "convergence" = iter != maxiter, "niter" = iter)
}
aridge_solver_graph <- function(exhaust, adj, sample_size, cuts,
                                pen_vect = 10 ^ seq(-8, 8, length.out = 50),
                                maxiter = 1000 * length(pen_vect)) {
  epsilon <- 1e-6
  K <- nrow(exhaust$O)
  J <- ncol(exhaust$O)
  sel_ls <- par_sel_ls <- vector('list', length(pen_vect))
  bic <- aic <- fd_bic <- ebic <- rep(NA, length(pen_vect))
  haz_sel_ls <- vector('list', length(pen_vect) * K)
  haz_df <- data_frame(
    time = numeric(0),
    hazard = numeric(0),
    region = character(0),
    pen = numeric(0)
  )
  dim(haz_sel_ls) <- c(K, length(pen_vect))
  weights <- NULL
  weights$spatial <- replicate(J, adj * outer(1:K, 1:K, function(a, b) a > b))
  weights$temporal <- matrix(1, K, J - 1)
  old_par <- rep(0, J * K)
  ind_pen <- 1
  for (iter in 1:maxiter) {
    cat("iter =", iter, '\n')
    pen <- list('spatial' = 0,
                'temporal' = pen_vect[ind_pen])
    par <- wridge_solver_graph(exhaust, adj, pen = pen, weights = weights,
                               old_par = old_par)$par
    eta <- matrix(par, K, J)
    temp1 <- replicate(K, eta) %>% aperm(c(3, 1, 2))
    temp2 <- replicate(K, eta) %>% aperm(c(1, 3, 2))
    weights$spatial <- replicate(J, adj) / ((temp1 - temp2) ^ 2 + epsilon ^ 2)
    weights$temporal[, ] <- 1 / (t(apply(eta, 1, diff)) ^ 2 + epsilon ^ 2)
    max(abs(par - old_par) / abs(old_par))
    sum((max(abs(par - old_par) / abs(old_par)) <= 1e-7))
    if ((max(abs(par - old_par) / abs(old_par)) <= 1e-7)) {
      valve_spatial <- weights$spatial * (temp1 - temp2) ^ 2
      valve_temporal <- weights$temporal * t(apply(eta, 1, diff)) ^ 2
      sel_ls[[ind_pen]] <- valve2sel_graph(valve_temporal, valve_spatial, exhaust, adj)
      exhaust_sel <- exhaustive_stat_sel(list("O" = exhaust$O, "R" = exhaust$R), sel_ls[[ind_pen]])
      par_sel_ls[[ind_pen]] <- log(exhaust_sel$O / exhaust_sel$R) %>%
        '[<-'(which(log(exhaust_sel$O / exhaust_sel$R) == -Inf), -500) %>%
        '[<-'(which(is.nan(log(exhaust_sel$O / exhaust_sel$R))), 0)
      par2df_sel(par_sel_ls[[ind_pen]], sel_ls[[ind_pen]], J, adj, cuts) %>%
        ggplot(., aes(time, hazard, color = region)) +
        geom_line()
      haz_df <- haz_df %>%
        bind_rows(., par2df_sel(par_sel_ls[[ind_pen]], sel_ls[[ind_pen]], J, adj, cuts) %>%
                    mutate(pen = pen_vect[ind_pen]))
      haz_ls <- par2haz_sel(par_sel_ls[[ind_pen]], sel_ls[[ind_pen]], J, adj, cuts)
      for (ind_k in 1:K) {
        haz_sel_ls[[ind_k, ind_pen]] <- haz_ls[[ind_k]]
      }
      bic[ind_pen] <- log(sample_size) * length(par_sel_ls[[ind_pen]]) +
        2 * loglik_sel_graph(par_sel_ls[[ind_pen]], exhaust_sel$O, exhaust_sel$R)
      ebic[ind_pen] <- bic[ind_pen] +
        2 * log(choose(K * J, length(par_sel_ls[[ind_pen]])))
      aic[ind_pen] <- 2 * length(par_sel_ls[[ind_pen]]) +
        2 * loglik_sel_graph(par_sel_ls[[ind_pen]], exhaust_sel$O, exhaust_sel$R)
      fd_bic[ind_pen] <- log(sample_size / (2 * pi)) * length(par_sel_ls[[ind_pen]]) +
        2 * loglik_sel_graph(par_sel_ls[[ind_pen]], exhaust_sel$O, exhaust_sel$R) -
        sum(par_sel_ls[[ind_pen]] + log(exhaust_sel$R))
      cat('progress:', ind_pen / length(pen_vect) * 100, '% \n')
      ind_pen <-  ind_pen + 1
      ## Le hot strat est désactivé avec les trois lignes suivantes
      par <- rep(0, K * J)
      weights$spatial <- replicate(J, adj * outer(1:K, 1:K, function(a, b) a > b))
      weights$temporal <- matrix(1, K, J - 1)
    }
    old_par <- par
    if (ind_pen == length(pen_vect) + 1) {
      break
    }
  }
  if (iter == maxiter) {
    warning("Warning: aridge did not converge")
  }
  return(list("sel" = sel_ls, "par" = par_sel_ls, "haz" = haz_sel_ls, 'df' = haz_df,
              "bic" = bic, "aic" = aic, "fd_bic" = fd_bic, 'ebic' = ebic))
}
valve2sel_graph <- function(valve_temporal, valve_spatial, exhaust, adj) {
  library(igraph)
  epsilon <- 1e-6
  valve_spatial <- (valve_spatial >= 1 - epsilon) * 1
  valve_temporal <- (valve_temporal >= 1 - epsilon) * 1
  valve_spatial_ls <- (replicate(J, adj) - valve_spatial) %>%
    lapply(seq(dim(.)[3]), asub, x = ., dims = 3)
  for (ind_j in 1:J) {
    dimnames(valve_spatial_ls[[ind_j]]) <- lapply(dimnames(adj), function(x) paste0(x, ind_j))
  }
  adj_spatial <- bdiag(valve_spatial_ls) %>% as.matrix()
  adj_temporal <- '[<-'(matrix(0, K * J, K * J),
                        cbind((K + 1):(K * J), 1:(K * J - K)),
                        as.vector(t(1 - valve_temporal)))
  adjacency <- adj_spatial + t(adj_spatial) +
    adj_temporal + t(adj_temporal)
  dimnames(adjacency) <- dimnames(adj_temporal)
  graph <- igraph::graph_from_adjacency_matrix(adjacency, mode = 'undirected')
  temp1 <- c(1, 1.5, 2, 2.75)
  temp2 <- c(2, 1.5, 2, 1.75)
  layout <- cbind(rep(temp1, J), rep(temp2, J) + rep(1:J, each = K))
  plot(graph, layout = layout, edge.color = 'black',
       vertex.label = NA, vertex.color = 'black',
       vertex.size = 1)
  selection <- matrix(clusters(graph)$membership, K, J) %>%
    'dimnames<-'(dimnames(valve_temporal))
  selection
}
exhaustive_stat_sel <- function(exhaust, sel) {
  if (any(dim(exhaust$O) != dim(sel))) {
    stop("error: dimensions of exhaust and selection must agree")
  }
  O <- lapply(levels(as.factor(sel)),
              function(sel_val) sum(exhaust$O[which(sel == sel_val, arr.ind = TRUE)])) %>%
    unlist()
  R <- lapply(levels(as.factor(sel)),
              function(sel_val) sum(exhaust$R[which(sel == sel_val, arr.ind = TRUE)])) %>%
    unlist()
  return(list("O" = O, "R" = R))
}
loglik_sel_graph <- function(par, O, R) {
  sum(exp(par) * R - "[<-"(par * O, which(is.nan(par * O), arr.ind = TRUE), 0))
}
par2df_sel <- function(par, sel, J, adj, cuts) {
  K <- nrow(adj)
  temp <- matrix(NA, K, J)
  for (level_ind in 1:nlevels(as.factor(sel))) {
    temp[sel == level_ind] <- par[level_ind]
  }
  alpha <- temp %>%
    matrix(., K, J) %>%
    split(., rep(1:K, J))
  hazard <- mapply(stepfun, rep(list(cuts), K), alpha)
  df <- data_frame(time = rep(time_seq, K),
                   hazard = lapply(
                     hazard, function(a) exp(a(time_seq))
                   ) %>% unlist(),
                   region = rep(colnames(adj),
                                each = length(time_seq)))
  df
}
par2haz_sel <- function(par, sel, J, adj, cuts) {
  K <- nrow(adj)
  temp <- matrix(NA, K, J)
  for (level_ind in 1:nlevels(as.factor(sel))) {
    temp[sel == level_ind] <- par[level_ind]
  }
  alpha <- temp %>%
    matrix(., K, J) %>%
    split(., rep(1:K, J))
  hazard <- mapply(stepfun, rep(list(cuts), K), alpha)
  hazard
}
