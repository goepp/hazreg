context("solver_graph_aridge.R")

library(igraph)
library(lattice)
library(mgcv)
library(tidyverse)
library(limSolve)

gamma <- c(2, 3, 2.5, -8, -7.5, 7, 7.5)
sigma_sq <- c(10, 20, 10, 10, 2, 10, 10)
temp <- matrix(c(1, 2, 2, 3, 3, 4, 4, 5, 4, 6, 5, 7, 7, 6, 6, 3, 3, 1) %>% as.character(),
               ncol = 2, byrow = TRUE)
graph <- graph_from_edgelist(temp, directed = FALSE)
plot(graph)
adj <- igraph::as_adjacency_matrix(graph, sparse = FALSE)
X <-  diag(length(gamma)) %>%
  "colnames<-"(colnames(adj))

pen <- 1.3
weight_adj <- adj
K <- diag(apply(weight_adj, 1, sum)) - weight_adj
set.seed(0)
par <- rnorm(length(gamma))


test_that("loglik corresponds to the log-likelihood of the penalized regression",  {
  expect_equal(optim(par, loglik, score, gamma = gamma, sigma_sq = sigma_sq,
                     K = K, pen = pen, control = list(maxit = 100000), method = "BFGS")$par %>% unname(),
               gam(gamma ~ X - 1, family = gaussian(),
                   data = data.frame(gamma = gamma, X = X),
                   drop.intercept = TRUE, weights = 1 / sigma_sq, H = 2 * pen * K)$coefficients %>% unname(),
               tolerance = 1e-5)
})

epsi <- 1e-9
test_that("score is the derivative of loglik", {
  expect_equal(score(par, gamma, sigma_sq, K, pen)[1],
               (loglik(par + epsi * (seq_along(par) == 1), gamma, sigma_sq, K, pen) -
                  loglik(par, gamma, sigma_sq, K, pen)) / epsi,
               tolerance = 1e-4)
  expect_equal(score(par, gamma, sigma_sq, K, pen)[2],
               (loglik(par + epsi * (seq_along(par) == 2), gamma, sigma_sq, K, pen) -
                  loglik(par, gamma, sigma_sq, K, pen)) / epsi,
               tolerance = 1e-4)
  expect_equal(score(par, gamma, sigma_sq, K, pen)[3],
               (loglik(par + epsi * (seq_along(par) == 3), gamma, sigma_sq, K, pen) -
                  loglik(par, gamma, sigma_sq, K, pen)) / epsi,
               tolerance = 1e-4)
  expect_equal(score(par, gamma, sigma_sq, K, pen)[4],
               (loglik(par + epsi * (seq_along(par) == 4), gamma, sigma_sq, K, pen) -
                  loglik(par, gamma, sigma_sq, K, pen)) / epsi,
               tolerance = 1e-4)
  expect_equal(score(par, gamma, sigma_sq, K, pen)[5],
               (loglik(par + epsi * (seq_along(par) == 5), gamma, sigma_sq, K, pen) -
                  loglik(par, gamma, sigma_sq, K, pen)) / epsi,
               tolerance = 1e-4)
  expect_equal(score(par, gamma, sigma_sq, K, pen)[6],
               (loglik(par + epsi * (seq_along(par) == 6), gamma, sigma_sq, K, pen) -
                  loglik(par, gamma, sigma_sq, K, pen)) / epsi,
               tolerance = 1e-4)
  expect_equal(score(par, gamma, sigma_sq, K, pen)[7],
               (loglik(par + epsi * (seq_along(par) == 7), gamma, sigma_sq, K, pen) -
                  loglik(par, gamma, sigma_sq, K, pen)) / epsi,
               tolerance = 1e-4)
})

test_that("hessian is the derivative of score", {
  expect_equal(hessian(par, gamma, sigma_sq, K, pen)[1, ] %>% unname(),
               (score(par + epsi * (seq_along(par) == 1), gamma, sigma_sq, K, pen) -
                  score(par, gamma, sigma_sq, K, pen)) / epsi,
               tolerance = 1e-4)
  expect_equal(hessian(par, gamma, sigma_sq, K, pen)[2, ] %>% unname(),
               (score(par + epsi * (seq_along(par) == 2), gamma, sigma_sq, K, pen) -
                  score(par, gamma, sigma_sq, K, pen)) / epsi,
               tolerance = 1e-4)
  expect_equal(hessian(par, gamma, sigma_sq, K, pen)[3, ] %>% unname(),
               (score(par + epsi * (seq_along(par) == 3), gamma, sigma_sq, K, pen) -
                  score(par, gamma, sigma_sq, K, pen)) / epsi,
               tolerance = 1e-4)
  expect_equal(hessian(par, gamma, sigma_sq, K, pen)[4, ] %>% unname(),
               (score(par + epsi * (seq_along(par) == 4), gamma, sigma_sq, K, pen) -
                  score(par, gamma, sigma_sq, K, pen)) / epsi,
               tolerance = 1e-4)
  expect_equal(hessian(par, gamma, sigma_sq, K, pen)[5, ] %>% unname(),
               (score(par + epsi * (seq_along(par) == 5), gamma, sigma_sq, K, pen) -
                  score(par, gamma, sigma_sq, K, pen)) / epsi,
               tolerance = 1e-4)
  expect_equal(hessian(par, gamma, sigma_sq, K, pen)[6, ] %>% unname(),
               (score(par + epsi * (seq_along(par) == 6), gamma, sigma_sq, K, pen) -
                  score(par, gamma, sigma_sq, K, pen)) / epsi,
               tolerance = 1e-4)
  expect_equal(hessian(par, gamma, sigma_sq, K, pen)[7, ] %>% unname(),
               (score(par + epsi * (seq_along(par) == 7), gamma, sigma_sq, K, pen) -
                  score(par, gamma, sigma_sq, K, pen)) / epsi,
               tolerance = 1e-4)
})

test_that("solver_nr fits the penalized regression", {
  expect_equal(solver_nr(gamma, sigma_sq, K, pen) %>% unname(),
               gam(gamma ~ X - 1, family = gaussian(),
                   data = data.frame(gamma = gamma, X = X),
                   drop.intercept = TRUE, weights = 1 / sigma_sq, H = 2 * pen * K)$coefficients %>% unname(),
               tolerance = 1e-6)
})
lambda <- 10 ^ seq(-5, 3, length = 10)

test_that("solvers with gam and with nr match", {
  expect_equal(solver_graph_aridge(gamma, sigma_sq, adj, lambda) %>% unlist(),
               solver_graph_aridge_2(gamma, sigma_sq, adj, lambda) %>% unlist())
})

##
gamma <- c(2, 3, 2.5, -8, -7.5, 7, 7.5)
sigma_sq <- c(10, 20, 10, 10, 2, 10, 10)
temp <- matrix(c(1, 2, 2, 3, 3, 4, 4, 5, 4, 6, 5, 7, 7, 6, 6, 3, 3, 1) %>% as.character(),
               ncol = 2, byrow = TRUE)
graph <- graph_from_edgelist(temp, directed = FALSE)
adj <- igraph::as_adjacency_matrix(graph, sparse = TRUE)
K <- diag(apply(adj, 1, sum)) - adj
pen <- 1.3
set.seed(0)
par <- norm(length(gamma))
mat <- hessian(par, gamma, sigma_sq, K, pen)
vect <- score(par, gamma, sigma_sq, K, pen)

test_that("solve_ssdp corresponds with solve", {
  expect_equal(solve(mat, vect) %>% as.vector(), solve_ssdp(mat, vect))
})
