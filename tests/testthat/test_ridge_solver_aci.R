context("ridge_solver_aci.R")

library(lattice)
library(tidyverse)
library(limSolve)

set.seed(0)
K <- 3
J <- 4
O <- matrix(rpois(K * J, 1), K, J)
R <- matrix(rpois(K * J, 10), K, J)
par <- rnorm(K * J)
pen <- 1000
weights_age <- matrix(rnorm((K - 1) * (J - 1)), K - 1, J - 1)
weights_cohort <- matrix(rnorm((K - 1) * (J - 1)), K - 1, J - 1)

par_mat <- ridge_solver_aci(O, R, pen, weights_age, weights_cohort, use_band = FALSE)$par
par_rot <- ridge_solver_aci(O, R, pen, weights_age, weights_cohort, use_band = TRUE)$par

test_that("ridge_solver_aci yields the same estimate with use_band = FALSE or TRUE", {
  expect_equal(par_mat %>% unname(), par_rot)
})
