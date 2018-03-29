context("test_function_aci_sel")

library(tidyverse)

set.seed(0)
K <- 3
J <- 4
O <- matrix(rpois(K * J, 1), K, J)
R <- matrix(rpois(K * J, 10), K, J)
par <- rnorm(K * J)
sel <- matrix(1:((K - 1) * (J - 1)), K - 1, J - 1)
L <- max(sel)
sel_array <- lapply(1:L, function(ind) sel == ind) %>%
  unlist() %>%
  array(., dim = c(K - 1, J - 1, L))

test_that("loglik_aci_sel consistent with loglik_aci", {
  expect_equal(loglik_aci(par, O, R, pen = 0),
               loglik_aci_sel(par, O, R, sel_array))
})

test_that("score_aci_sel consistent with score_aci", {
  expect_equal(score_aci(par, O, R, pen = 0),
               score_aci_sel(par, O, R, sel_array))
})

test_that("hessian_aci_sel consistent with hessian_aci", {
  expect_equal(hessian_aci(par, O, R, pen = 0),
               hessian_aci_sel(par, O, R, sel_array))
})

