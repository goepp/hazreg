context("hessian_aci.R")

test_that("O, R, and pen must be positive", {
  expect_error(hessian_aci(1, -1, 1, 0),
               'Error: O, R, and pen must have non-negative values')
  expect_error(hessian_aci(1, 1, -1, 0),
               'Error: O, R, and pen must have non-negative values')
  expect_error(hessian_aci(1, 1, 1, -1),
               'Error: O, R, and pen must have non-negative values')
})

test_that("dimensions of par, O, and R must agree", {
  expect_error(
    hessian_aci(rnorm(100), matrix(1, 10, 10), matrix(1, 10, 15), 0),
    'Error: dimensions of O, R, and pen must agree'
  )
  expect_error(
    hessian_aci(rnorm(100), matrix(1, 10, 10), matrix(1, 10, 15), 0),
    'Error: dimensions of O, R, and pen must agree'
  )
})

set.seed(0)
K <- 2
J <- 3
O <- matrix(rpois(K * J, 1), K, J)
R <- matrix(rpois(K * J, 10), K, J)
par <- rnorm(K * J)
pen <- 1000
weights_age <- matrix(rnorm((K - 1) * (J - 1)), K - 1, J - 1)
weights_cohort <- matrix(rnorm((K - 1) * (J - 1)), K - 1, J - 1)
epsi <- 1e-9

test_that("hessian_aci is the derivative of loglik_aci", {
  expect_equal(hessian_aci(par, O, R, pen, weights_age, weights_cohort)[1],
               (loglik_aci(par + epsi * (1:(K * J) == 1), O, R, pen, weights_age, weights_cohort) -
                  loglik_aci(par, O, R, pen, weights_age, weights_cohort)) / epsi, tol = 1e-5)
  expect_equal(hessian_aci(par, O, R, pen, weights_age, weights_cohort)[2],
               (loglik_aci(par + epsi * (1:(K * J) == 2), O, R, pen, weights_age, weights_cohort) -
                  loglik_aci(par, O, R, pen, weights_age, weights_cohort)) / epsi, tol = 1e-5)
  expect_equal(hessian_aci(par, O, R, pen, weights_age, weights_cohort)[3],
               (loglik_aci(par + epsi * (1:(K * J) == 3), O, R, pen, weights_age, weights_cohort) -
                  loglik_aci(par, O, R, pen, weights_age, weights_cohort)) / epsi, tol = 1e-5)
  expect_equal(hessian_aci(par, O, R, pen, weights_age, weights_cohort)[4],
               (loglik_aci(par + epsi * (1:(K * J) == 4), O, R, pen, weights_age, weights_cohort) -
                  loglik_aci(par, O, R, pen, weights_age, weights_cohort)) / epsi, tol = 1e-5)
  expect_equal(hessian_aci(par, O, R, pen, weights_age, weights_cohort)[5],
               (loglik_aci(par + epsi * (1:(K * J) == 5), O, R, pen, weights_age, weights_cohort) -
                  loglik_aci(par, O, R, pen, weights_age, weights_cohort)) / epsi, tol = 1e-5)
  expect_equal(hessian_aci(par, O, R, pen, weights_age, weights_cohort)[6],
               (loglik_aci(par + epsi * (1:(K * J) == 6), O, R, pen, weights_age, weights_cohort) -
                  loglik_aci(par, O, R, pen, weights_age, weights_cohort)) / epsi, tol = 1e-5)
})
