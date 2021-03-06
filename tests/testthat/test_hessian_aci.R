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
K <- 3
J <- 4
O <- matrix(rpois(K * J, 1), K, J)
R <- matrix(rpois(K * J, 10), K, J)
par <- rnorm(K * J)
pen <- 1000
weights_age <- matrix(rnorm((K - 1) * (J - 1)), K - 1, J - 1)
weights_cohort <- matrix(rnorm((K - 1) * (J - 1)), K - 1, J - 1)
epsi <- 1e-9

test_that("hessian_aci is the derivative of score_aci", {
  expect_equal(hessian_aci(par, O, R, pen, weights_age, weights_cohort)[1, ],
               (score_aci(par + epsi * (1:(K * J) == 1), O, R, pen, weights_age, weights_cohort) -
                  score_aci(par, O, R, pen, weights_age, weights_cohort)) / epsi, tol = 1e-3)
  expect_equal(hessian_aci(par, O, R, pen, weights_age, weights_cohort)[2, ],
               (score_aci(par + epsi * (1:(K * J) == 2), O, R, pen, weights_age, weights_cohort) -
                  score_aci(par, O, R, pen, weights_age, weights_cohort)) / epsi, tol = 1e-3)
  expect_equal(hessian_aci(par, O, R, pen, weights_age, weights_cohort)[3, ],
               (score_aci(par + epsi * (1:(K * J) == 3), O, R, pen, weights_age, weights_cohort) -
                  score_aci(par, O, R, pen, weights_age, weights_cohort)) / epsi, tol = 1e-3)
  expect_equal(hessian_aci(par, O, R, pen, weights_age, weights_cohort)[4, ],
               (score_aci(par + epsi * (1:(K * J) == 4), O, R, pen, weights_age, weights_cohort) -
                  score_aci(par, O, R, pen, weights_age, weights_cohort)) / epsi, tol = 1e-3)
  expect_equal(hessian_aci(par, O, R, pen, weights_age, weights_cohort)[5, ],
               (score_aci(par + epsi * (1:(K * J) == 5), O, R, pen, weights_age, weights_cohort) -
                  score_aci(par, O, R, pen, weights_age, weights_cohort)) / epsi, tol = 1e-3)
  expect_equal(hessian_aci(par, O, R, pen, weights_age, weights_cohort)[6, ],
               (score_aci(par + epsi * (1:(K * J) == 6), O, R, pen, weights_age, weights_cohort) -
                  score_aci(par, O, R, pen, weights_age, weights_cohort)) / epsi, tol = 1e-3)
})

set.seed(0)
K <- 3
J <- 4
O <- matrix(rpois(K * J, 1), K, J)
R <- matrix(rpois(K * J, 10), K, J)
par <- rnorm(K * J)
pen <- 1000
weights_age <- matrix(rnorm((K - 1) * (J - 1)), K - 1, J - 1)
weights_cohort <- matrix(rnorm((K - 1) * (J - 1)), K - 1, J - 1)
lim <- J + K - 1
rot <- hessian_aci(par, O, R, pen, weights_age, weights_cohort, use_band = TRUE)
mat <- hessian_aci(par, O, R, pen, weights_age, weights_cohort, use_band = FALSE)

test_that("hessian_aci yiels the same result with and without use_band", {
  expect_equal(rot[[1]], mat[1:lim, 1:lim])
  expect_equal(rot[[2]], mat[1:lim, -(1:lim)])
  expect_equal(rot2mat(rot[[3]]), mat[-(1:lim), -(1:lim)] %>% unname())
})
