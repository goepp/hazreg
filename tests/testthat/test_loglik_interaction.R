context("loglik_interaction.R")

test_that("O, R, and pen must be positive", {
  expect_error(loglik_interaction(1, -1, 1, 0),
               'Error: O, R, and pen must have non-negative values')
  expect_error(loglik_interaction(1, 1, -1, 0),
               'Error: O, R, and pen must have non-negative values')
  expect_error(loglik_interaction(1, 1, 1, -1),
               'Error: O, R, and pen must have non-negative values')
})

test_that("dimensions of par, O, and R must agree", {
  expect_error(
    loglik_interaction(rnorm(100), matrix(1, 10, 10), matrix(1, 10, 15), 0),
    'Error: dimensions of O, R, and pen must agree')
  expect_error(
    loglik_interaction(rnorm(100), matrix(1, 10, 10), matrix(1, 10, 15), 0),
    'Error: dimensions of O, R, and pen must agree')
})
