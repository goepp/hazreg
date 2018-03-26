context("block_bandsolve.R")

library(bandsolve)
library(limSolve)

set.seed(0)
A <- matrix(rnorm(36), 6, 6)
B <- matrix(rnorm(24), 6, 4)
D <- cbind(rnorm(4), c(rnorm(3), 0))
mat_rotated <- list('A' = A + t(A),
                    'B' = B,
                    'D' = D)
mat <- rbind(cbind(mat_rotated$A, mat_rotated$B),
             cbind(t(mat_rotated$B), rot2mat(mat_rotated$D)))
vect <- rnorm(10)
lim <- 6

test_that("block_bandsolve is equivalent to solve", {
  expect_equal(block_bandsolve(mat, vect, FALSE, lim = lim), solve(mat, vect))
  expect_equal(block_bandsolve(mat_rotated, vect, mat_as_rotated = TRUE), solve(mat, vect))
})

set.seed(0)
A <- matrix(rnorm(36), 6, 6)
A[, 1] <- 0
A[1, ] <- 0
B <- matrix(rnorm(24), 6, 4)
D <- diag(rnorm(4)) + rbind(cbind(0, diag(rnorm(3))), 0)
mat <- rbind(cbind(A + t(A), B),
             cbind(t(B), D + t(D)))
vect <- rnorm(10)
lim <- 6

test_that("block bandsolve works with A has zero-valued columns and rows", {
  expect_equal(block_bandsolve(mat, vect, FALSE, lim), Solve(mat, vect))
})



