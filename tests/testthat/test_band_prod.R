context("band_prod.R")

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
test_that("band_prod is equivalent to solve", {
  expect_equal(band_prod(mat, vect, lim = lim), solve(mat, vect))
  expect_equal(band_prod(mat_rotated, vect, mat_as_rotated = TRUE), solve(mat, vect))
})




