context("aridge_solver_aci.R")

library(lattice)
library(tidyverse)
library(reshape2)
library(limSolve)
library(igraph)
library(bandsolve)
library(progress)

set.seed(0)
K <- 3
J <- 4
O <- matrix(rpois(K * J, 1), K, J)
R <- matrix(rpois(K * J, 10), K, J)
par <- rnorm(K * J)


test_that("aridge_solver_aci yields same result with and without use_band", {
  pen_vect <- 1
  mat <- aridge_solver_aci(O, R, pen_vect, sample_size = 10, use_band = FALSE)
  rot <- aridge_solver_aci(O, R, pen_vect, sample_size = 10, use_band = TRUE)
  expect_equal(mat, rot)
  skip_on_cran()
  pen_vect <- 10 ^ (-2:3)
  mat <- aridge_solver_aci(O, R, pen_vect, sample_size = 10, use_band = FALSE)
  rot <- aridge_solver_aci(O, R, pen_vect, sample_size = 10, use_band = TRUE)
  expect_equal(mat, rot)
})




