context("exhaustive_stat_2d")

library(tidyverse)
library(pryr)

load_all()

age_length <- 10
age_seq <- seq(0, 100, length.out = age_length + 1)
cohort_seq <- seq(1900, 1900 + age_seq[age_length + 1], length.out = age_length + 1)
mu <- log(1e-2)
age_coef <- seq(0, 2.5, length.out = age_length - 1)
cohort_coef <- seq(0, 0.30, length.out = age_length - 1)
ground <- map_maker(mu, age_coef, cohort_coef, age_seq, cohort_seq)
cuts_age_ground <- age_seq[-c(1, length(age_seq))]
cuts_cohort_ground <- cohort_seq[-c(1, length(cohort_seq))]

cuts_age <- seq(10, 90, 10)
cuts_cohort <- seq(1910, 1990, 10)
J <- length(cuts_age) + 1; K <- length(cuts_cohort) + 1
sample_size <- 10

censoring <- runif(sample_size, 75, 100)
surv_data <- map_surv(ground, sample_size = sample_size,
                      dob_dist = "unif", dob_dist_par = c(1900, 2000)) %>%
  mutate(delta = as.numeric(censoring >= age)) %>%
  mutate(age = pmin(age, censoring))

ggplot(surv_data, aes(cohort, age)) +
  geom_point() +
  geom_vline(xintercept = cuts_cohort, color = "green") +
  geom_hline(yintercept = cuts_age, color = "green") +
  ylim(0, 100) +
  xlim(1900, 2000)

test_that("exhaustive_stat work with empty cohort intervals", {
  expect_error(exhaustive_stat_2d(surv_data, cuts_age, cuts_cohort), NA)
})

exhaust <- exhaustive_stat_2d(surv_data, cuts_age, cuts_cohort)

