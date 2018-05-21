#' Generate simulatied data
#'
#' @family simu_wrapper
#' @param sample_size_vect vector of sample sizes
#' @param n_rep number of times the data are simulated for each sample size
#' @param censoring_lim interval of the censoring. Only cohort-indenpendant censoring is
#' available.
#' @return sample_df simulated events
#' @return exhaust_ls list of exhaustive statistics computed for each simulation
#' @export
design_simu <- function(sample_size_vect, n_rep, censoring_lim = c(75, 100)) {
  design <- expand.grid('sample_size' = sample_size_vect, 'rep' = 1:n_rep)
  sample_df <- data_frame()
  exhaust_ls <- vector('list', nrow(design))
  set.seed(0)
  for (ind in 1:nrow(design)) {
    sample_size <- design$sample_size[ind]; ind_rep <- design$rep[ind]
    censoring <- runif(sample_size, censoring_lim[1], censoring_lim[2])
    surv_data <- map_surv(ground, sample_size = sample_size,
                          dob_dist = "unif", dob_dist_par = c(1900, 2000)) %>%
      mutate(delta = as.numeric(censoring >= age)) %>%
      mutate(age = pmin(age, censoring)) %>%
      mutate(sample_size = sample_size, rep = ind_rep)
    sample_df <- bind_rows(sample_df, surv_data)
    exhaust <- exhaustive_stat_2d(surv_data, cuts_age, cuts_cohort)
    exhaust_ls[[ind]] <- exhaust
  }
  list("sample_df" = sample_df, "exhaust_ls" = exhaust_ls)
}
#' Compare different selection criteria using simulated data
#'
#' @family simu_wrapper
#' @param ind_wrapper index of the simulation design. Between 1 and \code{nrow(design)}.
#' @param design simulation design matrix
#' @param exhaust_ls list of exhaustive statistics computed from simulated data
#' @param pen_vect the vector of penalties
#' @param scaled_ground_haz the true hazard used to simulate data
#' @return The estimated models from each criterion. The results are given in
#' the form of data frames. \code{nlevel} gives the number of selected areas,
#' \code{mse} gives the mean square error, \code{sel} gives the estimated regions,
#' and \©ode{haz} gives the estimated hazard.
#' @export
criterion_simu_wrapper <- function(ind_wrapper, design, exhaust_ls, pen_vect, scaled_ground_haz) {
  sample_size <- design$sample_size[ind_wrapper]; ind_rep <- design$rep[ind_wrapper]
  exhaust <- exhaust_ls[[ind_wrapper]]
  aridge_est <- aridge_solver_interaction(exhaust$O, exhaust$R, pen_vect, sample_size)
  aic_min <- which.min(unlist(aridge_est$aic))
  bic_min <- which.min(unlist(aridge_est$bic))
  ebic_min <- which.min(unlist(aridge_est$ebic))
  # nlevel
  aic_nlevel <- length(aridge_est$par[[aic_min]])
  bic_nlevel <- length(aridge_est$par[[bic_min]])
  ebic_nlevel <- length(aridge_est$par[[ebic_min]])
  nlevel_df <- data_frame('ind_rep' = ind_rep, 'sample_size' = sample_size, variable = 'nlevel',
                          'aic' = aic_nlevel, 'bic' = bic_nlevel, 'ebic' = ebic_nlevel)
  # mse
  aic_mse <- sum((aridge_est$haz[[aic_min]] - scaled_ground_haz) ^ 2) / (J * K)
  bic_mse <- sum((aridge_est$haz[[bic_min]] - scaled_ground_haz) ^ 2) / (J * K)
  ebic_mse <- sum((aridge_est$haz[[ebic_min]] - scaled_ground_haz) ^ 2) / (J * K)
  mse_df <- data_frame('ind_rep' = ind_rep, 'sample_size' = sample_size, variable = 'mse',
                       'aic' = aic_mse, 'bic' = bic_mse, 'ebic' = ebic_mse)
  # sel
  sel_df <- lapply(aridge_est$sel[c(aic_min, bic_min, ebic_min)],
                   melt, varnames = c('cohort', 'age')) %>%
    bind_rows() %>%
    dplyr::as_data_frame() %>%
    mutate(criterion = rep(c('aic', 'bic', 'ebic'), each = K * J)) %>%
    mutate(ind_rep = ind_rep, sample_size = sample_size, variable = 'sel')
  # haz
  haz_df <- lapply(aridge_est$haz[c(aic_min, bic_min, ebic_min)],
                   melt, varnames = c('cohort', 'age')) %>%
    bind_rows() %>%
    dplyr::as_data_frame() %>%
    mutate(criterion = rep(c('aic', 'bic', 'ebic'), each = K * J)) %>%
    mutate(ind_rep = ind_rep, sample_size = sample_size, variable = 'haz')

  list(scalar = bind_rows(nlevel_df, mse_df), sel_haz = bind_rows(sel_df, haz_df))
}
#' Run the interaction model with cross-validation on simulated data
#'
#' @family simu_wrapper
#' @param nfold constant \code{K} in K-fold cross-validation
#' @export
cv_simu_wrapper <- function(ind_wrapper, design, sample_df, pen_vect, scaled_ground_haz, nfold) {
  ind_sample_size <- design$sample_size[ind_wrapper]
  ind_rep <- design$rep[ind_wrapper]
  surv_data <- sample_df %>%
    filter(sample_size == ind_sample_size, rep == ind_rep)
  exhaust <- exhaustive_stat_2d(surv_data, cuts_age, cuts_cohort)
  score_aridge <- cv_aridge(pen_vect, nfold = nfold, surv_data, cuts_age, cuts_cohort)
  aridge_est <- aridge_solver_interaction(exhaust$O, exhaust$R,
                                          pen_vect[which.min(score_aridge)],
                                          ind_sample_size,
                                          use_band = TRUE,
                                          progress = FALSE)
  nlevel_df <- dplyr::data_frame('ind_rep' = ind_rep, 'sample_size' = ind_sample_size,
                                 variable = 'nlevel', 'cv' = aridge_est$par[[1]] %>% length())
  mse_df <- dplyr::data_frame('ind_rep' = ind_rep, 'sample_size' = ind_sample_size, 'variable' = 'mse',
                              'cv' = sum((aridge_est$haz[[1]] - scaled_ground_haz) ^ 2) / (J * K))
  sel_df <- melt(aridge_est$sel[[1]], varnames = c('cohort', 'age')) %>%
    bind_rows() %>%
    dplyr::as_data_frame() %>%
    mutate(criterion = rep('cv', each = K * J)) %>%
    mutate('ind_rep' = ind_rep, 'sample_size' = ind_sample_size, 'variable' = 'sel')
  haz_df <- melt(aridge_est$haz[[1]], varnames = c('cohort', 'age')) %>%
    bind_rows() %>%
    dplyr::as_data_frame() %>%
    mutate(criterion = rep('cv', each = K * J)) %>%
    mutate('ind_rep' = ind_rep, 'sample_size' = ind_sample_size, 'variable' = 'haz')
  return(list(scalar = bind_rows(nlevel_df, mse_df),
              sel_haz = bind_rows(sel_df, haz_df)))
}
#' Run the interaction model with ridge regularization on simulated data
#'
#' @family simu_wrapper
#' @export
ridge_simu_wrapper <- function(ind_wrapper, design, sample_df, pen_vect, scaled_ground_haz) {
  ind_sample_size <- design$sample_size[ind_wrapper]
  ind_rep <- design$rep[ind_wrapper]
  surv_data <- sample_df %>%
    filter(sample_size == ind_sample_size, rep == ind_rep)
  exhaust <- exhaustive_stat_2d(surv_data, cuts_age, cuts_cohort)
  score_ridge <- cv_ridge(nfold = 10, pen_vect, surv_data, cuts_age, cuts_cohort)
  ridge_est <- ridge_solver_interaction(exhaust$O, exhaust$R,
                                        pen_vect[which.min(score_ridge)])
  mse_df <- dplyr::data_frame(
    'ind_rep' = ind_rep, 'sample_size' = ind_sample_size, variable = 'mse',
    'cv_ridge' = sum((ridge_est$par %>% par2haz_interaction(J, K) - scaled_ground_haz) ^ 2) / (J * K))
  haz_df <- ridge_est$par %>%
    par2grid_interaction(cuts_age, cuts_cohort) %>%
    as_tibble() %>%
    mutate(criterion = rep('cv_ridge', each = K * J), 'ind_rep' = ind_rep, 'sample_size' = ind_sample_size, variable = 'haz')
  return(list(mse = mse_df, haz = haz_df))
}
