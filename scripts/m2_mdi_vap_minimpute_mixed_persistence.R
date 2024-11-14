#' ###########################################
#' MIXED EFFECT VAP ~ MDI LONGITUDINAL models
#' ###########################################


#' load libraries and set seed / control variables
library(tidyverse)
library(tidybayes)
library(brms)
set.seed(16)
run_models = TRUE
save_models = TRUE
run_plots = TRUE




#' #######################################################################
#' #######################################################################
#' VAP as a function of lower respiratory tract diversity
#' #######################################################################
#' #######################################################################


#' #######################################################################
#' SUSPECTED VAP
#' #######################################################################

if(run_models) {
dat_i %>%
  filter(suspected_vap == TRUE) %>%
  bind_rows(dat) %>%
  arrange(subject_id, date) %>%
  # filter to include dates only up to VAP (or last date, if no VAP)
  group_by(subject_id) %>%
  mutate(first_suspected_vap = ifelse(is.na(first_suspected_vap), max(date), first_suspected_vap)) %>%
  filter(date <= first_suspected_vap) %>%
  ungroup() %>%
  # # 48-hour window for case period
  # group_by(subject_id) %>%
  # arrange(date) %>%
  # mutate(suspected_vap = date == first_suspected_vap | date == first_suspected_vap - 1) %>%
  # ungroup() %>%
  # #
  select(subject_id, subject_day, shannon, suspected_vap) %>%
  filter(complete.cases(.)) %>%
  #summarise(summary(shannon))
  mutate(shannon_less_than_2 = shannon < 2) %>%
  group_by(subject_id) %>%
  arrange(subject_day) %>%
  mutate(shannon_improved = shannon_less_than_2 == FALSE & lag(shannon_less_than_2, 1) == TRUE) %>%
  replace_na(replace = list("shannon_improved" = FALSE)) %>%
  mutate(shannon_flip_minimpute_count = as.character(cumsum(shannon_improved)),
  ) %>%
  tidyr::fill(shannon_flip_minimpute_count, .direction = "down") %>% #arrange(desc(shannon_flip_minimpute_count)) %>%
  group_by(subject_id, shannon_flip_minimpute_count) %>%
  arrange(subject_day) %>%
  #mutate(period_day = subject_day - min(subject_day, na.rm = TRUE)) %>%
  mutate(shannon_low_count = cumsum(shannon_less_than_2),
         #shannon_low_days = ifelse(shannon_low_count > 0, max(period_day, na.rm = TRUE) + 1, 0)
  ) %>%
  ungroup() %>% #qplot(data = ., x = shannon_low_count)
  #qplot(data = ., x = shannon_low_count, shannon_low_days)
  distinct() %>%
  filter(complete.cases(.)) %>%
  identity() -> d_vap

d_vap


#' run mixed effects model with data
d_vap %>%
  brm(suspected_vap ~ shannon_low_count + 0 + (1|subject_id),
      data = .,
      family = "bernoulli",
      iter = 2000,
      warmup = 1000,
      chains = 4,
      cores = 4,
      control = list("adapt_delta" = 0.99, max_treedepth = 16),
      backend = "cmdstanr",
      seed = 16) -> m_minimpute_sv_shannon_low_count_brms
m_minimpute_sv_shannon_low_count_brms

rstan::check_hmc_diagnostics(m_minimpute_sv_shannon_low_count_brms$fit)


#' model results
m_minimpute_sv_shannon_low_count_brms %>%
  brms::posterior_summary() %>%
  as_tibble(rownames = "param") %>%
  mutate_if(is.numeric, list("OR" = ~ exp(.x))) %>%
  identity() -> m_minimpute_sv_shannon_low_count_param_est
m_minimpute_sv_shannon_low_count_param_est
}



#' mixed model plot
if(run_plots) {
m_minimpute_sv_shannon_low_count_brms$data %>%
  as_tibble() %>%
  expand(shannon_low_count = modelr::seq_range(shannon_low_count, n = 100), subject_id = unique(subject_id)) %>%
  tidybayes::add_fitted_draws(newdata = ., model = m_minimpute_sv_shannon_low_count_brms) %>%
  ungroup() %>%
  ggplot(data = ., aes(x = shannon_low_count, y = .value)) +
  tidybayes::stat_lineribbon(.width = c(0.5, 0.8, 0.95), alpha = 0.7, color = colorspace::sequential_hcl(5, palette = "Blues 3")[1]) +
  colorspace::scale_fill_discrete_sequential(palette = "Blues 3", nmax = 5, order = 2:4) +
  facet_wrap(~ subject_id)
}








#' #######################################################################
#' SUSPECTED VAP + CX + SIGNS (same as SUSPECTED VAP + CX)
#' #######################################################################

if(run_models) {
dat_i %>%
  filter(suspected_vap == TRUE) %>%
  bind_rows(dat) %>%
  arrange(subject_id, date) %>%
  # filter to include dates only up to VAP (or last date, if no VAP)
  group_by(subject_id) %>%
  mutate(first_suspected_vap_plus_pathogen = ifelse(is.na(first_suspected_vap_plus_pathogen), max(date), first_suspected_vap_plus_pathogen)) %>%
  filter(date <= first_suspected_vap_plus_pathogen) %>%
  ungroup() %>%
  # # 48-hour window for case period
  # group_by(subject_id) %>%
  # arrange(date) %>%
  # mutate(suspected_vap_plus_pathogen = date == first_suspected_vap_plus_pathogen | date == first_suspected_vap_plus_pathogen - 1) %>%
  # ungroup() %>%
  # #
  select(subject_id, subject_day, shannon, suspected_vap_plus_signs_and_pathogen) %>%
  filter(complete.cases(.)) %>%
  mutate(shannon_less_than_2 = shannon < 2) %>%
  group_by(subject_id) %>%
  arrange(subject_day) %>%
  mutate(shannon_improved = shannon_less_than_2 == FALSE & lag(shannon_less_than_2, 1) == TRUE) %>%
  replace_na(replace = list("shannon_improved" = FALSE)) %>%
  mutate(shannon_flip_minimpute_count = as.character(cumsum(shannon_improved)),
  ) %>%
  tidyr::fill(shannon_flip_minimpute_count, .direction = "down") %>% #arrange(desc(shannon_flip_minimpute_count)) %>%
  group_by(subject_id, shannon_flip_minimpute_count) %>%
  arrange(subject_day) %>%
  #mutate(period_day = subject_day - min(subject_day, na.rm = TRUE)) %>%
  mutate(shannon_low_count = cumsum(shannon_less_than_2),
         #shannon_low_days = ifelse(shannon_low_count > 0, max(period_day[shannon_less_than_2], na.rm = TRUE) + 1, 0)
  ) %>%
  ungroup() %>% #qplot(data = ., x = shannon_low_count)
  #qplot(data = ., x = shannon_low_count, shannon_low_days)
  distinct() %>%
  filter(complete.cases(.)) %>%
  identity() -> d_vap

d_vap


#' run mixed effects model with data
d_vap %>%
  brm(suspected_vap_plus_signs_and_pathogen ~ shannon_low_count + 0 + (1|subject_id),
      data = .,
      family = "bernoulli",
      iter = 2000,
      warmup = 1000,
      chains = 4,
      cores = 4,
      control = list("adapt_delta" = 0.99, max_treedepth = 16),
      backend = "cmdstanr",
      seed = 16) -> m_minimpute_svp_minimpute_shannon_low_count_brms
m_minimpute_svp_minimpute_shannon_low_count_brms

rstan::check_hmc_diagnostics(m_minimpute_svp_minimpute_shannon_low_count_brms$fit)

#' model results
m_minimpute_svp_minimpute_shannon_low_count_brms %>%
  brms::posterior_summary() %>%
  as_tibble(rownames = "param") %>%
  mutate_if(is.numeric, list("OR" = ~ exp(.x))) %>%
  identity() -> m_minimpute_svp_minimpute_shannon_low_count_param_est
m_minimpute_svp_minimpute_shannon_low_count_param_est
}



#' mixed model plot
if(run_plots) {
m_minimpute_svp_minimpute_shannon_low_count_brms$data %>%
  as_tibble() %>%
  expand(shannon_low_count = modelr::seq_range(shannon_low_count, n = 100), subject_id = unique(subject_id)) %>%
  tidybayes::add_fitted_draws(newdata = ., model = m_minimpute_svp_minimpute_shannon_low_count_brms) %>%
  ungroup() %>%
  ggplot(data = ., aes(x = shannon_low_count, y = .value)) +
  tidybayes::stat_lineribbon(.width = c(0.5, 0.8, 0.95), alpha = 0.7, color = colorspace::sequential_hcl(5, palette = "Blues 3")[1]) +
  colorspace::scale_fill_discrete_sequential(palette = "Blues 3", nmax = 5, order = 2:4) +
  facet_wrap(~ subject_id)
}






#' #######################################################################
#' SUSPECTED VAP + ABX and PATHOGEN and SIGNS
#' #######################################################################

if(run_models) {
dat_i %>%
  filter(suspected_vap == TRUE) %>%
  bind_rows(dat) %>%
  arrange(subject_id, date) %>%
  # filter to include dates only up to VAP (or last date, if no VAP)
  group_by(subject_id) %>%
  mutate(first_suspected_vap_plus_abx_and_pathogen_and_signs = ifelse(is.na(first_suspected_vap_plus_abx_and_pathogen_and_signs), max(date), first_suspected_vap_plus_abx_and_pathogen_and_signs)) %>%
  filter(date <= first_suspected_vap_plus_abx_and_pathogen_and_signs) %>%
  ungroup() %>%
  # # 48-hour window for case period
  # group_by(subject_id) %>%
  # arrange(date) %>%
  # mutate(suspected_vap_plus_abx_and_pathogen_and_signs = date == first_suspected_vap_plus_abx_and_pathogen_and_signs | date == first_suspected_vap_plus_abx_and_pathogen_and_signs - 1) %>%
  # ungroup() %>%
  # #
  select(subject_id, subject_day, shannon, suspected_vap_plus_abx_and_pathogen_and_signs) %>%
  filter(complete.cases(.)) %>%
  mutate(shannon_less_than_2 = shannon < 2) %>%
  group_by(subject_id) %>%
  arrange(subject_day) %>%
  mutate(shannon_improved = shannon_less_than_2 == FALSE & lag(shannon_less_than_2, 1) == TRUE) %>%
  replace_na(replace = list("shannon_improved" = FALSE)) %>%
  mutate(shannon_flip_minimpute_count = as.character(cumsum(shannon_improved)),
  ) %>%
  tidyr::fill(shannon_flip_minimpute_count, .direction = "down") %>% #arrange(desc(shannon_flip_minimpute_count)) %>%
  group_by(subject_id, shannon_flip_minimpute_count) %>%
  arrange(subject_day) %>%
  #mutate(period_day = subject_day - min(subject_day, na.rm = TRUE)) %>%
  mutate(shannon_low_count = cumsum(shannon_less_than_2),
         #shannon_low_days = ifelse(shannon_low_count > 0, max(period_day[shannon_less_than_2], na.rm = TRUE) + 1, 0)
  ) %>%
  ungroup() %>% #qplot(data = ., x = shannon_low_count)
  #qplot(data = ., x = shannon_low_count, shannon_low_days)
  distinct() %>%
  filter(complete.cases(.)) %>%
  identity() -> d_vap

d_vap

#' run mixed effects model with data
d_vap %>%
  brm(suspected_vap_plus_abx_and_pathogen_and_signs ~ shannon_low_count + 0 + (1|subject_id),
      data = .,
      family = "bernoulli",
      iter = 2000,
      warmup = 1000,
      chains = 4,
      cores = 4,
      control = list("adapt_delta" = 0.99, max_treedepth = 16),
      backend = "cmdstanr",
      seed = 16) -> m_minimpute_svaps_shannon_low_count_brms
m_minimpute_svaps_shannon_low_count_brms

rstan::check_hmc_diagnostics(m_minimpute_svaps_shannon_low_count_brms$fit)


#' model results
m_minimpute_svaps_shannon_low_count_brms %>%
  brms::posterior_summary() %>%
  as_tibble(rownames = "param") %>%
  mutate_if(is.numeric, list("OR" = ~ exp(.x))) %>%
  identity() -> m_minimpute_svaps_shannon_low_count_param_est
m_minimpute_svaps_shannon_low_count_param_est
}



#' mixed model plot
if(run_plots) {
m_minimpute_svaps_shannon_low_count_brms$data %>%
  as_tibble() %>%
  expand(shannon_low_count = modelr::seq_range(shannon_low_count, n = 100), subject_id = unique(subject_id)) %>%
  tidybayes::add_fitted_draws(newdata = ., model = m_minimpute_svaps_shannon_low_count_brms) %>%
  ungroup() %>%
  ggplot(data = ., aes(x = shannon_low_count, y = .value)) +
  tidybayes::stat_lineribbon(.width = c(0.5, 0.8, 0.95), alpha = 0.7, color = colorspace::sequential_hcl(5, palette = "Blues 3")[1]) +
  colorspace::scale_fill_discrete_sequential(palette = "Blues 3", nmax = 5, order = 2:4) +
  facet_wrap(~ subject_id)
}





#' #######################################################################
#' #######################################################################
#' VAP as a function of lower respiratory tract absolute bacterial abundance (16S QPCR)
#' #######################################################################
#' #######################################################################


#' #######################################################################
#' SUSPECTED VAP
#' #######################################################################

if(run_models) {
dat_i %>%
  filter(suspected_vap == TRUE) %>%
  bind_rows(dat) %>%
  arrange(subject_id, date) %>%
  # filter to include dates only up to VAP (or last date, if no VAP)
  group_by(subject_id) %>%
  mutate(first_suspected_vap = ifelse(is.na(first_suspected_vap), max(date), first_suspected_vap)) %>%
  filter(date <= first_suspected_vap) %>%
  ungroup() %>%
  # # 48-hour window for case period
  # group_by(subject_id) %>%
  # arrange(date) %>%
  # mutate(suspected_vap = date == first_suspected_vap | date == first_suspected_vap - 1) %>%
  # ungroup() %>%
  # #
  select(subject_id, subject_day, log_copy_16S_per_ml_sputum, suspected_vap) %>%
  filter(complete.cases(.)) %>%
  #summarise(summary(log_copy_16S_per_ml_sputum))
  mutate(log_copy_16S_per_ml_sputum_greater_than_9 = log_copy_16S_per_ml_sputum > 9) %>%
  group_by(subject_id) %>%
  arrange(subject_day) %>%
  mutate(log_copy_16S_per_ml_sputum_improved = log_copy_16S_per_ml_sputum_greater_than_9 == FALSE & lag(log_copy_16S_per_ml_sputum_greater_than_9, 1) == TRUE) %>%
  replace_na(replace = list("log_copy_16S_per_ml_sputum_improved" = FALSE)) %>%
  mutate(log_copy_16S_per_ml_sputum_flip_minimpute_count = as.character(cumsum(log_copy_16S_per_ml_sputum_improved)),
  ) %>%
  tidyr::fill(log_copy_16S_per_ml_sputum_flip_minimpute_count, .direction = "down") %>% #arrange(desc(log_copy_16S_per_ml_sputum_flip_minimpute_count)) %>%
  group_by(subject_id, log_copy_16S_per_ml_sputum_flip_minimpute_count) %>%
  arrange(subject_day) %>%
  mutate(log_copy_16S_per_ml_sputum_high_count = cumsum(log_copy_16S_per_ml_sputum_greater_than_9)) %>%
  ungroup() %>% #qplot(data = ., x = log_copy_16S_per_ml_sputum_high_count)
  distinct() %>%
  filter(complete.cases(.)) %>%
  identity() -> d_vap

d_vap





#' run mixed effects model with data
d_vap %>%
  brm(suspected_vap ~ log_copy_16S_per_ml_sputum_high_count + 0 + (1|subject_id),
      data = .,
      family = "bernoulli",
      iter = 2000,
      warmup = 1000,
      chains = 4,
      cores = 4,
      control = list("adapt_delta" = 0.99, max_treedepth = 16),
      backend = "cmdstanr",
      seed = 16) -> m_minimpute_sv_log_copy_16S_per_ml_sputum_high_count_brms
m_minimpute_sv_log_copy_16S_per_ml_sputum_high_count_brms

rstan::check_hmc_diagnostics(m_minimpute_sv_log_copy_16S_per_ml_sputum_high_count_brms$fit)


#' model results
m_minimpute_sv_log_copy_16S_per_ml_sputum_high_count_brms %>%
  brms::posterior_summary() %>%
  as_tibble(rownames = "param") %>%
  mutate_if(is.numeric, list("OR" = ~ exp(.x))) %>%
  identity() -> m_minimpute_sv_log_copy_16S_per_ml_sputum_high_count_param_est
m_minimpute_sv_log_copy_16S_per_ml_sputum_high_count_param_est
}



#' mixed model plot
if(run_plots) {
m_minimpute_sv_log_copy_16S_per_ml_sputum_high_count_brms$data %>%
  as_tibble() %>%
  expand(log_copy_16S_per_ml_sputum_high_count = modelr::seq_range(log_copy_16S_per_ml_sputum_high_count, n = 100), subject_id = unique(subject_id)) %>%
  tidybayes::add_fitted_draws(newdata = ., model = m_minimpute_sv_log_copy_16S_per_ml_sputum_high_count_brms) %>%
  ungroup() %>%
  ggplot(data = ., aes(x = log_copy_16S_per_ml_sputum_high_count, y = .value)) +
  tidybayes::stat_lineribbon(.width = c(0.5, 0.8, 0.95), alpha = 0.7, color = colorspace::sequential_hcl(5, palette = "Blues 3")[1]) +
  colorspace::scale_fill_discrete_sequential(palette = "Blues 3", nmax = 5, order = 2:4) +
  facet_wrap(~ subject_id)
}









#' #######################################################################
#' SUSPECTED VAP + CX + SIGNS (same as SUSPECTED VAP + CX)
#' #######################################################################

if(run_models) {
dat_i %>%
  filter(suspected_vap == TRUE) %>%
  bind_rows(dat) %>%
  arrange(subject_id, date) %>%
  # filter to include dates only up to VAP (or last date, if no VAP)
  group_by(subject_id) %>%
  mutate(first_suspected_vap_plus_pathogen = ifelse(is.na(first_suspected_vap_plus_pathogen), max(date), first_suspected_vap_plus_pathogen)) %>%
  filter(date <= first_suspected_vap_plus_pathogen) %>%
  ungroup() %>%
  # # 48-hour window for case period
  # group_by(subject_id) %>%
  # arrange(date) %>%
  # mutate(suspected_vap_plus_signs_and_pathogen = date == first_suspected_vap_plus_pathogen | date == first_suspected_vap_plus_pathogen - 1) %>%
  # ungroup() %>%
  # #
  select(subject_id, subject_day, log_copy_16S_per_ml_sputum, suspected_vap_plus_signs_and_pathogen) %>%
  filter(complete.cases(.)) %>%
  mutate(log_copy_16S_per_ml_sputum_greater_than_9 = log_copy_16S_per_ml_sputum > 9) %>%
  group_by(subject_id) %>%
  arrange(subject_day) %>%
  mutate(log_copy_16S_per_ml_sputum_improved = log_copy_16S_per_ml_sputum_greater_than_9 == FALSE & lag(log_copy_16S_per_ml_sputum_greater_than_9, 1) == TRUE) %>%
  replace_na(replace = list("log_copy_16S_per_ml_sputum_improved" = FALSE)) %>%
  mutate(log_copy_16S_per_ml_sputum_flip_minimpute_count = as.character(cumsum(log_copy_16S_per_ml_sputum_improved)),
  ) %>%
  tidyr::fill(log_copy_16S_per_ml_sputum_flip_minimpute_count, .direction = "down") %>% #arrange(desc(log_copy_16S_per_ml_sputum_flip_minimpute_count)) %>%
  group_by(subject_id, log_copy_16S_per_ml_sputum_flip_minimpute_count) %>%
  arrange(subject_day) %>%
  mutate(log_copy_16S_per_ml_sputum_high_count = cumsum(log_copy_16S_per_ml_sputum_greater_than_9)) %>%
  ungroup() %>% #qplot(data = ., x = log_copy_16S_per_ml_sputum_high_count)
  distinct() %>%
  filter(complete.cases(.)) %>%
  identity() -> d_vap

d_vap



#' run mixed effects model with data
d_vap %>%
  brm(suspected_vap_plus_signs_and_pathogen ~ log_copy_16S_per_ml_sputum_high_count + 0 + (1|subject_id),
      data = .,
      family = "bernoulli",
      iter = 2000,
      warmup = 1000,
      chains = 4,
      cores = 4,
      control = list("adapt_delta" = 0.99, max_treedepth = 16),
      backend = "cmdstanr",
      seed = 16) -> m_minimpute_svp_minimpute_log_copy_16S_per_ml_sputum_high_count_brms
m_minimpute_svp_minimpute_log_copy_16S_per_ml_sputum_high_count_brms

rstan::check_hmc_diagnostics(m_minimpute_svp_minimpute_log_copy_16S_per_ml_sputum_high_count_brms$fit)



#' model results
m_minimpute_svp_minimpute_log_copy_16S_per_ml_sputum_high_count_brms %>%
  brms::posterior_summary() %>%
  as_tibble(rownames = "param") %>%
  mutate_if(is.numeric, list("OR" = ~ exp(.x))) %>%
  identity() -> m_minimpute_svp_minimpute_log_copy_16S_per_ml_sputum_high_count_param_est
m_minimpute_svp_minimpute_log_copy_16S_per_ml_sputum_high_count_param_est
}



#' mixed model plot
if(run_plots) {
m_minimpute_svp_minimpute_log_copy_16S_per_ml_sputum_high_count_brms$data %>%
  as_tibble() %>%
  expand(log_copy_16S_per_ml_sputum_high_count = modelr::seq_range(log_copy_16S_per_ml_sputum_high_count, n = 100), subject_id = unique(subject_id)) %>%
  tidybayes::add_fitted_draws(newdata = ., model = m_minimpute_svp_minimpute_log_copy_16S_per_ml_sputum_high_count_brms) %>%
  ungroup() %>%
  ggplot(data = ., aes(x = log_copy_16S_per_ml_sputum_high_count, y = .value)) +
  tidybayes::stat_lineribbon(.width = c(0.5, 0.8, 0.95), alpha = 0.7, color = colorspace::sequential_hcl(5, palette = "Blues 3")[1]) +
  colorspace::scale_fill_discrete_sequential(palette = "Blues 3", nmax = 5, order = 2:4) +
  facet_wrap(~ subject_id)
}





#' #######################################################################
#' SUSPECTED VAP + ABX and PATHOGEN and SIGNS
#' #######################################################################

if(run_models) {
dat_i %>%
  filter(suspected_vap == TRUE) %>%
  bind_rows(dat) %>%
  arrange(subject_id, date) %>%
  # filter to include dates only up to VAP (or last date, if no VAP)
  group_by(subject_id) %>%
  mutate(first_suspected_vap_plus_abx_and_pathogen_and_signs = ifelse(is.na(first_suspected_vap_plus_abx_and_pathogen_and_signs), max(date), first_suspected_vap_plus_abx_and_pathogen_and_signs)) %>%
  filter(date <= first_suspected_vap_plus_abx_and_pathogen_and_signs) %>%
  ungroup() %>%
  # # 48-hour window for case period
  # group_by(subject_id) %>%
  # arrange(date) %>%
  # mutate(suspected_vap_plus_abx_and_pathogen_and_signs = date == first_suspected_vap_plus_abx_and_pathogen_and_signs | date == first_suspected_vap_plus_abx_and_pathogen_and_signs - 1) %>%
  # ungroup() %>%
  # #
  select(subject_id, subject_day, log_copy_16S_per_ml_sputum, suspected_vap_plus_abx_and_pathogen_and_signs) %>%
  filter(complete.cases(.)) %>%
  mutate(log_copy_16S_per_ml_sputum_greater_than_9 = log_copy_16S_per_ml_sputum > 9) %>%
  group_by(subject_id) %>%
  arrange(subject_day) %>%
  mutate(log_copy_16S_per_ml_sputum_improved = log_copy_16S_per_ml_sputum_greater_than_9 == FALSE & lag(log_copy_16S_per_ml_sputum_greater_than_9, 1) == TRUE) %>%
  replace_na(replace = list("log_copy_16S_per_ml_sputum_improved" = FALSE)) %>%
  mutate(log_copy_16S_per_ml_sputum_flip_minimpute_count = as.character(cumsum(log_copy_16S_per_ml_sputum_improved)),
  ) %>%
  tidyr::fill(log_copy_16S_per_ml_sputum_flip_minimpute_count, .direction = "down") %>% #arrange(desc(log_copy_16S_per_ml_sputum_flip_minimpute_count)) %>%
  group_by(subject_id, log_copy_16S_per_ml_sputum_flip_minimpute_count) %>%
  arrange(subject_day) %>%
  mutate(log_copy_16S_per_ml_sputum_high_count = cumsum(log_copy_16S_per_ml_sputum_greater_than_9)) %>%
  ungroup() %>% #qplot(data = ., x = log_copy_16S_per_ml_sputum_high_count)
  distinct() %>%
  filter(complete.cases(.)) %>%
  identity() -> d_vap

d_vap



#' run mixed effects model with data
d_vap %>%
  brm(suspected_vap_plus_abx_and_pathogen_and_signs ~ log_copy_16S_per_ml_sputum_high_count + 0 + (1|subject_id),
      data = .,
      family = "bernoulli",
      iter = 2000,
      warmup = 1000,
      chains = 4,
      cores = 4,
      control = list("adapt_delta" = 0.99, max_treedepth = 16),
      backend = "cmdstanr",
      seed = 16) -> m_minimpute_svaps_log_copy_16S_per_ml_sputum_high_count_brms
m_minimpute_svaps_log_copy_16S_per_ml_sputum_high_count_brms

rstan::check_hmc_diagnostics(m_minimpute_svaps_log_copy_16S_per_ml_sputum_high_count_brms$fit)




#' model results
m_minimpute_svaps_log_copy_16S_per_ml_sputum_high_count_brms %>%
  brms::posterior_summary() %>%
  as_tibble(rownames = "param") %>%
  mutate_if(is.numeric, list("OR" = ~ exp(.x))) %>%
  identity() -> m_minimpute_svaps_log_copy_16S_per_ml_sputum_high_count_param_est
m_minimpute_svaps_log_copy_16S_per_ml_sputum_high_count_param_est
}



#' mixed model plot
if(run_plots) {
m_minimpute_svaps_log_copy_16S_per_ml_sputum_high_count_brms$data %>%
  as_tibble() %>%
  expand(log_copy_16S_per_ml_sputum_high_count = modelr::seq_range(log_copy_16S_per_ml_sputum_high_count, n = 100), subject_id = unique(subject_id)) %>%
  tidybayes::add_fitted_draws(newdata = ., model = m_minimpute_svaps_log_copy_16S_per_ml_sputum_high_count_brms) %>%
  ungroup() %>%
  ggplot(data = ., aes(x = log_copy_16S_per_ml_sputum_high_count, y = .value)) +
  tidybayes::stat_lineribbon(.width = c(0.5, 0.8, 0.95), alpha = 0.7, color = colorspace::sequential_hcl(5, palette = "Blues 3")[1]) +
  colorspace::scale_fill_discrete_sequential(palette = "Blues 3", nmax = 5, order = 2:4) +
  facet_wrap(~ subject_id)
}







#' #######################################################################
#' #######################################################################
#' VAP as a function of maximum proportional abundance 
#' #######################################################################
#' #######################################################################


#' #######################################################################
#' SUSPECTED VAP
#' #######################################################################

if(run_models) {
dat_i %>%
  filter(suspected_vap == TRUE) %>%
  bind_rows(dat) %>%
  arrange(subject_id, date) %>%
  # filter to include dates only up to VAP (or last date, if no VAP)
  group_by(subject_id) %>%
  mutate(first_suspected_vap = ifelse(is.na(first_suspected_vap), max(date), first_suspected_vap)) %>%
  filter(date <= first_suspected_vap) %>%
  ungroup() %>%
  # # 48-hour window for case period
  # group_by(subject_id) %>%
  # arrange(date) %>%
  # mutate(suspected_vap = date == first_suspected_vap | date == first_suspected_vap - 1) %>%
  # ungroup() %>%
  # #
  select(subject_id, subject_day, max_read_prop, suspected_vap) %>%
  filter(complete.cases(.)) %>%
  #summarise(summary(max_read_prop))
  mutate(max_read_prop_minimpute_greater_than_0.5 = max_read_prop > 0.5) %>%
  group_by(subject_id) %>%
  arrange(subject_day) %>%
  mutate(max_read_prop_minimpute_improved = max_read_prop_minimpute_greater_than_0.5 == FALSE & lag(max_read_prop_minimpute_greater_than_0.5, 1) == TRUE) %>%
  replace_na(replace = list("max_read_prop_minimpute_improved" = FALSE)) %>%
  mutate(max_read_prop_minimpute_flip_minimpute_count = as.character(cumsum(max_read_prop_minimpute_improved)),
  ) %>%
  tidyr::fill(max_read_prop_minimpute_flip_minimpute_count, .direction = "down") %>% #arrange(desc(max_read_prop_minimpute_flip_minimpute_count)) %>%
  group_by(subject_id, max_read_prop_minimpute_flip_minimpute_count) %>%
  arrange(subject_day) %>%
  mutate(max_read_prop_minimpute_high_count = cumsum(max_read_prop_minimpute_greater_than_0.5)) %>%
  ungroup() %>% #qplot(data = ., x = max_read_prop_minimpute_high_count)
  distinct() %>%
  filter(complete.cases(.)) %>%
  identity() -> d_vap

d_vap




#' run mixed effects model with data
d_vap %>%
  brm(suspected_vap ~ max_read_prop_minimpute_high_count + 0 + (1|subject_id),
      data = .,
      family = "bernoulli",
      iter = 2000,
      warmup = 1000,
      chains = 4,
      cores = 4,
      control = list("adapt_delta" = 0.99, max_treedepth = 16),
      backend = "cmdstanr",
      seed = 16) -> m_minimpute_sv_max_read_prop_minimpute_high_count_brms
m_minimpute_sv_max_read_prop_minimpute_high_count_brms

rstan::check_hmc_diagnostics(m_minimpute_sv_max_read_prop_minimpute_high_count_brms$fit)


#' model results
m_minimpute_sv_max_read_prop_minimpute_high_count_brms %>%
  brms::posterior_summary() %>%
  as_tibble(rownames = "param") %>%
  mutate_if(is.numeric, list("OR" = ~ exp(.x))) %>%
  identity() -> m_minimpute_sv_max_read_prop_minimpute_high_count_param_est
m_minimpute_sv_max_read_prop_minimpute_high_count_param_est
}



#' mixed model plot
if(run_plots) {
m_minimpute_sv_max_read_prop_minimpute_high_count_brms$data %>%
  as_tibble() %>%
  expand(max_read_prop_minimpute_high_count = modelr::seq_range(max_read_prop_minimpute_high_count, n = 100), subject_id = unique(subject_id)) %>%
  tidybayes::add_fitted_draws(newdata = ., model = m_minimpute_sv_max_read_prop_minimpute_high_count_brms) %>%
  ungroup() %>%
  ggplot(data = ., aes(x = max_read_prop_minimpute_high_count, y = .value)) +
  tidybayes::stat_lineribbon(.width = c(0.5, 0.8, 0.95), alpha = 0.7, color = colorspace::sequential_hcl(5, palette = "Blues 3")[1]) +
  colorspace::scale_fill_discrete_sequential(palette = "Blues 3", nmax = 5, order = 2:4) +
  facet_wrap(~ subject_id)
}









#' #######################################################################
#' SUSPECTED VAP + CX + SIGNS (same as SUSPECTED VAP + CX)
#' #######################################################################

if(run_models) {
dat_i %>%
  filter(suspected_vap == TRUE) %>%
  bind_rows(dat) %>%
  arrange(subject_id, date) %>%
  # filter to include dates only up to VAP (or last date, if no VAP)
  group_by(subject_id) %>%
  mutate(first_suspected_vap_plus_pathogen = ifelse(is.na(first_suspected_vap_plus_pathogen), max(date), first_suspected_vap_plus_pathogen)) %>%
  filter(date <= first_suspected_vap_plus_pathogen) %>%
  ungroup() %>%
  # # 48-hour window for case period
  # group_by(subject_id) %>%
  # arrange(date) %>%
  # mutate(suspected_vap_plus_pathogen = date == first_suspected_vap_plus_pathogen | date == first_suspected_vap_plus_pathogen - 1) %>%
  # ungroup() %>%
  # #
  select(subject_id, subject_day, max_read_prop, suspected_vap_plus_signs_and_pathogen) %>%
  filter(complete.cases(.)) %>%
  mutate(max_read_prop_minimpute_greater_than_0.5 = max_read_prop > 0.5) %>%
  group_by(subject_id) %>%
  arrange(subject_day) %>%
  mutate(max_read_prop_minimpute_improved = max_read_prop_minimpute_greater_than_0.5 == FALSE & lag(max_read_prop_minimpute_greater_than_0.5, 1) == TRUE) %>%
  replace_na(replace = list("max_read_prop_minimpute_improved" = FALSE)) %>%
  mutate(max_read_prop_minimpute_flip_minimpute_count = as.character(cumsum(max_read_prop_minimpute_improved)),
  ) %>%
  tidyr::fill(max_read_prop_minimpute_flip_minimpute_count, .direction = "down") %>% #arrange(desc(max_read_prop_minimpute_flip_minimpute_count)) %>%
  group_by(subject_id, max_read_prop_minimpute_flip_minimpute_count) %>%
  arrange(subject_day) %>%
  mutate(max_read_prop_minimpute_high_count = cumsum(max_read_prop_minimpute_greater_than_0.5)) %>%
  ungroup() %>% #qplot(data = ., x = max_read_prop_minimpute_high_count)
  distinct() %>%
  filter(complete.cases(.)) %>%
  identity() -> d_vap

d_vap




#' run mixed effects model with data
d_vap %>%
  brm(suspected_vap_plus_signs_and_pathogen ~ max_read_prop_minimpute_high_count + 0 + (1|subject_id),
      data = .,
      family = "bernoulli",
      iter = 2000,
      warmup = 1000,
      chains = 4,
      cores = 4,
      control = list("adapt_delta" = 0.99, max_treedepth = 16),
      backend = "cmdstanr",
      seed = 16) -> m_minimpute_svp_minimpute_max_read_prop_minimpute_high_count_brms
m_minimpute_svp_minimpute_max_read_prop_minimpute_high_count_brms

rstan::check_hmc_diagnostics(m_minimpute_svp_minimpute_max_read_prop_minimpute_high_count_brms$fit)



#' model results
m_minimpute_svp_minimpute_max_read_prop_minimpute_high_count_brms %>%
  brms::posterior_summary() %>%
  as_tibble(rownames = "param") %>%
  mutate_if(is.numeric, list("OR" = ~ exp(.x))) %>%
  identity() -> m_minimpute_svp_minimpute_max_read_prop_minimpute_high_count_param_est
m_minimpute_svp_minimpute_max_read_prop_minimpute_high_count_param_est
}



#' mixed model plot
if(run_plots) {
m_minimpute_svp_minimpute_max_read_prop_minimpute_high_count_brms$data %>%
  as_tibble() %>%
  expand(max_read_prop_minimpute_high_count = modelr::seq_range(max_read_prop_minimpute_high_count, n = 100), subject_id = unique(subject_id)) %>%
  tidybayes::add_fitted_draws(newdata = ., model = m_minimpute_svp_minimpute_max_read_prop_minimpute_high_count_brms) %>%
  ungroup() %>%
  ggplot(data = ., aes(x = max_read_prop_minimpute_high_count, y = .value)) +
  tidybayes::stat_lineribbon(.width = c(0.5, 0.8, 0.95), alpha = 0.7, color = colorspace::sequential_hcl(5, palette = "Blues 3")[1]) +
  colorspace::scale_fill_discrete_sequential(palette = "Blues 3", nmax = 5, order = 2:4) +
  facet_wrap(~ subject_id)
}




#' #######################################################################
#' SUSPECTED VAP + ABX and PATHOGEN and SIGNS
#' #######################################################################

if(run_models) {
dat_i %>%
  filter(suspected_vap == TRUE) %>%
  bind_rows(dat) %>%
  arrange(subject_id, date) %>%
  # filter to include dates only up to VAP (or last date, if no VAP)
  group_by(subject_id) %>%
  mutate(first_suspected_vap_plus_abx_and_pathogen_and_signs = ifelse(is.na(first_suspected_vap_plus_abx_and_pathogen_and_signs), max(date), first_suspected_vap_plus_abx_and_pathogen_and_signs)) %>%
  filter(date <= first_suspected_vap_plus_abx_and_pathogen_and_signs) %>%
  ungroup() %>%
  # # 48-hour window for case period
  # group_by(subject_id) %>%
  # arrange(date) %>%
  # mutate(suspected_vap_plus_abx_and_pathogen_and_signs = date == first_suspected_vap_plus_abx_and_pathogen_and_signs | date == first_suspected_vap_plus_abx_and_pathogen_and_signs - 1) %>%
  # ungroup() %>%
  # #
  select(subject_id, subject_day, max_read_prop, suspected_vap_plus_abx_and_pathogen_and_signs) %>%
  filter(complete.cases(.)) %>%
  mutate(max_read_prop_minimpute_greater_than_0.5 = max_read_prop > 0.5) %>%
  group_by(subject_id) %>%
  arrange(subject_day) %>%
  mutate(max_read_prop_minimpute_improved = max_read_prop_minimpute_greater_than_0.5 == FALSE & lag(max_read_prop_minimpute_greater_than_0.5, 1) == TRUE) %>%
  replace_na(replace = list("max_read_prop_minimpute_improved" = FALSE)) %>%
  mutate(max_read_prop_minimpute_flip_minimpute_count = as.character(cumsum(max_read_prop_minimpute_improved)),
  ) %>%
  tidyr::fill(max_read_prop_minimpute_flip_minimpute_count, .direction = "down") %>% #arrange(desc(max_read_prop_minimpute_flip_minimpute_count)) %>%
  group_by(subject_id, max_read_prop_minimpute_flip_minimpute_count) %>%
  arrange(subject_day) %>%
  mutate(max_read_prop_minimpute_high_count = cumsum(max_read_prop_minimpute_greater_than_0.5)) %>%
  ungroup() %>% #qplot(data = ., x = max_read_prop_minimpute_high_count)
  distinct() %>%
  filter(complete.cases(.)) %>%
  identity() -> d_vap

d_vap



#' run mixed effects model with data
d_vap %>%
  brm(suspected_vap_plus_abx_and_pathogen_and_signs ~ max_read_prop_minimpute_high_count + 0 + (1|subject_id),
      data = .,
      family = "bernoulli",
      iter = 2000,
      warmup = 1000,
      chains = 4,
      cores = 4,
      control = list("adapt_delta" = 0.99, max_treedepth = 16),
      backend = "cmdstanr",
      seed = 16) -> m_minimpute_svaps_max_read_prop_minimpute_high_count_brms
m_minimpute_svaps_max_read_prop_minimpute_high_count_brms

rstan::check_hmc_diagnostics(m_minimpute_svaps_max_read_prop_minimpute_high_count_brms$fit)




#' model results
m_minimpute_svaps_max_read_prop_minimpute_high_count_brms %>%
  brms::posterior_summary() %>%
  as_tibble(rownames = "param") %>%
  mutate_if(is.numeric, list("OR" = ~ exp(.x))) %>%
  identity() -> m_minimpute_svaps_max_read_prop_minimpute_high_count_param_est
m_minimpute_svaps_max_read_prop_minimpute_high_count_param_est
}



#' mixed model plot
if(run_plots) {
m_minimpute_svaps_max_read_prop_minimpute_high_count_brms$data %>%
  as_tibble() %>%
  expand(max_read_prop_minimpute_high_count = modelr::seq_range(max_read_prop_minimpute_high_count, n = 100), subject_id = unique(subject_id)) %>%
  tidybayes::add_fitted_draws(newdata = ., model = m_minimpute_svaps_max_read_prop_minimpute_high_count_brms) %>%
  ungroup() %>%
  ggplot(data = ., aes(x = max_read_prop_minimpute_high_count, y = .value)) +
  tidybayes::stat_lineribbon(.width = c(0.5, 0.8, 0.95), alpha = 0.7, color = colorspace::sequential_hcl(5, palette = "Blues 3")[1]) +
  colorspace::scale_fill_discrete_sequential(palette = "Blues 3", nmax = 5, order = 2:4) +
  facet_wrap(~ subject_id)
}







#' #################################################
#' #################################################
#' 
#' COMBINE PARAMETER ESTIMATES FROM MULTIPLE MODELS
#' 
#' #################################################
#' #################################################

tibble(value = c("m_minimpute_sv_shannon_low_count_param_est",
                 "m_minimpute_svp_minimpute_shannon_low_count_param_est",
                 "m_minimpute_svaps_shannon_low_count_param_est",
                 "m_minimpute_sv_max_read_prop_minimpute_high_count_param_est",
                 "m_minimpute_svp_minimpute_max_read_prop_minimpute_high_count_param_est",
                 "m_minimpute_svaps_max_read_prop_minimpute_high_count_param_est",
                 "m_minimpute_sv_log_copy_16S_per_ml_sputum_high_count_param_est",
                 "m_minimpute_svp_minimpute_log_copy_16S_per_ml_sputum_high_count_param_est",
                 "m_minimpute_svaps_log_copy_16S_per_ml_sputum_high_count_param_est")) %>%
  mutate(parameter = stringr::str_extract(string = value, pattern = "log_copy_16S_per_ml_sputum|shannon|max_read_prop")) %>%
  mutate(tib = map(.x = value, .f = ~ get(.x))) %>%
  unnest(cols = c(tib)) %>% #View()
  filter(map2_lgl(.x = parameter, .y = param, .f = ~ grepl(pattern = .x, x = .y))) %>%
  identity() -> m_minimpute_mixed_longitudinal_combined_parameter_estimates
m_minimpute_mixed_longitudinal_combined_parameter_estimates

if(save_models) {
m_minimpute_mixed_longitudinal_combined_parameter_estimates %>%
  write_csv(file = "./models/longitudinal/m_minimpute_mixed_longitudinal_combined_parameter_estimates.csv")
}

#' compare paramaters on OR scale in a table
m_minimpute_mixed_longitudinal_combined_parameter_estimates %>%
  mutate(outcome = stringr::str_extract(string = value, pattern = "sv_|svp|svaps")) %>%
  mutate(outcome = case_when(outcome == "sv_" ~ "Suspected VA-LRTI",
                             outcome == "svp" ~ "Suspected VA-LRTI\n+Culture\n+Signs",
                             outcome == "svaps" ~ "Suspected VA-LRTI\n+Culture\n+Signs\n+Antibiotic Change")) %>%
  select(outcome, parameter, contains("_OR")) %>%
  gt::gt() %>%
  gt::fmt_number(columns = 3:6, n_sigfig = 3) %>%
  gt::tab_header(title = "MDI Persistence Models")



#' #################################################
#' #################################################
#' 
#' WRITE MODELS TO FILE
#' 
#' #################################################
#' #################################################

if(save_models) {
  tibble(value = c("m_minimpute_sv_shannon_low_count_brms",
                   "m_minimpute_svp_minimpute_shannon_low_count_brms",
                   "m_minimpute_svaps_shannon_low_count_brms",
                   "m_minimpute_sv_max_read_prop_minimpute_high_count_brms",
                   "m_minimpute_svp_minimpute_max_read_prop_minimpute_high_count_brms",
                   "m_minimpute_svaps_max_read_prop_minimpute_high_count_brms",
                   "m_minimpute_sv_log_copy_16S_per_ml_sputum_high_count_brms",
                   "m_minimpute_svp_minimpute_log_copy_16S_per_ml_sputum_high_count_brms",
                   "m_minimpute_svaps_log_copy_16S_per_ml_sputum_high_count_brms")) %>%
  mutate(model = map(.x = value, .f = ~ get(.x))) %>%
  map2(.x = .$value, .y = .$model, .f = ~ write_rds(x = .y, file = paste0("./models/longitudinal/",.x,".rds.gz"), compress = "gz"))
}



#' #################################################
#' #################################################
#' 
#' READ MODELS FROM FILE
#' 
#' #################################################
#' #################################################

tibble(model_files = list.files(path = "./models/longitudinal/", pattern = "^m_minimpute_.*_brms\\.rds\\.gz", full.names = TRUE),
       model_names = stringr::str_remove(list.files(path = "./models/longitudinal/", pattern = "^m_minimpute_.*_brms\\.rds\\.gz", full.names = FALSE), pattern = "\\.rds\\.gz"),
) %>%
  mutate(models = map(.x = model_files, .f = ~ read_rds(.x))) %>%
  identity() -> persistence_models
persistence_models

persistence_models %>%
  map2(.x = .$model_names, .y = .$models, .f = ~ assign(.x, .y, envir = .GlobalEnv))



#' #################################################
#' #################################################
#' 
#' PLOT ENSEMBLE MODEL PARAMETERS (PROBABILITY & ODDS RATIO SCALES)
#' 
#' #################################################
#' #################################################

#' probability scale
tibble(value = c("m_minimpute_sv_shannon_low_count_brms",
                 "m_minimpute_svp_minimpute_shannon_low_count_brms",
                 "m_minimpute_svaps_shannon_low_count_brms",
                 "m_minimpute_sv_max_read_prop_minimpute_high_count_brms",
                 "m_minimpute_svp_minimpute_max_read_prop_minimpute_high_count_brms",
                 "m_minimpute_svaps_max_read_prop_minimpute_high_count_brms",
                 "m_minimpute_sv_log_copy_16S_per_ml_sputum_high_count_brms",
                 "m_minimpute_svp_minimpute_log_copy_16S_per_ml_sputum_high_count_brms",
                 "m_minimpute_svaps_log_copy_16S_per_ml_sputum_high_count_brms")) %>%
  mutate(parameter = stringr::str_extract(string = value, pattern = "log_copy_16S_per_ml_sputum|shannon|max_read_prop")) %>%
  mutate(model = map(.x = value, .f = ~ get(.x))) %>%
  mutate(fit = map(.x = model, .f = ~ .x$fit)) %>%
  mutate(param = gsub("m_minimpute_svp_minimpute_","b_",gsub("_brms","",value))) %>%
  mutate(posterior = map(.x = fit, .f = ~ posterior::as_draws_df(.x) %>% as_tibble() %>% select(contains("b_")) %>% pull())) %>%
  select(value, parameter, posterior) %>%
  unnest(cols = c(posterior)) %>%
  mutate(posterior_OR = exp(posterior)) %>%
  ggplot(data = ., aes(y = parameter, x = posterior, fill = stat(x > 0))) +
  tidybayes::stat_halfeye(alpha = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("gray80","skyblue"))



#' odds ratio scale
tibble(value = c("m_minimpute_sv_shannon_low_count_brms",
                 "m_minimpute_svp_minimpute_shannon_low_count_brms",
                 "m_minimpute_svaps_shannon_low_count_brms",
                 "m_minimpute_sv_max_read_prop_minimpute_high_count_brms",
                 "m_minimpute_svp_minimpute_max_read_prop_minimpute_high_count_brms",
                 "m_minimpute_svaps_max_read_prop_minimpute_high_count_brms",
                 "m_minimpute_sv_log_copy_16S_per_ml_sputum_high_count_brms",
                 "m_minimpute_svp_minimpute_log_copy_16S_per_ml_sputum_high_count_brms",
                 "m_minimpute_svaps_log_copy_16S_per_ml_sputum_high_count_brms")) %>%
  mutate(parameter = stringr::str_extract(string = value, pattern = "log_copy_16S_per_ml_sputum|shannon|max_read_prop")) %>%
  mutate(model = map(.x = value, .f = ~ get(.x))) %>%
  mutate(fit = map(.x = model, .f = ~ .x$fit)) %>%
  mutate(param = gsub("m_minimpute_svp_minimpute_","b_",gsub("_brms","",value))) %>%
  mutate(posterior = map(.x = fit, .f = ~ posterior::as_draws_df(.x) %>% as_tibble() %>% select(contains("b_")) %>% pull())) %>%
  select(value, parameter, posterior) %>%
  unnest(cols = c(posterior)) %>%
  mutate(posterior_OR = exp(posterior)) %>%
  ggplot(data = ., aes(y = parameter, x = posterior_OR, fill = stat(x > 1))) +
  tidybayes::stat_halfeye(alpha = 0.8) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_fill_manual(values = c("gray80","skyblue"))
  


#' #################################################
#' #################################################
#' 
#' PLOT SVP MODEL PARAMETERS (PROBABILITY & ODDS RATIO SCALES)
#' 
#' #################################################
#' #################################################


#' probability scale
tibble(value = c("m_minimpute_sv_shannon_low_count_brms",
                 "m_minimpute_svp_minimpute_shannon_low_count_brms",
                 "m_minimpute_svaps_shannon_low_count_brms",
                 "m_minimpute_sv_max_read_prop_minimpute_high_count_brms",
                 "m_minimpute_svp_minimpute_max_read_prop_minimpute_high_count_brms",
                 "m_minimpute_svaps_max_read_prop_minimpute_high_count_brms",
                 "m_minimpute_sv_log_copy_16S_per_ml_sputum_high_count_brms",
                 "m_minimpute_svp_minimpute_log_copy_16S_per_ml_sputum_high_count_brms",
                 "m_minimpute_svaps_log_copy_16S_per_ml_sputum_high_count_brms")) %>%
  mutate(parameter = stringr::str_extract(string = value, pattern = "log_copy_16S_per_ml_sputum|shannon|max_read_prop")) %>%
  mutate(model = map(.x = value, .f = ~ get(.x))) %>%
  mutate(fit = map(.x = model, .f = ~ .x$fit)) %>%
  mutate(param = gsub("m_minimpute_svp_minimpute_","b_",gsub("_brms","",value))) %>%
  mutate(posterior = map(.x = fit, .f = ~ posterior::as_draws_df(.x) %>% as_tibble() %>% select(contains("b_")) %>% pull())) %>%
  select(value, parameter, posterior) %>%
  unnest(cols = c(posterior)) %>%
  mutate(posterior_OR = exp(posterior)) %>%
  identity() %>%
  # choose primary outcome definition
  filter(grepl("m_minimpute_svp_minimpute_",value)) %>%
  #
  ggplot(data = ., aes(y = parameter, x = posterior, fill = stat(x > 0))) +
  tidybayes::stat_halfeye(alpha = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("gray80","skyblue"))



#' odds ratio scale
tibble(value = c("m_minimpute_sv_shannon_low_count_brms",
                 "m_minimpute_svp_minimpute_shannon_low_count_brms",
                 "m_minimpute_svaps_shannon_low_count_brms",
                 "m_minimpute_sv_max_read_prop_minimpute_high_count_brms",
                 "m_minimpute_svp_minimpute_max_read_prop_minimpute_high_count_brms",
                 "m_minimpute_svaps_max_read_prop_minimpute_high_count_brms",
                 "m_minimpute_sv_log_copy_16S_per_ml_sputum_high_count_brms",
                 "m_minimpute_svp_minimpute_log_copy_16S_per_ml_sputum_high_count_brms",
                 "m_minimpute_svaps_log_copy_16S_per_ml_sputum_high_count_brms")) %>%
  mutate(parameter = stringr::str_extract(string = value, pattern = "log_copy_16S_per_ml_sputum|shannon|max_read_prop")) %>%
  mutate(model = map(.x = value, .f = ~ get(.x))) %>%
  mutate(fit = map(.x = model, .f = ~ .x$fit)) %>%
  mutate(param = gsub("m_minimpute_svp_minimpute_","b_",gsub("_brms","",value))) %>%
  mutate(posterior = map(.x = fit, .f = ~ posterior::as_draws_df(.x) %>% as_tibble() %>% select(contains("b_")) %>% pull())) %>%
  select(value, parameter, posterior) %>%
  unnest(cols = c(posterior)) %>%
  mutate(posterior_OR = exp(posterior)) %>%
  identity() %>%
  # choose primary outcome definition
  filter(grepl("m_minimpute_svp_minimpute_",value)) %>%
  #
  ggplot(data = ., aes(y = parameter, x = posterior_OR, fill = stat(x > 1))) +
  tidybayes::stat_halfeye(alpha = 0.8) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_fill_manual(values = c("gray80","skyblue"))



#' format plots

outcome_names <- c("Suspected VA-LRTI",
                   "VA-LRTI with Clinical Signs<br>& Positive Culture",
                   "VA-LRTI with Clinical Signs<br>& Positive Culture<br>& Antibiotic Change")

exposure_names <- c("Low Shannon Diversity",
                    "Single Dominant ASV by<br>Proportional Abundance",
                    "High Total<br>Bacterial Abundance<br>(16S rRNA qPCR)")






#' odds ratio scale
tibble(value = c("m_minimpute_sv_shannon_low_count_brms",
                 "m_minimpute_svp_minimpute_shannon_low_count_brms",
                 "m_minimpute_svaps_shannon_low_count_brms",
                 "m_minimpute_sv_max_read_prop_minimpute_high_count_brms",
                 "m_minimpute_svp_minimpute_max_read_prop_minimpute_high_count_brms",
                 "m_minimpute_svaps_max_read_prop_minimpute_high_count_brms",
                 "m_minimpute_sv_log_copy_16S_per_ml_sputum_high_count_brms",
                 "m_minimpute_svp_minimpute_log_copy_16S_per_ml_sputum_high_count_brms",
                 "m_minimpute_svaps_log_copy_16S_per_ml_sputum_high_count_brms")) %>%
  mutate(parameter = stringr::str_extract(string = value, pattern = "log_copy_16S_per_ml_sputum|shannon|max_read_prop")) %>%
  mutate(model = map(.x = value, .f = ~ get(.x))) %>%
  mutate(fit = map(.x = model, .f = ~ .x$fit)) %>%
  mutate(param = gsub("m_minimpute_svp_minimpute_","b_",gsub("_brms","",value))) %>%
  mutate(posterior = map(.x = fit, .f = ~ posterior::as_draws_df(.x) %>% as_tibble() %>% select(contains("b_")) %>% pull())) %>%
  select(value, parameter, posterior) %>%
  unnest(cols = c(posterior)) %>%
  mutate(posterior_OR = exp(posterior)) %>%
  #group_by(parameter) %>%
  #mutate(posterior_OR_median = c(median(posterior_OR, na.rm = TRUE),rep(NA, n() - 1))) %>%
  #ungroup() %>%
  mutate(parameter = case_when(parameter == "shannon" ~ exposure_names[1],
                               parameter == "max_read_prop" ~ exposure_names[2],
                               parameter == "log_copy_16S_per_ml_sputum" ~ exposure_names[3])) %>%
  mutate(parameter = factor(parameter, level = exposure_names)) %>%
  identity() %>%
  # choose primary outcome definition
  filter(grepl("m_minimpute_",value)) %>%
  #
  ggplot(data = ., aes(y = parameter, x = posterior_OR)) +
  #geom_point(aes(y = parameter, x = posterior_OR_median)) +
  tidybayes::stat_interval(.width = c(0.5,0.8,0.95)) +
  tidybayes::stat_pointinterval(.width = c(0)) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  colorspace::scale_color_discrete_sequential(palette = "Blues 3", nmax = 5, order = 2:4) +
  theme_bw() +
  theme(axis.text.x = ggtext::element_markdown(color = "black"),
        axis.text.y = ggtext::element_markdown(color = "black"),
        axis.title.x = ggtext::element_markdown(color = "black"),
        legend.title = ggtext::element_markdown(color = "black"),
        legend.background = element_rect(fill = "white", color = "black", size = 0.25),
        legend.position = c(0.86,0.78)) +
  labs(x = "Odds Ratio of VA-LRTI per Day<br>of Persistent Microbiome Disruption",
       y = "",
       color = "Posterior<br>Credible<br>Interval") -> p_minimpute_persistent_mdi
p_minimpute_persistent_mdi



# p_minimpute_persistent_mdi %>%
#   ggsave(filename = "./figs/p_minimpute_persistent_mdi.pdf", height = 4, width = 5, units = "in")
# p_minimpute_persistent_mdi %>%
#   ggsave(filename = "./figs/p_minimpute_persistent_mdi.svg", height = 4, width = 5, units = "in")
# p_minimpute_persistent_mdi %>%
#   ggsave(filename = "./figs/p_minimpute_persistent_mdi.png", height = 4, width = 5, units = "in", dpi = 600)


tibble(value = c("m_minimpute_sv_shannon_low_count_brms",
                 "m_minimpute_svp_minimpute_shannon_low_count_brms",
                 "m_minimpute_svaps_shannon_low_count_brms",
                 "m_minimpute_sv_max_read_prop_minimpute_high_count_brms",
                 "m_minimpute_svp_minimpute_max_read_prop_minimpute_high_count_brms",
                 "m_minimpute_svaps_max_read_prop_minimpute_high_count_brms",
                 "m_minimpute_sv_log_copy_16S_per_ml_sputum_high_count_brms",
                 "m_minimpute_svp_minimpute_log_copy_16S_per_ml_sputum_high_count_brms",
                 "m_minimpute_svaps_log_copy_16S_per_ml_sputum_high_count_brms")) %>%
  mutate(parameter = stringr::str_extract(string = value, pattern = "log_copy_16S_per_ml_sputum|shannon|max_read_prop")) %>%
  mutate(model = map(.x = value, .f = ~ get(.x))) %>%
  mutate(fit = map(.x = model, .f = ~ .x$fit)) %>%
  mutate(param = gsub("m_minimpute_svp_minimpute_","b_",gsub("_brms","",value))) %>%
  mutate(posterior = map(.x = fit, .f = ~ posterior::as_draws_df(.x) %>% as_tibble() %>% select(contains("b_")) %>% pull())) %>%
  select(value, parameter, posterior) %>%
  unnest(cols = c(posterior)) %>%
  mutate(posterior_OR = exp(posterior)) %>%
  #group_by(parameter) %>%
  #mutate(posterior_OR_median = c(median(posterior_OR, na.rm = TRUE),rep(NA, n() - 1))) %>%
  #ungroup() %>%
  mutate(parameter = case_when(parameter == "shannon" ~ exposure_names[1],
                               parameter == "max_read_prop" ~ exposure_names[2],
                               parameter == "log_copy_16S_per_ml_sputum" ~ exposure_names[3])) %>%
  mutate(parameter = factor(parameter, level = exposure_names)) %>%
  identity() %>%
  # choose primary outcome definition
  filter(grepl("m_minimpute_",value)) %>%
  group_by(parameter) %>%
  tidybayes::median_qi(posterior_OR) %>%
  ungroup() %>%
  gt() %>%
  gt::fmt_number(columns = c(posterior_OR, .lower, .upper), n_sigfig = 3) %>%
  gt::tab_header(title = "Minimal Imputation")



#' #################################################
#' #################################################
#' 
#' COMPARE MODEL PARAMETERS (PROBABILITY & ODDS RATIO SCALES)
#' ACROSS OUTCOME DEFINITIONS
#' 
#' #################################################
#' #################################################

#' format plots

outcome_names <- c("Suspected VA-LRTI",
                   "VA-LRTI with Clinical Signs<br>& Positive Culture",
                   "VA-LRTI with Clinical Signs<br>& Positive Culture<br>& Antibiotic Change")

exposure_names <- c("Low Shannon Diversity",
                    "Single Dominant ASV by<br>Proportional Abundance",
                    "High Total<br>Bacterial Abundance<br>(16S rRNA qPCR)")



#' odds ratio scale
tibble(value = c("m_minimpute_sv_shannon_low_count_brms",
                 "m_minimpute_svp_minimpute_shannon_low_count_brms",
                 "m_minimpute_svaps_shannon_low_count_brms",
                 "m_minimpute_sv_max_read_prop_minimpute_high_count_brms",
                 "m_minimpute_svp_minimpute_max_read_prop_minimpute_high_count_brms",
                 "m_minimpute_svaps_max_read_prop_minimpute_high_count_brms",
                 "m_minimpute_sv_log_copy_16S_per_ml_sputum_high_count_brms",
                 "m_minimpute_svp_minimpute_log_copy_16S_per_ml_sputum_high_count_brms",
                 "m_minimpute_svaps_log_copy_16S_per_ml_sputum_high_count_brms")) %>%
  mutate(parameter = stringr::str_extract(string = value, pattern = "log_copy_16S_per_ml_sputum|shannon|max_read_prop")) %>%
  mutate(model = map(.x = value, .f = ~ get(.x))) %>%
  mutate(fit = map(.x = model, .f = ~ .x$fit)) %>%
  mutate(param = gsub("m_","b_",gsub("_brms","",gsub("sv_|svp_minimpute_|svap_|svaps_","",value)))) %>%
  mutate(posterior = map(.x = fit, .f = ~ posterior::as_draws_df(.x) %>% as_tibble() %>% select(contains("b_")) %>% pull())) %>%
  mutate(outcome_definition = stringr::str_extract(string = value, pattern = "sv_|svp_minimpute_|svap_|svaps_"),
         outcome_definition = case_when(outcome_definition == "sv_" ~ outcome_names[1],
                                        outcome_definition == "svp_minimpute_" ~ outcome_names[2],
                                        outcome_definition == "svaps_" ~ outcome_names[3])) %>%
  mutate(parameter = case_when(parameter == "shannon" ~ exposure_names[1],
                               parameter == "max_read_prop" ~ exposure_names[2],
                               parameter == "log_copy_16S_per_ml_sputum" ~ exposure_names[3])) %>%
  mutate(parameter = factor(parameter, level = exposure_names)) %>%
  select(value, parameter, posterior, outcome_definition) %>%
  unnest(cols = c(posterior)) %>%
  mutate(posterior_OR = exp(posterior)) %>%
  ggplot(data = ., aes(y = parameter, x = posterior_OR, fill = stat(x > 1))) +
  #tidybayes::stat_halfeye(alpha = 0.9, height = 0.5) +
  tidybayes::stat_halfeye(alpha = 0.9, height = 0.5) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_fill_manual(values = c("gray80","skyblue")) +
  facet_wrap(facets = ~ outcome_definition, scales = "fixed") +
  theme_bw() +
  theme(legend.position = "none",
        panel.spacing = unit(2, "lines"),
        axis.text.x = ggtext::element_markdown(color = "black"),
        axis.text.y = ggtext::element_markdown(color = "black"),
        axis.title.x = ggtext::element_markdown(color = "black"),
        axis.title.y = ggtext::element_markdown(color = "black"),
        strip.text.x = ggtext::element_markdown(color = "black"),
        strip.text.y = ggtext::element_markdown(color = "black"),
        strip.background = element_blank(),
        ) +
  labs(x = "Odds Ratio of VA-LRTI per Day of Persistent Microbiome Disruption",
       y = "") -> p_minimpute_compare_outcome_definitions_persistence
p_minimpute_compare_outcome_definitions_persistence


# p_minimpute_compare_outcome_definitions_persistence %>%
#   ggsave(filename = "./figs/p_minimpute_compare_outcome_definitions_persistence.pdf", height = 3, width = 7, units = "in")
# p_minimpute_compare_outcome_definitions_persistence %>%
#   ggsave(filename = "./figs/p_minimpute_compare_outcome_definitions_persistence.svg", height = 3, width = 7, units = "in")
# p_minimpute_compare_outcome_definitions_persistence %>%
#   ggsave(filename = "./figs/p_minimpute_compare_outcome_definitions_persistence.png", height = 3, width = 7, units = "in", dpi = 600)



