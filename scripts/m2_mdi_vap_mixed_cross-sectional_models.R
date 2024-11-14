#' ###########################################
#' MIXED EFFECT VAP ~ CROSS-SECTIONAL MDI models
#' ###########################################



#' load libraries and set seed / control variables
library(tidyverse)
library(tidybayes)
library(brms)
set.seed(16)
run_models = FALSE
save_models = FALSE
run_plots = FALSE




#' #######################################################################
#' #######################################################################
#' VAP as a function of lower respiratory tract diversity: SHANNON
#' #######################################################################
#' #######################################################################


#' #######################################################################
#' SUSPECTED VAP
#' #######################################################################

dat %>%
  # filter to include dates only up to VAP (or last date, if no VAP)
  group_by(subject_id) %>%
  mutate(first_suspected_vap = ifelse(is.na(first_suspected_vap), max(date), first_suspected_vap)) %>%
  filter(date <= first_suspected_vap) %>%
  ungroup() %>%
  select(subject_id, subject_day, shannon, max_read_prop, log_copy_16S_per_ml_sputum, suspected_vap) %>%
  mutate_at(.vars = vars(-subject_id, -subject_day, -suspected_vap), .funs = list("scaled" = ~ as.vector(scale(.x)))) %>%
  distinct() %>%
  filter(complete.cases(.)) %>%
  identity() -> d_vap

d_vap


if(run_models) {
d_vap %>%
  get_prior(data = ., family = bernoulli,
            suspected_vap ~ (1 | subject_id) + shannon_scaled)


d_vap %>%
  brm(data = ., family = bernoulli,
      suspected_vap ~ (1 | subject_id) + shannon_scaled,
      prior = c(prior(student_t(3,0,2.5), class = Intercept),
                prior(student_t(3,0,2.5), class = b)
      ),
      sample_prior = "only", # prior predictive
      iter = 2000,
      warmup = 1000,
      chains = 4,
      cores = 4,
      control = list("adapt_delta" = 0.99, max_treedepth = 16),
      backend = "cmdstanr",
      seed = 16) -> d_prior


d_prior %>%
  posterior_epred() %>%
  as_tibble() %>%
  gather(key = "observations", value = "pp_samples") %>%
  qplot(data = ., x = pp_samples, geom = "density")



#' brms
d_vap %>%
  brm(data = ., family = bernoulli,
      suspected_vap ~ (1 | subject_id) + shannon_scaled,
      prior = c(prior(student_t(3,0,2.5), class = Intercept),
                prior(student_t(3,0,2.5), class = b)
      ),
      #sample_prior = "only", # prior predictive
      iter = 2000,
      warmup = 1000,
      chains = 4,
      cores = 4,
      control = list("adapt_delta" = 0.99, max_treedepth = 16),
      backend = "cmdstanr",
      seed = 16) -> m_suspected_vap_shannon_brms
}

if(save_models) {
m_suspected_vap_shannon_brms %>% write_rds(file = "./models/m_suspected_vap_shannon_brms_mixed.rds.gz", compress = "gz")
m_suspected_vap_shannon_brms$fit %>% write_rds(file = "./models/m_suspected_vap_shannon_brms_stanfit_mixed.rds.gz", compress = "gz")
}

m_suspected_vap_shannon_brms <- read_rds(file = "./models/m_suspected_vap_shannon_brms_mixed.rds.gz")

m_suspected_vap_shannon_brms$formula

pp_check(m_suspected_vap_shannon_brms)
rstan::check_hmc_diagnostics(m_suspected_vap_shannon_brms$fit)


#' model results
m_suspected_vap_shannon_brms %>%
  brms::posterior_summary() %>%
  as_tibble(rownames = "param") %>%
  mutate_if(is.numeric, list("OR" = ~ exp(.x))) %>%
  identity() -> m_suspected_vap_shannon_brms_param_est
m_suspected_vap_shannon_brms_param_est

if(run_plots) {
m_suspected_vap_shannon_brms_param_est %>%
  filter(grepl("b_",param)) %>%
  ggplot(data = .) + 
  geom_segment(aes(x = Q2.5_OR, xend = Q97.5_OR, y = param, yend = param)) +
  geom_point(aes(x = Estimate_OR, y = param)) +
  geom_vline(xintercept = 1, linetype = 2)
}


#' mixed model plot: by subject
if(run_plots) {
  m_suspected_vap_shannon_brms$data %>%
    as_tibble() %>%
    expand(shannon_scaled = modelr::seq_range(shannon_scaled, n = 100), subject_id = unique(subject_id)) %>%
    mutate(shannon = shannon_scaled * sd(d_vap$shannon, na.rm = TRUE) + mean(d_vap$shannon, na.rm = TRUE)) %>%
    tidybayes::add_fitted_draws(newdata = ., model = m_suspected_vap_shannon_brms) %>%
    ungroup() %>%
    ggplot(data = ., aes(x = shannon_scaled, y = .value)) +
    tidybayes::stat_lineribbon(.width = c(0.5, 0.8, 0.95), alpha = 0.7, color = colorspace::sequential_hcl(5, palette = "Blues 3")[1]) +
    colorspace::scale_fill_discrete_sequential(palette = "Blues 3", nmax = 5, order = 2:4) +
    facet_wrap(~ subject_id)
}


#' mixed model plot: marginal
if(run_plots) {
  m_suspected_vap_shannon_brms$data %>%
    as_tibble() %>%
    expand(shannon_scaled = modelr::seq_range(shannon_scaled, n = 100), subject_id = unique(subject_id)) %>%
    mutate(shannon = shannon_scaled * sd(d_vap$shannon, na.rm = TRUE) + mean(d_vap$shannon, na.rm = TRUE)) %>%
    tidybayes::add_fitted_draws(newdata = ., model = m_suspected_vap_shannon_brms, re_formula = NULL, seed = 16) %>%
    ungroup() %>%
    ggplot(data = ., aes(x = shannon, y = .value)) +
    tidybayes::stat_lineribbon(.width = c(0.5, 0.8, 0.95), alpha = 0.7, color = colorspace::sequential_hcl(5, palette = "Blues 3")[1]) +
    colorspace::scale_fill_discrete_sequential(palette = "Blues 3", nmax = 5, order = 2:4)
}







#' #######################################################################
#' SUSPECTED VAP + CX + SIGNS (same as SUSPECTED VAP + CX)
#' #######################################################################

dat %>%
  # filter to include dates only up to VAP (or last date, if no VAP)
  group_by(subject_id) %>%
  mutate(first_suspected_vap_plus_pathogen = ifelse(is.na(first_suspected_vap_plus_pathogen), max(date), first_suspected_vap_plus_pathogen)) %>%
  filter(date <= first_suspected_vap_plus_pathogen) %>%
  ungroup() %>%
  select(subject_id, subject_day, shannon, max_read_prop, log_copy_16S_per_ml_sputum, suspected_vap_plus_signs_and_pathogen) %>%
  mutate_at(.vars = vars(-subject_id, -subject_day, -suspected_vap_plus_signs_and_pathogen), .funs = list("scaled" = ~ as.vector(scale(.x)))) %>%
  distinct() %>%
  filter(complete.cases(.)) %>%
  identity() -> d_vap

d_vap


if (run_models) {
d_vap %>%
  get_prior(data = ., family = bernoulli,
            suspected_vap_plus_signs_and_pathogen ~ (1 | subject_id) + shannon_scaled)


d_vap %>%
  brm(data = ., family = bernoulli,
      suspected_vap_plus_signs_and_pathogen ~ (1 | subject_id) + shannon_scaled,
      prior = c(prior(student_t(3,0,2.5), class = Intercept),
        prior(student_t(3,0,2.5), class = b)
      ),
      sample_prior = "only", # prior predictive
      iter = 2000,
      warmup = 1000,
      chains = 4,
      cores = 4,
      control = list("adapt_delta" = 0.99, max_treedepth = 16),
      backend = "cmdstanr",
      seed = 16) -> d_prior


d_prior %>%
  posterior_epred() %>%
  as_tibble() %>%
  gather(key = "observations", value = "pp_samples") %>%
  qplot(data = ., x = pp_samples, geom = "density")



#' brms
d_vap %>%
  brm(data = ., family = bernoulli,
      suspected_vap_plus_signs_and_pathogen ~ (1 | subject_id) + shannon_scaled,
      prior = c(prior(student_t(3,0,2.5), class = Intercept),
        prior(student_t(3,0,2.5), class = b)
      ),
      #sample_prior = "only", # prior predictive
      iter = 2000,
      warmup = 1000,
      chains = 4,
      cores = 4,
      control = list("adapt_delta" = 0.99, max_treedepth = 16),
      backend = "cmdstanr",
      seed = 16) -> m_suspected_vap_plus_signs_and_pathogen_shannon_brms
}

if(save_models) {
m_suspected_vap_plus_signs_and_pathogen_shannon_brms %>% write_rds(file = "./models/m_suspected_vap_plus_signs_and_pathogen_shannon_brms_mixed.rds.gz", compress = "gz")
m_suspected_vap_plus_signs_and_pathogen_shannon_brms$fit %>% write_rds(file = "./models/m_suspected_vap_plus_signs_and_pathogen_shannon_brms_stanfit_mixed.rds.gz", compress = "gz")
}

m_suspected_vap_plus_signs_and_pathogen_shannon_brms <- read_rds(file = "./models/m_suspected_vap_plus_signs_and_pathogen_shannon_brms_mixed.rds.gz")

m_suspected_vap_plus_signs_and_pathogen_shannon_brms$formula

pp_check(m_suspected_vap_plus_signs_and_pathogen_shannon_brms)
rstan::check_hmc_diagnostics(m_suspected_vap_plus_signs_and_pathogen_shannon_brms$fit)


#' model results
m_suspected_vap_plus_signs_and_pathogen_shannon_brms %>%
  brms::posterior_summary() %>%
  as_tibble(rownames = "param") %>%
  mutate_if(is.numeric, list("OR" = ~ exp(.x))) %>%
  identity() -> m_suspected_vap_plus_signs_and_pathogen_shannon_brms_param_est
m_suspected_vap_plus_signs_and_pathogen_shannon_brms_param_est

if(run_plots) {
  m_suspected_vap_plus_signs_and_pathogen_shannon_brms_param_est %>%
    filter(grepl("b_",param)) %>%
    ggplot(data = .) + 
    geom_segment(aes(x = Q2.5_OR, xend = Q97.5_OR, y = param, yend = param)) +
    geom_point(aes(x = Estimate_OR, y = param)) +
    geom_vline(xintercept = 1, linetype = 2)
}


#' mixed model plot: by subject
if(run_plots) {
  m_suspected_vap_plus_signs_and_pathogen_shannon_brms$data %>%
    as_tibble() %>%
    expand(shannon_scaled = modelr::seq_range(shannon_scaled, n = 100), subject_id = unique(subject_id)) %>%
    mutate(shannon = shannon_scaled * sd(d_vap$shannon, na.rm = TRUE) + mean(d_vap$shannon, na.rm = TRUE)) %>%
    tidybayes::add_fitted_draws(newdata = ., model = m_suspected_vap_plus_signs_and_pathogen_shannon_brms) %>%
    ungroup() %>%
    ggplot(data = ., aes(x = shannon_scaled, y = .value)) +
    tidybayes::stat_lineribbon(.width = c(0.5, 0.8, 0.95), alpha = 0.7, color = colorspace::sequential_hcl(5, palette = "Blues 3")[1]) +
    colorspace::scale_fill_discrete_sequential(palette = "Blues 3", nmax = 5, order = 2:4) +
    facet_wrap(~ subject_id)
}


#' mixed model plot: marginal
if(run_plots) {
  m_suspected_vap_plus_signs_and_pathogen_shannon_brms$data %>%
    as_tibble() %>%
    expand(shannon_scaled = modelr::seq_range(shannon_scaled, n = 100), subject_id = unique(subject_id)) %>%
    mutate(shannon = shannon_scaled * sd(d_vap$shannon, na.rm = TRUE) + mean(d_vap$shannon, na.rm = TRUE)) %>%
    tidybayes::add_fitted_draws(newdata = ., model = m_suspected_vap_plus_signs_and_pathogen_shannon_brms, re_formula = NULL, seed = 16) %>%
    ungroup() %>%
    ggplot(data = ., aes(x = shannon, y = .value)) +
    tidybayes::stat_lineribbon(.width = c(0.5, 0.8, 0.95), alpha = 0.7, color = colorspace::sequential_hcl(5, palette = "Blues 3")[1]) +
    colorspace::scale_fill_discrete_sequential(palette = "Blues 3", nmax = 5, order = 2:4)
}




#' #######################################################################
#' SUSPECTED VAP + CX + SIGNS + ABX
#' #######################################################################

dat %>%
  # filter to include dates only up to VAP (or last date, if no VAP)
  group_by(subject_id) %>%
  mutate(first_suspected_vap_plus_abx_and_pathogen_and_signs = ifelse(is.na(first_suspected_vap_plus_abx_and_pathogen_and_signs), max(date), first_suspected_vap_plus_abx_and_pathogen_and_signs)) %>%
  filter(date <= first_suspected_vap_plus_abx_and_pathogen_and_signs) %>%
  ungroup() %>%
  select(subject_id, subject_day, shannon, max_read_prop, log_copy_16S_per_ml_sputum, suspected_vap_plus_abx_and_pathogen_and_signs) %>%
  mutate_at(.vars = vars(-subject_id, -subject_day, -suspected_vap_plus_abx_and_pathogen_and_signs), .funs = list("scaled" = ~ as.vector(scale(.x)))) %>%
  distinct() %>%
  filter(complete.cases(.)) %>%
  identity() -> d_vap

d_vap


if (run_models) {
  d_vap %>%
    get_prior(data = ., family = bernoulli,
              suspected_vap_plus_abx_and_pathogen_and_signs ~ (1 | subject_id) + shannon_scaled)
  
  
  d_vap %>%
    brm(data = ., family = bernoulli,
        suspected_vap_plus_abx_and_pathogen_and_signs ~ (1 | subject_id) + shannon_scaled,
        prior = c(prior(student_t(3,0,2.5), class = Intercept),
                  prior(student_t(3,0,2.5), class = b)
        ),
        sample_prior = "only", # prior predictive
        iter = 2000,
        warmup = 1000,
        chains = 4,
        cores = 4,
        control = list("adapt_delta" = 0.99, max_treedepth = 16),
        backend = "cmdstanr",
        seed = 16) -> d_prior
  
  
  d_prior %>%
    posterior_epred() %>%
    as_tibble() %>%
    gather(key = "observations", value = "pp_samples") %>%
    qplot(data = ., x = pp_samples, geom = "density")
  
  
  
  #' brms
  d_vap %>%
    brm(data = ., family = bernoulli,
        suspected_vap_plus_abx_and_pathogen_and_signs ~ (1 | subject_id) + shannon_scaled,
        prior = c(prior(student_t(3,0,2.5), class = Intercept),
                  prior(student_t(3,0,2.5), class = b)
        ),
        #sample_prior = "only", # prior predictive
        iter = 2000,
        warmup = 1000,
        chains = 4,
        cores = 4,
        control = list("adapt_delta" = 0.99, max_treedepth = 16),
        backend = "cmdstanr",
        seed = 16) -> m_suspected_vap_plus_abx_and_pathogen_and_signs_shannon_brms
}

if(save_models) {
  m_suspected_vap_plus_abx_and_pathogen_and_signs_shannon_brms %>% write_rds(file = "./models/m_suspected_vap_plus_abx_and_pathogen_and_signs_shannon_brms_mixed.rds.gz", compress = "gz")
  m_suspected_vap_plus_abx_and_pathogen_and_signs_shannon_brms$fit %>% write_rds(file = "./models/m_suspected_vap_plus_abx_and_pathogen_and_signs_shannon_brms_stanfit_mixed.rds.gz", compress = "gz")
}

m_suspected_vap_plus_abx_and_pathogen_and_signs_shannon_brms <- read_rds(file = "./models/m_suspected_vap_plus_abx_and_pathogen_and_signs_shannon_brms_mixed.rds.gz")

m_suspected_vap_plus_abx_and_pathogen_and_signs_shannon_brms$formula

pp_check(m_suspected_vap_plus_abx_and_pathogen_and_signs_shannon_brms)
rstan::check_hmc_diagnostics(m_suspected_vap_plus_abx_and_pathogen_and_signs_shannon_brms$fit)


#' model results
m_suspected_vap_plus_abx_and_pathogen_and_signs_shannon_brms %>%
  brms::posterior_summary() %>%
  as_tibble(rownames = "param") %>%
  mutate_if(is.numeric, list("OR" = ~ exp(.x))) %>%
  identity() -> m_suspected_vap_plus_abx_and_pathogen_and_signs_shannon_brms_param_est
m_suspected_vap_plus_abx_and_pathogen_and_signs_shannon_brms_param_est

if(run_plots) {
  m_suspected_vap_plus_abx_and_pathogen_and_signs_shannon_brms_param_est %>%
    filter(grepl("b_",param)) %>%
    ggplot(data = .) + 
    geom_segment(aes(x = Q2.5_OR, xend = Q97.5_OR, y = param, yend = param)) +
    geom_point(aes(x = Estimate_OR, y = param)) +
    geom_vline(xintercept = 1, linetype = 2)
}


#' mixed model plot: by subject
if(run_plots) {
  m_suspected_vap_plus_abx_and_pathogen_and_signs_shannon_brms$data %>%
    as_tibble() %>%
    expand(shannon_scaled = modelr::seq_range(shannon_scaled, n = 100), subject_id = unique(subject_id)) %>%
    mutate(shannon = shannon_scaled * sd(d_vap$shannon, na.rm = TRUE) + mean(d_vap$shannon, na.rm = TRUE)) %>%
    tidybayes::add_fitted_draws(newdata = ., model = m_suspected_vap_plus_abx_and_pathogen_and_signs_shannon_brms) %>%
    ungroup() %>%
    ggplot(data = ., aes(x = shannon_scaled, y = .value)) +
    tidybayes::stat_lineribbon(.width = c(0.5, 0.8, 0.95), alpha = 0.7, color = colorspace::sequential_hcl(5, palette = "Blues 3")[1]) +
    colorspace::scale_fill_discrete_sequential(palette = "Blues 3", nmax = 5, order = 2:4) +
    facet_wrap(~ subject_id)
}


#' mixed model plot: marginal
if(run_plots) {
  m_suspected_vap_plus_abx_and_pathogen_and_signs_shannon_brms$data %>%
    as_tibble() %>%
    expand(shannon_scaled = modelr::seq_range(shannon_scaled, n = 100), subject_id = unique(subject_id)) %>%
    mutate(shannon = shannon_scaled * sd(d_vap$shannon, na.rm = TRUE) + mean(d_vap$shannon, na.rm = TRUE)) %>%
    tidybayes::add_fitted_draws(newdata = ., model = m_suspected_vap_plus_abx_and_pathogen_and_signs_shannon_brms, re_formula = NULL, seed = 16) %>%
    ungroup() %>%
    ggplot(data = ., aes(x = shannon, y = .value)) +
    tidybayes::stat_lineribbon(.width = c(0.5, 0.8, 0.95), alpha = 0.7, color = colorspace::sequential_hcl(5, palette = "Blues 3")[1]) +
    colorspace::scale_fill_discrete_sequential(palette = "Blues 3", nmax = 5, order = 2:4)
}











#' #######################################################################
#' #######################################################################
#' VAP as a function of dominant ASV proportional abundance: MAX_READ_PROP
#' #######################################################################
#' #######################################################################


#' #######################################################################
#' SUSPECTED VAP
#' #######################################################################

dat %>%
  # filter to include dates only up to VAP (or last date, if no VAP)
  group_by(subject_id) %>%
  mutate(first_suspected_vap = ifelse(is.na(first_suspected_vap), max(date), first_suspected_vap)) %>%
  filter(date <= first_suspected_vap) %>%
  ungroup() %>%
  select(subject_id, subject_day, max_read_prop, max_read_prop, log_copy_16S_per_ml_sputum, suspected_vap) %>%
  mutate_at(.vars = vars(-subject_id, -subject_day, -suspected_vap), .funs = list("scaled" = ~ as.vector(scale(.x)))) %>%
  distinct() %>%
  filter(complete.cases(.)) %>%
  identity() -> d_vap

d_vap


if(run_models) {
  d_vap %>%
    get_prior(data = ., family = bernoulli,
              suspected_vap ~ (1 | subject_id) + max_read_prop_scaled)
  
  
  d_vap %>%
    brm(data = ., family = bernoulli,
        suspected_vap ~ (1 | subject_id) + max_read_prop_scaled,
        prior = c(prior(student_t(3,0,2.5), class = Intercept),
                  prior(student_t(3,0,2.5), class = b)
        ),
        sample_prior = "only", # prior predictive
        iter = 2000,
        warmup = 1000,
        chains = 4,
        cores = 4,
        control = list("adapt_delta" = 0.99, max_treedepth = 16),
        backend = "cmdstanr",
        seed = 16) -> d_prior
  
  
  d_prior %>%
    posterior_epred() %>%
    as_tibble() %>%
    gather(key = "observations", value = "pp_samples") %>%
    qplot(data = ., x = pp_samples, geom = "density")
  
  
  
  #' brms
  d_vap %>%
    brm(data = ., family = bernoulli,
        suspected_vap ~ (1 | subject_id) + max_read_prop_scaled,
        prior = c(prior(student_t(3,0,2.5), class = Intercept),
                  prior(student_t(3,0,2.5), class = b)
        ),
        #sample_prior = "only", # prior predictive
        iter = 2000,
        warmup = 1000,
        chains = 4,
        cores = 4,
        control = list("adapt_delta" = 0.99, max_treedepth = 16),
        backend = "cmdstanr",
        seed = 16) -> m_suspected_vap_max_read_prop_brms
}

if(save_models) {
  m_suspected_vap_max_read_prop_brms %>% write_rds(file = "./models/m_suspected_vap_max_read_prop_brms_mixed.rds.gz", compress = "gz")
  m_suspected_vap_max_read_prop_brms$fit %>% write_rds(file = "./models/m_suspected_vap_max_read_prop_brms_stanfit_mixed.rds.gz", compress = "gz")
}

m_suspected_vap_max_read_prop_brms <- read_rds(file = "./models/m_suspected_vap_max_read_prop_brms_mixed.rds.gz")

m_suspected_vap_max_read_prop_brms$formula

pp_check(m_suspected_vap_max_read_prop_brms)
rstan::check_hmc_diagnostics(m_suspected_vap_max_read_prop_brms$fit)


#' model results
m_suspected_vap_max_read_prop_brms %>%
  brms::posterior_summary() %>%
  as_tibble(rownames = "param") %>%
  mutate_if(is.numeric, list("OR" = ~ exp(.x))) %>%
  identity() -> m_suspected_vap_max_read_prop_brms_param_est
m_suspected_vap_max_read_prop_brms_param_est

if(run_plots) {
  m_suspected_vap_max_read_prop_brms_param_est %>%
    filter(grepl("b_",param)) %>%
    ggplot(data = .) + 
    geom_segment(aes(x = Q2.5_OR, xend = Q97.5_OR, y = param, yend = param)) +
    geom_point(aes(x = Estimate_OR, y = param)) +
    geom_vline(xintercept = 1, linetype = 2)
}


#' mixed model plot: by subject
if(run_plots) {
  m_suspected_vap_max_read_prop_brms$data %>%
    as_tibble() %>%
    expand(max_read_prop_scaled = modelr::seq_range(max_read_prop_scaled, n = 100), subject_id = unique(subject_id)) %>%
    mutate(max_read_prop = max_read_prop_scaled * sd(d_vap$max_read_prop, na.rm = TRUE) + mean(d_vap$max_read_prop, na.rm = TRUE)) %>%
    tidybayes::add_fitted_draws(newdata = ., model = m_suspected_vap_max_read_prop_brms) %>%
    ungroup() %>%
    ggplot(data = ., aes(x = max_read_prop_scaled, y = .value)) +
    tidybayes::stat_lineribbon(.width = c(0.5, 0.8, 0.95), alpha = 0.7, color = colorspace::sequential_hcl(5, palette = "Blues 3")[1]) +
    colorspace::scale_fill_discrete_sequential(palette = "Blues 3", nmax = 5, order = 2:4) +
    facet_wrap(~ subject_id)
}


#' mixed model plot: marginal
if(run_plots) {
  m_suspected_vap_max_read_prop_brms$data %>%
    as_tibble() %>%
    expand(max_read_prop_scaled = modelr::seq_range(max_read_prop_scaled, n = 100), subject_id = unique(subject_id)) %>%
    mutate(max_read_prop = max_read_prop_scaled * sd(d_vap$max_read_prop, na.rm = TRUE) + mean(d_vap$max_read_prop, na.rm = TRUE)) %>%
    tidybayes::add_fitted_draws(newdata = ., model = m_suspected_vap_max_read_prop_brms, re_formula = NULL, seed = 16) %>%
    ungroup() %>%
    ggplot(data = ., aes(x = max_read_prop, y = .value)) +
    tidybayes::stat_lineribbon(.width = c(0.5, 0.8, 0.95), alpha = 0.7, color = colorspace::sequential_hcl(5, palette = "Blues 3")[1]) +
    colorspace::scale_fill_discrete_sequential(palette = "Blues 3", nmax = 5, order = 2:4)
}







#' #######################################################################
#' SUSPECTED VAP + CX + SIGNS (same as SUSPECTED VAP + CX)
#' #######################################################################

dat %>%
  # filter to include dates only up to VAP (or last date, if no VAP)
  group_by(subject_id) %>%
  mutate(first_suspected_vap_plus_pathogen = ifelse(is.na(first_suspected_vap_plus_pathogen), max(date), first_suspected_vap_plus_pathogen)) %>%
  filter(date <= first_suspected_vap_plus_pathogen) %>%
  ungroup() %>%
  select(subject_id, subject_day, max_read_prop, max_read_prop, log_copy_16S_per_ml_sputum, suspected_vap_plus_signs_and_pathogen) %>%
  mutate_at(.vars = vars(-subject_id, -subject_day, -suspected_vap_plus_signs_and_pathogen), .funs = list("scaled" = ~ as.vector(scale(.x)))) %>%
  distinct() %>%
  filter(complete.cases(.)) %>%
  identity() -> d_vap

d_vap


if (run_models) {
  d_vap %>%
    get_prior(data = ., family = bernoulli,
              suspected_vap_plus_signs_and_pathogen ~ (1 | subject_id) + max_read_prop_scaled)
  
  
  d_vap %>%
    brm(data = ., family = bernoulli,
        suspected_vap_plus_signs_and_pathogen ~ (1 | subject_id) + max_read_prop_scaled,
        prior = c(prior(student_t(3,0,2.5), class = Intercept),
                  prior(student_t(3,0,2.5), class = b)
        ),
        sample_prior = "only", # prior predictive
        iter = 2000,
        warmup = 1000,
        chains = 4,
        cores = 4,
        control = list("adapt_delta" = 0.99, max_treedepth = 16),
        backend = "cmdstanr",
        seed = 16) -> d_prior
  
  
  d_prior %>%
    posterior_epred() %>%
    as_tibble() %>%
    gather(key = "observations", value = "pp_samples") %>%
    qplot(data = ., x = pp_samples, geom = "density")
  
  
  
  #' brms
  d_vap %>%
    brm(data = ., family = bernoulli,
        suspected_vap_plus_signs_and_pathogen ~ (1 | subject_id) + max_read_prop_scaled,
        prior = c(prior(student_t(3,0,2.5), class = Intercept),
                  prior(student_t(3,0,2.5), class = b)
        ),
        #sample_prior = "only", # prior predictive
        iter = 2000,
        warmup = 1000,
        chains = 4,
        cores = 4,
        control = list("adapt_delta" = 0.99, max_treedepth = 16),
        backend = "cmdstanr",
        seed = 16) -> m_suspected_vap_plus_signs_and_pathogen_max_read_prop_brms
}

if(save_models) {
  m_suspected_vap_plus_signs_and_pathogen_max_read_prop_brms %>% write_rds(file = "./models/m_suspected_vap_plus_signs_and_pathogen_max_read_prop_brms_mixed.rds.gz", compress = "gz")
  m_suspected_vap_plus_signs_and_pathogen_max_read_prop_brms$fit %>% write_rds(file = "./models/m_suspected_vap_plus_signs_and_pathogen_max_read_prop_brms_stanfit_mixed.rds.gz", compress = "gz")
}

m_suspected_vap_plus_signs_and_pathogen_max_read_prop_brms <- read_rds(file = "./models/m_suspected_vap_plus_signs_and_pathogen_max_read_prop_brms_mixed.rds.gz")

m_suspected_vap_plus_signs_and_pathogen_max_read_prop_brms$formula

pp_check(m_suspected_vap_plus_signs_and_pathogen_max_read_prop_brms)
rstan::check_hmc_diagnostics(m_suspected_vap_plus_signs_and_pathogen_max_read_prop_brms$fit)


#' model results
m_suspected_vap_plus_signs_and_pathogen_max_read_prop_brms %>%
  brms::posterior_summary() %>%
  as_tibble(rownames = "param") %>%
  mutate_if(is.numeric, list("OR" = ~ exp(.x))) %>%
  identity() -> m_suspected_vap_plus_signs_and_pathogen_max_read_prop_brms_param_est
m_suspected_vap_plus_signs_and_pathogen_max_read_prop_brms_param_est

if(run_plots) {
  m_suspected_vap_plus_signs_and_pathogen_max_read_prop_brms_param_est %>%
    filter(grepl("b_",param)) %>%
    ggplot(data = .) + 
    geom_segment(aes(x = Q2.5_OR, xend = Q97.5_OR, y = param, yend = param)) +
    geom_point(aes(x = Estimate_OR, y = param)) +
    geom_vline(xintercept = 1, linetype = 2)
}


#' mixed model plot: by subject
if(run_plots) {
  m_suspected_vap_plus_signs_and_pathogen_max_read_prop_brms$data %>%
    as_tibble() %>%
    expand(max_read_prop_scaled = modelr::seq_range(max_read_prop_scaled, n = 100), subject_id = unique(subject_id)) %>%
    mutate(max_read_prop = max_read_prop_scaled * sd(d_vap$max_read_prop, na.rm = TRUE) + mean(d_vap$max_read_prop, na.rm = TRUE)) %>%
    tidybayes::add_fitted_draws(newdata = ., model = m_suspected_vap_plus_signs_and_pathogen_max_read_prop_brms) %>%
    ungroup() %>%
    ggplot(data = ., aes(x = max_read_prop_scaled, y = .value)) +
    tidybayes::stat_lineribbon(.width = c(0.5, 0.8, 0.95), alpha = 0.7, color = colorspace::sequential_hcl(5, palette = "Blues 3")[1]) +
    colorspace::scale_fill_discrete_sequential(palette = "Blues 3", nmax = 5, order = 2:4) +
    facet_wrap(~ subject_id)
}


#' mixed model plot: marginal
if(run_plots) {
  m_suspected_vap_plus_signs_and_pathogen_max_read_prop_brms$data %>%
    as_tibble() %>%
    expand(max_read_prop_scaled = modelr::seq_range(max_read_prop_scaled, n = 100), subject_id = unique(subject_id)) %>%
    mutate(max_read_prop = max_read_prop_scaled * sd(d_vap$max_read_prop, na.rm = TRUE) + mean(d_vap$max_read_prop, na.rm = TRUE)) %>%
    tidybayes::add_fitted_draws(newdata = ., model = m_suspected_vap_plus_signs_and_pathogen_max_read_prop_brms, re_formula = NULL, seed = 16) %>%
    ungroup() %>%
    ggplot(data = ., aes(x = max_read_prop, y = .value)) +
    tidybayes::stat_lineribbon(.width = c(0.5, 0.8, 0.95), alpha = 0.7, color = colorspace::sequential_hcl(5, palette = "Blues 3")[1]) +
    colorspace::scale_fill_discrete_sequential(palette = "Blues 3", nmax = 5, order = 2:4)
}




#' #######################################################################
#' SUSPECTED VAP + CX + SIGNS + ABX
#' #######################################################################

dat %>%
  # filter to include dates only up to VAP (or last date, if no VAP)
  group_by(subject_id) %>%
  mutate(first_suspected_vap_plus_abx_and_pathogen_and_signs = ifelse(is.na(first_suspected_vap_plus_abx_and_pathogen_and_signs), max(date), first_suspected_vap_plus_abx_and_pathogen_and_signs)) %>%
  filter(date <= first_suspected_vap_plus_abx_and_pathogen_and_signs) %>%
  ungroup() %>%
  select(subject_id, subject_day, max_read_prop, max_read_prop, log_copy_16S_per_ml_sputum, suspected_vap_plus_abx_and_pathogen_and_signs) %>%
  mutate_at(.vars = vars(-subject_id, -subject_day, -suspected_vap_plus_abx_and_pathogen_and_signs), .funs = list("scaled" = ~ as.vector(scale(.x)))) %>%
  distinct() %>%
  filter(complete.cases(.)) %>%
  identity() -> d_vap

d_vap


if (run_models) {
  d_vap %>%
    get_prior(data = ., family = bernoulli,
              suspected_vap_plus_abx_and_pathogen_and_signs ~ (1 | subject_id) + max_read_prop_scaled)
  
  
  d_vap %>%
    brm(data = ., family = bernoulli,
        suspected_vap_plus_abx_and_pathogen_and_signs ~ (1 | subject_id) + max_read_prop_scaled,
        prior = c(prior(student_t(3,0,2.5), class = Intercept),
                  prior(student_t(3,0,2.5), class = b)
        ),
        sample_prior = "only", # prior predictive
        iter = 2000,
        warmup = 1000,
        chains = 4,
        cores = 4,
        control = list("adapt_delta" = 0.99, max_treedepth = 16),
        backend = "cmdstanr",
        seed = 16) -> d_prior
  
  
  d_prior %>%
    posterior_epred() %>%
    as_tibble() %>%
    gather(key = "observations", value = "pp_samples") %>%
    qplot(data = ., x = pp_samples, geom = "density")
  
  
  
  #' brms
  d_vap %>%
    brm(data = ., family = bernoulli,
        suspected_vap_plus_abx_and_pathogen_and_signs ~ (1 | subject_id) + max_read_prop_scaled,
        prior = c(prior(student_t(3,0,2.5), class = Intercept),
                  prior(student_t(3,0,2.5), class = b)
        ),
        #sample_prior = "only", # prior predictive
        iter = 2000,
        warmup = 1000,
        chains = 4,
        cores = 4,
        control = list("adapt_delta" = 0.999, max_treedepth = 16),
        backend = "cmdstanr",
        seed = 16) -> m_suspected_vap_plus_abx_and_pathogen_and_signs_max_read_prop_brms
}

if(save_models) {
  m_suspected_vap_plus_abx_and_pathogen_and_signs_max_read_prop_brms %>% write_rds(file = "./models/m_suspected_vap_plus_abx_and_pathogen_and_signs_max_read_prop_brms_mixed.rds.gz", compress = "gz")
  m_suspected_vap_plus_abx_and_pathogen_and_signs_max_read_prop_brms$fit %>% write_rds(file = "./models/m_suspected_vap_plus_abx_and_pathogen_and_signs_max_read_prop_brms_stanfit_mixed.rds.gz", compress = "gz")
}

m_suspected_vap_plus_abx_and_pathogen_and_signs_max_read_prop_brms <- read_rds(file = "./models/m_suspected_vap_plus_abx_and_pathogen_and_signs_max_read_prop_brms_mixed.rds.gz")

m_suspected_vap_plus_abx_and_pathogen_and_signs_max_read_prop_brms$formula

pp_check(m_suspected_vap_plus_abx_and_pathogen_and_signs_max_read_prop_brms)
rstan::check_hmc_diagnostics(m_suspected_vap_plus_abx_and_pathogen_and_signs_max_read_prop_brms$fit)


#' model results
m_suspected_vap_plus_abx_and_pathogen_and_signs_max_read_prop_brms %>%
  brms::posterior_summary() %>%
  as_tibble(rownames = "param") %>%
  mutate_if(is.numeric, list("OR" = ~ exp(.x))) %>%
  identity() -> m_suspected_vap_plus_abx_and_pathogen_and_signs_max_read_prop_brms_param_est
m_suspected_vap_plus_abx_and_pathogen_and_signs_max_read_prop_brms_param_est

if(run_plots) {
  m_suspected_vap_plus_abx_and_pathogen_and_signs_max_read_prop_brms_param_est %>%
    filter(grepl("b_",param)) %>%
    ggplot(data = .) + 
    geom_segment(aes(x = Q2.5_OR, xend = Q97.5_OR, y = param, yend = param)) +
    geom_point(aes(x = Estimate_OR, y = param)) +
    geom_vline(xintercept = 1, linetype = 2)
}


#' mixed model plot: by subject
if(run_plots) {
  m_suspected_vap_plus_abx_and_pathogen_and_signs_max_read_prop_brms$data %>%
    as_tibble() %>%
    expand(max_read_prop_scaled = modelr::seq_range(max_read_prop_scaled, n = 100), subject_id = unique(subject_id)) %>%
    mutate(max_read_prop = max_read_prop_scaled * sd(d_vap$max_read_prop, na.rm = TRUE) + mean(d_vap$max_read_prop, na.rm = TRUE)) %>%
    tidybayes::add_fitted_draws(newdata = ., model = m_suspected_vap_plus_abx_and_pathogen_and_signs_max_read_prop_brms) %>%
    ungroup() %>%
    ggplot(data = ., aes(x = max_read_prop_scaled, y = .value)) +
    tidybayes::stat_lineribbon(.width = c(0.5, 0.8, 0.95), alpha = 0.7, color = colorspace::sequential_hcl(5, palette = "Blues 3")[1]) +
    colorspace::scale_fill_discrete_sequential(palette = "Blues 3", nmax = 5, order = 2:4) +
    facet_wrap(~ subject_id)
}


#' mixed model plot: marginal
if(run_plots) {
  m_suspected_vap_plus_abx_and_pathogen_and_signs_max_read_prop_brms$data %>%
    as_tibble() %>%
    expand(max_read_prop_scaled = modelr::seq_range(max_read_prop_scaled, n = 100), subject_id = unique(subject_id)) %>%
    mutate(max_read_prop = max_read_prop_scaled * sd(d_vap$max_read_prop, na.rm = TRUE) + mean(d_vap$max_read_prop, na.rm = TRUE)) %>%
    tidybayes::add_fitted_draws(newdata = ., model = m_suspected_vap_plus_abx_and_pathogen_and_signs_max_read_prop_brms, re_formula = NULL, seed = 16) %>%
    ungroup() %>%
    ggplot(data = ., aes(x = max_read_prop, y = .value)) +
    tidybayes::stat_lineribbon(.width = c(0.5, 0.8, 0.95), alpha = 0.7, color = colorspace::sequential_hcl(5, palette = "Blues 3")[1]) +
    colorspace::scale_fill_discrete_sequential(palette = "Blues 3", nmax = 5, order = 2:4)
}










#' #######################################################################
#' #######################################################################
#' VAP as a function of lower respiratory tract total bacterial abundance: log_copy_16S_per_ml_sputum
#' #######################################################################
#' #######################################################################


#' #######################################################################
#' SUSPECTED VAP
#' #######################################################################

dat %>%
  # filter to include dates only up to VAP (or last date, if no VAP)
  group_by(subject_id) %>%
  mutate(first_suspected_vap = ifelse(is.na(first_suspected_vap), max(date), first_suspected_vap)) %>%
  filter(date <= first_suspected_vap) %>%
  ungroup() %>%
  select(subject_id, subject_day, log_copy_16S_per_ml_sputum, max_read_prop, log_copy_16S_per_ml_sputum, suspected_vap) %>%
  mutate_at(.vars = vars(-subject_id, -subject_day, -suspected_vap), .funs = list("scaled" = ~ as.vector(scale(.x)))) %>%
  distinct() %>%
  filter(complete.cases(.)) %>%
  identity() -> d_vap

d_vap


if(run_models) {
  d_vap %>%
    get_prior(data = ., family = bernoulli,
              suspected_vap ~ (1 | subject_id) + log_copy_16S_per_ml_sputum_scaled)
  
  
  d_vap %>%
    brm(data = ., family = bernoulli,
        suspected_vap ~ (1 | subject_id) + log_copy_16S_per_ml_sputum_scaled,
        prior = c(prior(student_t(3,0,2.5), class = Intercept),
                  prior(student_t(3,0,2.5), class = b)
        ),
        sample_prior = "only", # prior predictive
        iter = 2000,
        warmup = 1000,
        chains = 4,
        cores = 4,
        control = list("adapt_delta" = 0.99, max_treedepth = 16),
        backend = "cmdstanr",
        seed = 16) -> d_prior
  
  
  d_prior %>%
    posterior_epred() %>%
    as_tibble() %>%
    gather(key = "observations", value = "pp_samples") %>%
    qplot(data = ., x = pp_samples, geom = "density")
  
  
  
  #' brms
  d_vap %>%
    brm(data = ., family = bernoulli,
        suspected_vap ~ (1 | subject_id) + log_copy_16S_per_ml_sputum_scaled,
        prior = c(prior(student_t(3,0,2.5), class = Intercept),
                  prior(student_t(3,0,2.5), class = b)
        ),
        #sample_prior = "only", # prior predictive
        iter = 2000,
        warmup = 1000,
        chains = 4,
        cores = 4,
        control = list("adapt_delta" = 0.99, max_treedepth = 16),
        backend = "cmdstanr",
        seed = 16) -> m_suspected_vap_log_copy_16S_per_ml_sputum_brms
}

if(save_models) {
  m_suspected_vap_log_copy_16S_per_ml_sputum_brms %>% write_rds(file = "./models/m_suspected_vap_log_copy_16S_per_ml_sputum_brms_mixed.rds.gz", compress = "gz")
  m_suspected_vap_log_copy_16S_per_ml_sputum_brms$fit %>% write_rds(file = "./models/m_suspected_vap_log_copy_16S_per_ml_sputum_brms_stanfit_mixed.rds.gz", compress = "gz")
}

m_suspected_vap_log_copy_16S_per_ml_sputum_brms <- read_rds(file = "./models/m_suspected_vap_log_copy_16S_per_ml_sputum_brms_mixed.rds.gz")

m_suspected_vap_log_copy_16S_per_ml_sputum_brms$formula

pp_check(m_suspected_vap_log_copy_16S_per_ml_sputum_brms)
rstan::check_hmc_diagnostics(m_suspected_vap_log_copy_16S_per_ml_sputum_brms$fit)


#' model results
m_suspected_vap_log_copy_16S_per_ml_sputum_brms %>%
  brms::posterior_summary() %>%
  as_tibble(rownames = "param") %>%
  mutate_if(is.numeric, list("OR" = ~ exp(.x))) %>%
  identity() -> m_suspected_vap_log_copy_16S_per_ml_sputum_brms_param_est
m_suspected_vap_log_copy_16S_per_ml_sputum_brms_param_est

if(run_plots) {
  m_suspected_vap_log_copy_16S_per_ml_sputum_brms_param_est %>%
    filter(grepl("b_",param)) %>%
    ggplot(data = .) + 
    geom_segment(aes(x = Q2.5_OR, xend = Q97.5_OR, y = param, yend = param)) +
    geom_point(aes(x = Estimate_OR, y = param)) +
    geom_vline(xintercept = 1, linetype = 2)
}


#' mixed model plot: by subject
if(run_plots) {
  m_suspected_vap_log_copy_16S_per_ml_sputum_brms$data %>%
    as_tibble() %>%
    expand(log_copy_16S_per_ml_sputum_scaled = modelr::seq_range(log_copy_16S_per_ml_sputum_scaled, n = 100), subject_id = unique(subject_id)) %>%
    mutate(log_copy_16S_per_ml_sputum = log_copy_16S_per_ml_sputum_scaled * sd(d_vap$log_copy_16S_per_ml_sputum, na.rm = TRUE) + mean(d_vap$log_copy_16S_per_ml_sputum, na.rm = TRUE)) %>%
    tidybayes::add_fitted_draws(newdata = ., model = m_suspected_vap_log_copy_16S_per_ml_sputum_brms) %>%
    ungroup() %>%
    ggplot(data = ., aes(x = log_copy_16S_per_ml_sputum_scaled, y = .value)) +
    tidybayes::stat_lineribbon(.width = c(0.5, 0.8, 0.95), alpha = 0.7, color = colorspace::sequential_hcl(5, palette = "Blues 3")[1]) +
    colorspace::scale_fill_discrete_sequential(palette = "Blues 3", nmax = 5, order = 2:4) +
    facet_wrap(~ subject_id)
}


#' mixed model plot: marginal
if(run_plots) {
  m_suspected_vap_log_copy_16S_per_ml_sputum_brms$data %>%
    as_tibble() %>%
    expand(log_copy_16S_per_ml_sputum_scaled = modelr::seq_range(log_copy_16S_per_ml_sputum_scaled, n = 100), subject_id = unique(subject_id)) %>%
    mutate(log_copy_16S_per_ml_sputum = log_copy_16S_per_ml_sputum_scaled * sd(d_vap$log_copy_16S_per_ml_sputum, na.rm = TRUE) + mean(d_vap$log_copy_16S_per_ml_sputum, na.rm = TRUE)) %>%
    tidybayes::add_fitted_draws(newdata = ., model = m_suspected_vap_log_copy_16S_per_ml_sputum_brms, re_formula = NULL, seed = 16) %>%
    ungroup() %>%
    ggplot(data = ., aes(x = log_copy_16S_per_ml_sputum, y = .value)) +
    tidybayes::stat_lineribbon(.width = c(0.5, 0.8, 0.95), alpha = 0.7, color = colorspace::sequential_hcl(5, palette = "Blues 3")[1]) +
    colorspace::scale_fill_discrete_sequential(palette = "Blues 3", nmax = 5, order = 2:4)
}







#' #######################################################################
#' SUSPECTED VAP + CX + SIGNS (same as SUSPECTED VAP + CX)
#' #######################################################################

dat %>%
  # filter to include dates only up to VAP (or last date, if no VAP)
  group_by(subject_id) %>%
  mutate(first_suspected_vap_plus_pathogen = ifelse(is.na(first_suspected_vap_plus_pathogen), max(date), first_suspected_vap_plus_pathogen)) %>%
  filter(date <= first_suspected_vap_plus_pathogen) %>%
  ungroup() %>%
  select(subject_id, subject_day, log_copy_16S_per_ml_sputum, max_read_prop, log_copy_16S_per_ml_sputum, suspected_vap_plus_signs_and_pathogen) %>%
  mutate_at(.vars = vars(-subject_id, -subject_day, -suspected_vap_plus_signs_and_pathogen), .funs = list("scaled" = ~ as.vector(scale(.x)))) %>%
  distinct() %>%
  filter(complete.cases(.)) %>%
  identity() -> d_vap

d_vap


if (run_models) {
  d_vap %>%
    get_prior(data = ., family = bernoulli,
              suspected_vap_plus_signs_and_pathogen ~ (1 | subject_id) + log_copy_16S_per_ml_sputum_scaled)
  
  
  d_vap %>%
    brm(data = ., family = bernoulli,
        suspected_vap_plus_signs_and_pathogen ~ (1 | subject_id) + log_copy_16S_per_ml_sputum_scaled,
        prior = c(prior(student_t(3,0,2.5), class = Intercept),
                  prior(student_t(3,0,2.5), class = b)
        ),
        sample_prior = "only", # prior predictive
        iter = 2000,
        warmup = 1000,
        chains = 4,
        cores = 4,
        control = list("adapt_delta" = 0.99, max_treedepth = 16),
        backend = "cmdstanr",
        seed = 16) -> d_prior
  
  
  d_prior %>%
    posterior_epred() %>%
    as_tibble() %>%
    gather(key = "observations", value = "pp_samples") %>%
    qplot(data = ., x = pp_samples, geom = "density")
  
  
  
  #' brms
  d_vap %>%
    brm(data = ., family = bernoulli,
        suspected_vap_plus_signs_and_pathogen ~ (1 | subject_id) + log_copy_16S_per_ml_sputum_scaled,
        prior = c(prior(student_t(3,0,2.5), class = Intercept),
                  prior(student_t(3,0,2.5), class = b)
        ),
        #sample_prior = "only", # prior predictive
        iter = 2000,
        warmup = 1000,
        chains = 4,
        cores = 4,
        control = list("adapt_delta" = 0.99, max_treedepth = 16),
        backend = "cmdstanr",
        seed = 16) -> m_suspected_vap_plus_signs_and_pathogen_log_copy_16S_per_ml_sputum_brms
}

if(save_models) {
  m_suspected_vap_plus_signs_and_pathogen_log_copy_16S_per_ml_sputum_brms %>% write_rds(file = "./models/m_suspected_vap_plus_signs_and_pathogen_log_copy_16S_per_ml_sputum_brms_mixed.rds.gz", compress = "gz")
  m_suspected_vap_plus_signs_and_pathogen_log_copy_16S_per_ml_sputum_brms$fit %>% write_rds(file = "./models/m_suspected_vap_plus_signs_and_pathogen_log_copy_16S_per_ml_sputum_brms_stanfit_mixed.rds.gz", compress = "gz")
}

m_suspected_vap_plus_signs_and_pathogen_log_copy_16S_per_ml_sputum_brms <- read_rds(file = "./models/m_suspected_vap_plus_signs_and_pathogen_log_copy_16S_per_ml_sputum_brms_mixed.rds.gz")

m_suspected_vap_plus_signs_and_pathogen_log_copy_16S_per_ml_sputum_brms$formula

pp_check(m_suspected_vap_plus_signs_and_pathogen_log_copy_16S_per_ml_sputum_brms)
rstan::check_hmc_diagnostics(m_suspected_vap_plus_signs_and_pathogen_log_copy_16S_per_ml_sputum_brms$fit)


#' model results
m_suspected_vap_plus_signs_and_pathogen_log_copy_16S_per_ml_sputum_brms %>%
  brms::posterior_summary() %>%
  as_tibble(rownames = "param") %>%
  mutate_if(is.numeric, list("OR" = ~ exp(.x))) %>%
  identity() -> m_suspected_vap_plus_signs_and_pathogen_log_copy_16S_per_ml_sputum_brms_param_est
m_suspected_vap_plus_signs_and_pathogen_log_copy_16S_per_ml_sputum_brms_param_est

if(run_plots) {
  m_suspected_vap_plus_signs_and_pathogen_log_copy_16S_per_ml_sputum_brms_param_est %>%
    filter(grepl("b_",param)) %>%
    ggplot(data = .) + 
    geom_segment(aes(x = Q2.5_OR, xend = Q97.5_OR, y = param, yend = param)) +
    geom_point(aes(x = Estimate_OR, y = param)) +
    geom_vline(xintercept = 1, linetype = 2)
}


#' mixed model plot: by subject
if(run_plots) {
  m_suspected_vap_plus_signs_and_pathogen_log_copy_16S_per_ml_sputum_brms$data %>%
    as_tibble() %>%
    expand(log_copy_16S_per_ml_sputum_scaled = modelr::seq_range(log_copy_16S_per_ml_sputum_scaled, n = 100), subject_id = unique(subject_id)) %>%
    mutate(log_copy_16S_per_ml_sputum = log_copy_16S_per_ml_sputum_scaled * sd(d_vap$log_copy_16S_per_ml_sputum, na.rm = TRUE) + mean(d_vap$log_copy_16S_per_ml_sputum, na.rm = TRUE)) %>%
    tidybayes::add_fitted_draws(newdata = ., model = m_suspected_vap_plus_signs_and_pathogen_log_copy_16S_per_ml_sputum_brms) %>%
    ungroup() %>%
    ggplot(data = ., aes(x = log_copy_16S_per_ml_sputum_scaled, y = .value)) +
    tidybayes::stat_lineribbon(.width = c(0.5, 0.8, 0.95), alpha = 0.7, color = colorspace::sequential_hcl(5, palette = "Blues 3")[1]) +
    colorspace::scale_fill_discrete_sequential(palette = "Blues 3", nmax = 5, order = 2:4) +
    facet_wrap(~ subject_id)
}


#' mixed model plot: marginal
if(run_plots) {
  m_suspected_vap_plus_signs_and_pathogen_log_copy_16S_per_ml_sputum_brms$data %>%
    as_tibble() %>%
    expand(log_copy_16S_per_ml_sputum_scaled = modelr::seq_range(log_copy_16S_per_ml_sputum_scaled, n = 100), subject_id = unique(subject_id)) %>%
    mutate(log_copy_16S_per_ml_sputum = log_copy_16S_per_ml_sputum_scaled * sd(d_vap$log_copy_16S_per_ml_sputum, na.rm = TRUE) + mean(d_vap$log_copy_16S_per_ml_sputum, na.rm = TRUE)) %>%
    tidybayes::add_fitted_draws(newdata = ., model = m_suspected_vap_plus_signs_and_pathogen_log_copy_16S_per_ml_sputum_brms, re_formula = NULL, seed = 16) %>%
    ungroup() %>%
    ggplot(data = ., aes(x = log_copy_16S_per_ml_sputum, y = .value)) +
    tidybayes::stat_lineribbon(.width = c(0.5, 0.8, 0.95), alpha = 0.7, color = colorspace::sequential_hcl(5, palette = "Blues 3")[1]) +
    colorspace::scale_fill_discrete_sequential(palette = "Blues 3", nmax = 5, order = 2:4)
}




#' #######################################################################
#' SUSPECTED VAP + CX + SIGNS + ABX
#' #######################################################################

dat %>%
  # filter to include dates only up to VAP (or last date, if no VAP)
  group_by(subject_id) %>%
  mutate(first_suspected_vap_plus_abx_and_pathogen_and_signs = ifelse(is.na(first_suspected_vap_plus_abx_and_pathogen_and_signs), max(date), first_suspected_vap_plus_abx_and_pathogen_and_signs)) %>%
  filter(date <= first_suspected_vap_plus_abx_and_pathogen_and_signs) %>%
  ungroup() %>%
  select(subject_id, subject_day, log_copy_16S_per_ml_sputum, max_read_prop, log_copy_16S_per_ml_sputum, suspected_vap_plus_abx_and_pathogen_and_signs) %>%
  mutate_at(.vars = vars(-subject_id, -subject_day, -suspected_vap_plus_abx_and_pathogen_and_signs), .funs = list("scaled" = ~ as.vector(scale(.x)))) %>%
  distinct() %>%
  filter(complete.cases(.)) %>%
  identity() -> d_vap

d_vap


if (run_models) {
  d_vap %>%
    get_prior(data = ., family = bernoulli,
              suspected_vap_plus_abx_and_pathogen_and_signs ~ (1 | subject_id) + log_copy_16S_per_ml_sputum_scaled)
  
  
  d_vap %>%
    brm(data = ., family = bernoulli,
        suspected_vap_plus_abx_and_pathogen_and_signs ~ (1 | subject_id) + log_copy_16S_per_ml_sputum_scaled,
        prior = c(prior(student_t(3,0,2.5), class = Intercept),
                  prior(student_t(3,0,2.5), class = b)
        ),
        sample_prior = "only", # prior predictive
        iter = 2000,
        warmup = 1000,
        chains = 4,
        cores = 4,
        control = list("adapt_delta" = 0.99, max_treedepth = 16),
        backend = "cmdstanr",
        seed = 16) -> d_prior
  
  
  d_prior %>%
    posterior_epred() %>%
    as_tibble() %>%
    gather(key = "observations", value = "pp_samples") %>%
    qplot(data = ., x = pp_samples, geom = "density")
  
  
  
  #' brms
  d_vap %>%
    brm(data = ., family = bernoulli,
        suspected_vap_plus_abx_and_pathogen_and_signs ~ (1 | subject_id) + log_copy_16S_per_ml_sputum_scaled,
        prior = c(prior(student_t(3,0,2.5), class = Intercept),
                  prior(student_t(3,0,2.5), class = b)
        ),
        #sample_prior = "only", # prior predictive
        iter = 2000,
        warmup = 1000,
        chains = 4,
        cores = 4,
        control = list("adapt_delta" = 0.999, max_treedepth = 16),
        backend = "cmdstanr",
        seed = 16) -> m_suspected_vap_plus_abx_and_pathogen_and_signs_log_copy_16S_per_ml_sputum_brms
}

if(save_models) {
  m_suspected_vap_plus_abx_and_pathogen_and_signs_log_copy_16S_per_ml_sputum_brms %>% write_rds(file = "./models/m_suspected_vap_plus_abx_and_pathogen_and_signs_log_copy_16S_per_ml_sputum_brms_mixed.rds.gz", compress = "gz")
  m_suspected_vap_plus_abx_and_pathogen_and_signs_log_copy_16S_per_ml_sputum_brms$fit %>% write_rds(file = "./models/m_suspected_vap_plus_abx_and_pathogen_and_signs_log_copy_16S_per_ml_sputum_brms_stanfit_mixed.rds.gz", compress = "gz")
}

m_suspected_vap_plus_abx_and_pathogen_and_signs_log_copy_16S_per_ml_sputum_brms <- read_rds(file = "./models/m_suspected_vap_plus_abx_and_pathogen_and_signs_log_copy_16S_per_ml_sputum_brms_mixed.rds.gz")

m_suspected_vap_plus_abx_and_pathogen_and_signs_log_copy_16S_per_ml_sputum_brms$formula

pp_check(m_suspected_vap_plus_abx_and_pathogen_and_signs_log_copy_16S_per_ml_sputum_brms)
rstan::check_hmc_diagnostics(m_suspected_vap_plus_abx_and_pathogen_and_signs_log_copy_16S_per_ml_sputum_brms$fit)


#' model results
m_suspected_vap_plus_abx_and_pathogen_and_signs_log_copy_16S_per_ml_sputum_brms %>%
  brms::posterior_summary() %>%
  as_tibble(rownames = "param") %>%
  mutate_if(is.numeric, list("OR" = ~ exp(.x))) %>%
  identity() -> m_suspected_vap_plus_abx_and_pathogen_and_signs_log_copy_16S_per_ml_sputum_brms_param_est
m_suspected_vap_plus_abx_and_pathogen_and_signs_log_copy_16S_per_ml_sputum_brms_param_est

if(run_plots) {
  m_suspected_vap_plus_abx_and_pathogen_and_signs_log_copy_16S_per_ml_sputum_brms_param_est %>%
    filter(grepl("b_",param)) %>%
    ggplot(data = .) + 
    geom_segment(aes(x = Q2.5_OR, xend = Q97.5_OR, y = param, yend = param)) +
    geom_point(aes(x = Estimate_OR, y = param)) +
    geom_vline(xintercept = 1, linetype = 2)
}


#' mixed model plot: by subject
if(run_plots) {
  m_suspected_vap_plus_abx_and_pathogen_and_signs_log_copy_16S_per_ml_sputum_brms$data %>%
    as_tibble() %>%
    expand(log_copy_16S_per_ml_sputum_scaled = modelr::seq_range(log_copy_16S_per_ml_sputum_scaled, n = 100), subject_id = unique(subject_id)) %>%
    mutate(log_copy_16S_per_ml_sputum = log_copy_16S_per_ml_sputum_scaled * sd(d_vap$log_copy_16S_per_ml_sputum, na.rm = TRUE) + mean(d_vap$log_copy_16S_per_ml_sputum, na.rm = TRUE)) %>%
    tidybayes::add_fitted_draws(newdata = ., model = m_suspected_vap_plus_abx_and_pathogen_and_signs_log_copy_16S_per_ml_sputum_brms) %>%
    ungroup() %>%
    ggplot(data = ., aes(x = log_copy_16S_per_ml_sputum_scaled, y = .value)) +
    tidybayes::stat_lineribbon(.width = c(0.5, 0.8, 0.95), alpha = 0.7, color = colorspace::sequential_hcl(5, palette = "Blues 3")[1]) +
    colorspace::scale_fill_discrete_sequential(palette = "Blues 3", nmax = 5, order = 2:4) +
    facet_wrap(~ subject_id)
}


#' mixed model plot: marginal
if(run_plots) {
  m_suspected_vap_plus_abx_and_pathogen_and_signs_log_copy_16S_per_ml_sputum_brms$data %>%
    as_tibble() %>%
    expand(log_copy_16S_per_ml_sputum_scaled = modelr::seq_range(log_copy_16S_per_ml_sputum_scaled, n = 100), subject_id = unique(subject_id)) %>%
    mutate(log_copy_16S_per_ml_sputum = log_copy_16S_per_ml_sputum_scaled * sd(d_vap$log_copy_16S_per_ml_sputum, na.rm = TRUE) + mean(d_vap$log_copy_16S_per_ml_sputum, na.rm = TRUE)) %>%
    tidybayes::add_fitted_draws(newdata = ., model = m_suspected_vap_plus_abx_and_pathogen_and_signs_log_copy_16S_per_ml_sputum_brms, re_formula = NULL, seed = 16) %>%
    ungroup() %>%
    ggplot(data = ., aes(x = log_copy_16S_per_ml_sputum, y = .value)) +
    tidybayes::stat_lineribbon(.width = c(0.5, 0.8, 0.95), alpha = 0.7, color = colorspace::sequential_hcl(5, palette = "Blues 3")[1]) +
    colorspace::scale_fill_discrete_sequential(palette = "Blues 3", nmax = 5, order = 2:4)
}









#' #######################################################################
#' COMPARE MODELS ACROSS MDIs & OUTCOME DEFINITIONS
#' #######################################################################

outcome_names <- c("Suspected VA-LRTI",
                   "VA-LRTI with Clinical Signs<br>& Positive Culture",
                   "VA-LRTI with Clinical Signs<br>& Positive Culture<br>& Antibiotic Change")

exposure_names <- c("Shannon Diversity<br>(base e)",
                    "Maximum ASV<br>Proportional Abundance",
                    "Total Bacterial Abundance<br>by 16S rRNA Gene qPCR")


tibble(model = c("m_suspected_vap_shannon_brms",
                 "m_suspected_vap_plus_signs_and_pathogen_shannon_brms",
                 "m_suspected_vap_plus_abx_and_pathogen_and_signs_shannon_brms",
                 "m_suspected_vap_max_read_prop_brms",
                 "m_suspected_vap_plus_signs_and_pathogen_max_read_prop_brms",
                 "m_suspected_vap_plus_abx_and_pathogen_and_signs_max_read_prop_brms",
                 "m_suspected_vap_log_copy_16S_per_ml_sputum_brms",
                 "m_suspected_vap_plus_signs_and_pathogen_log_copy_16S_per_ml_sputum_brms",
                 "m_suspected_vap_plus_abx_and_pathogen_and_signs_log_copy_16S_per_ml_sputum_brms"
                 )
       ) %>%
  mutate(outcome_definition = gsub("m_|_brms","", model)) %>%
  mutate(outcome_definition = case_when(grepl("suspected_vap_plus_signs_and_pathogen",model) ~ outcome_names[2],
                                        grepl("suspected_vap_plus_abx_and_pathogen_and_signs",model) ~ outcome_names[3],
                                        grepl("suspected_vap_plus_signs_and_pathogen|suspected_vap_plus_abx_and_pathogen_and_signs",model) == FALSE ~ outcome_names[1])) %>%
  mutate(outcome_definition = factor(outcome_definition, levels = outcome_names)) %>%
  mutate(parameter = stringr::str_extract(string = model, pattern = "log_copy_16S_per_ml_sputum|shannon|max_read_prop")) %>%
  mutate(parameter = case_when(parameter == "shannon" ~ exposure_names[1],
                               parameter == "max_read_prop" ~ exposure_names[2],
                               parameter == "log_copy_16S_per_ml_sputum" ~ exposure_names[3])) %>%
  mutate(parameter = factor(parameter, levels = exposure_names)) %>%
  mutate(brms_model = map(.x = model, .f = ~ get(.x))) %>%
  mutate(fit = map(.x = brms_model, .f = ~ .x$fit)) %>%
  mutate(posterior = map(.x = fit, .f = ~ posterior::as_draws_df(.x) %>% as_tibble() %>% select(matches("b_shannon_scaled|b_max_read_prop_scaled|b_log_copy_16S_per_ml_sputum_scaled")) %>% pull())) %>%
  identity() -> crossectional_models
crossectional_models



crossectional_models %>%
  select(outcome_definition, parameter, posterior) %>%
  unnest(cols = c(posterior)) %>%
  ungroup() %>%
  mutate(posterior_OR = exp(posterior)) %>% 
  ggplot(data = ., aes(y = outcome_definition, x = posterior_OR, fill = stat(x > 1))) +
  tidybayes::stat_halfeye() +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_fill_manual(values = c("gray80","skyblue")) +
  facet_wrap(facets = ~ parameter, scales = "free_x") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = ggtext::element_markdown(color = "black"),
        axis.text.y = ggtext::element_markdown(color = "black"),
        axis.title.x = ggtext::element_markdown(color = "black"),
        axis.title.y = ggtext::element_markdown(color = "black"),
        strip.text.x = ggtext::element_markdown(color = "black"),
        strip.text.y = ggtext::element_markdown(color = "black"),
        strip.background = element_blank(),
  ) +
  labs(x = "Odds Ratio VA-LRTI",
       y = "")





#' #######################################################################
#' COMPARE MDIs with PRIMARY OUTCOME DEFINITION
#' #######################################################################

if(run_plots) {
  
counterfact <- function(brms_mod, exp_var_name = "shannon_scaled", draw_number = NULL) {
  brms_mod[["data"]] %>%
    expand(exp_var = modelr::seq_range(brms_mod[["data"]][[exp_var_name]], n = 100),
           subject_id = unique(subject_id)) %>%
    rename_with(.cols = exp_var, .fn = ~ exp_var_name) %>%
    tidybayes::add_fitted_draws(model = brms_mod, re_formula = NULL, n = draw_number, seed = 16) %>% # n = NULL for all draws
    identity() -> post_fitted
  return(post_fitted)
}


outcome_names <- c("Suspected VA-LRTI",
                   "VA-LRTI with Clinical Signs<br>& Positive Culture",
                   "VA-LRTI with Clinical Signs<br>& Positive Culture<br>& Antibiotic Change")

exposure_names <- c("Shannon Diversity (base e)",
                    "Maximum ASV<br>Proportional Abundance",
                    "Total Bacterial Abundance<br>by 16S rRNA Gene qPCR<br>(log copies per mL sputum)")


tibble(model = c("m_suspected_vap_shannon_brms",
                 "m_suspected_vap_plus_signs_and_pathogen_shannon_brms",
                 "m_suspected_vap_plus_abx_and_pathogen_and_signs_shannon_brms",
                 "m_suspected_vap_max_read_prop_brms",
                 "m_suspected_vap_plus_signs_and_pathogen_max_read_prop_brms",
                 "m_suspected_vap_plus_abx_and_pathogen_and_signs_max_read_prop_brms",
                 "m_suspected_vap_log_copy_16S_per_ml_sputum_brms",
                 "m_suspected_vap_plus_signs_and_pathogen_log_copy_16S_per_ml_sputum_brms",
                 "m_suspected_vap_plus_abx_and_pathogen_and_signs_log_copy_16S_per_ml_sputum_brms"
)
) %>%
  mutate(outcome_definition = gsub("m_|_brms","", model)) %>%
  mutate(outcome_definition = case_when(grepl("suspected_vap_plus_signs_and_pathogen",model) ~ outcome_names[2],
                                        grepl("suspected_vap_plus_abx_and_pathogen_and_signs",model) ~ outcome_names[3],
                                        grepl("suspected_vap_plus_signs_and_pathogen|suspected_vap_plus_abx_and_pathogen_and_signs",model) == FALSE ~ outcome_names[1])) %>%
  mutate(outcome_definition = factor(outcome_definition, levels = outcome_names)) %>%
  mutate(parameter = stringr::str_extract(string = model, pattern = "log_copy_16S_per_ml_sputum|shannon|max_read_prop")) %>%
  mutate(param = paste0(parameter,"_scaled")) %>%
  mutate(parameter = case_when(parameter == "shannon" ~ exposure_names[1],
                               parameter == "max_read_prop" ~ exposure_names[2],
                               parameter == "log_copy_16S_per_ml_sputum" ~ exposure_names[3])) %>%
  mutate(parameter = factor(parameter, levels = exposure_names)) %>%
  mutate(brms_model = map(.x = model, .f = ~ get(.x))) %>%
  mutate(fit = map(.x = brms_model, .f = ~ .x$fit)) %>%
  mutate(dat = map(.x = brms_model, .f = ~ .$dat %>% as_tibble())) %>%
  mutate(post_fit = map2(.x = brms_model, .y = param, .f = ~ counterfact(brms_mod = .x, exp_var_name = .y, draw_number = 100))) %>%
  mutate(post_fit = map(.x = post_fit, .f = ~ rename_at(.tbl = .x, .vars = vars(matches("scaled")), .funs = ~ "exp_var"))) %>%
  identity() -> crossectional_post_fit
crossectional_post_fit


#' review data for scale correction
dat %>%
  # filter to include dates only up to VAP (or last date, if no VAP)
  group_by(subject_id) %>%
  mutate(first_suspected_vap = ifelse(is.na(first_suspected_vap), max(date), first_suspected_vap)) %>%
  filter(date <= first_suspected_vap) %>%
  ungroup() %>%
  distinct() %>%
  #filter(complete.cases(.)) %>%
  identity() -> d_vap1

dat %>%
  # filter to include dates only up to VAP (or last date, if no VAP)
  group_by(subject_id) %>%
  mutate(first_suspected_vap_plus_pathogen = ifelse(is.na(first_suspected_vap_plus_pathogen), max(date), first_suspected_vap_plus_pathogen)) %>%
  filter(date <= first_suspected_vap_plus_pathogen) %>%
  ungroup() %>%
  distinct() %>%
  #filter(complete.cases(.)) %>%
  identity() -> d_vap2

dat %>%
  # filter to include dates only up to VAP (or last date, if no VAP)
  group_by(subject_id) %>%
  mutate(first_suspected_vap_plus_abx_and_pathogen_and_signs = ifelse(is.na(first_suspected_vap_plus_abx_and_pathogen_and_signs), max(date), first_suspected_vap_plus_abx_and_pathogen_and_signs)) %>%
  filter(date <= first_suspected_vap_plus_abx_and_pathogen_and_signs) %>%
  ungroup() %>%
  distinct() %>%
  #filter(complete.cases(.)) %>%
  identity() -> d_vap3





#' integrate data to rescale exposure
crossectional_post_fit %>%
  mutate(dat = case_when(outcome_definition == outcome_names[1] ~ list(d_vap1),
                         outcome_definition == outcome_names[2] ~ list(d_vap2),
                         outcome_definition == outcome_names[3] ~ list(d_vap3))) %>%
  mutate(d_mean_shannon = map_dbl(.x = dat, .f = ~ summarise_at(.tbl = .x, .vars = vars(matches("shannon$")), .funs = function(d) mean(d, na.rm = TRUE)) %>% pull()),
         d_sd_shannon = map_dbl(.x = dat, .f = ~ summarise_at(.tbl = .x, .vars = vars(matches("shannon$")), .funs = function(d) sd(d, na.rm = TRUE)) %>% pull()),
         d_mean_asv = map_dbl(.x = dat, .f = ~ summarise_at(.tbl = .x, .vars = vars(matches("max_read_prop$")), .funs = function(d) mean(d, na.rm = TRUE)) %>% pull()),
         d_sd_asv = map_dbl(.x = dat, .f = ~ summarise_at(.tbl = .x, .vars = vars(matches("max_read_prop$")), .funs = function(d) sd(d, na.rm = TRUE)) %>% pull()),
         d_mean_qpcr = map_dbl(.x = dat, .f = ~ summarise_at(.tbl = .x, .vars = vars(matches("log_copy_16S_per_ml_sputum$")), .funs = function(d) mean(d, na.rm = TRUE)) %>% pull()),
         d_sd_qpcr = map_dbl(.x = dat, .f = ~ summarise_at(.tbl = .x, .vars = vars(matches("log_copy_16S_per_ml_sputum$")), .funs = function(d) sd(d, na.rm = TRUE)) %>% pull()),
         ) %>%
  select(outcome_definition, parameter, post_fit, contains("d_mean"), contains("d_sd")) %>%
  unnest(cols = c(post_fit)) %>%
  ungroup() %>%
  mutate(exp_var_raw = case_when(parameter == exposure_names[1] ~ (exp_var * d_sd_shannon + d_mean_shannon),
                                 parameter == exposure_names[2] ~ (exp_var * d_sd_asv + d_mean_asv),
                                 parameter == exposure_names[3] ~ (exp_var * d_sd_qpcr + d_mean_qpcr)
         )) %>%
  identity() %>% {. -> crossectional_post_fit_rescale} %>% 
  # select primary outcome definition
  filter(outcome_definition == outcome_names[2]) %>%
  #
  ggplot(data = ., aes(x = exp_var_raw, y = .value)) +
  tidybayes::stat_lineribbon(.width = c(0.5,0.8,0.95),
                             alpha = 0.7,
                             color = colorspace::sequential_hcl(5, palette = "Blues 3")[1]) +
  colorspace::scale_fill_discrete_sequential(palette = "Blues 3", nmax = 5, order = 2:4) +
  facet_wrap(facets = ~ parameter, scales = "free_x") +
  theme_bw() +
  theme(legend.position = c(0.08,0.82),
        legend.direction = "vertical",
        axis.text.x = ggtext::element_markdown(color = "black"),
        axis.text.y = ggtext::element_markdown(color = "black"),
        axis.title.y = ggtext::element_markdown(),
        strip.text.y = ggtext::element_markdown(),
        strip.text.x = ggtext::element_markdown(),
        legend.title = ggtext::element_markdown(size = 10),
        legend.text = ggtext::element_markdown(size = 10),
        #legend.background = element_blank(),
        legend.background = element_rect(color = "black", fill = "white", size = 0.25),
        strip.background = element_blank()) +
  labs(x = "Microbiome Disruption Index", y = "Probability of VA-LRTI", fill = "Posterior<br>Credible<br>Interval") -> p_combined_mdi_crosssection_counterfact
p_combined_mdi_crosssection_counterfact

}


# p_combined_mdi_crosssection_counterfact %>%
#   ggsave(plot = ., filename = "./figs/p_combined_mdi_crosssection_counterfact.pdf", height = 5, width = 6, units = "in")
# p_combined_mdi_crosssection_counterfact %>%
#   ggsave(plot = ., filename = "./figs/p_combined_mdi_crosssection_counterfact.svg", height = 5, width = 6, units = "in")
# p_combined_mdi_crosssection_counterfact %>%
#   ggsave(plot = ., filename = "./figs/p_combined_mdi_crosssection_counterfact.png", height = 5, width = 6, units = "in", dpi = 600)

  
#' #######################################################################
#' POSTERIOR CONTRASTS & MODEL FIT BY LOO-CV
#' #######################################################################

outcome_names <- c("Suspected VA-LRTI",
                   "VA-LRTI with Clinical Signs<br>& Positive Culture",
                   "VA-LRTI with Clinical Signs<br>& Positive Culture<br>& Antibiotic Change")

exposure_names <- c("Shannon Diversity (base e)",
                    "Maximum ASV<br>Proportional Abundance",
                    "Total Bacterial Abundance<br>by 16S rRNA Gene qPCR<br>(log copies per mL sputum)")


#' summaries
tibble(model = c("m_suspected_vap_shannon_brms",
                 "m_suspected_vap_plus_signs_and_pathogen_shannon_brms",
                 "m_suspected_vap_plus_abx_and_pathogen_and_signs_shannon_brms",
                 "m_suspected_vap_max_read_prop_brms",
                 "m_suspected_vap_plus_signs_and_pathogen_max_read_prop_brms",
                 "m_suspected_vap_plus_abx_and_pathogen_and_signs_max_read_prop_brms",
                 "m_suspected_vap_log_copy_16S_per_ml_sputum_brms",
                 "m_suspected_vap_plus_signs_and_pathogen_log_copy_16S_per_ml_sputum_brms",
                 "m_suspected_vap_plus_abx_and_pathogen_and_signs_log_copy_16S_per_ml_sputum_brms"
)
) %>%
  mutate(outcome_definition = gsub("m_|_brms","", model)) %>%
  mutate(outcome_definition = case_when(grepl("suspected_vap_plus_signs_and_pathogen",model) ~ outcome_names[2],
                                        grepl("suspected_vap_plus_abx_and_pathogen_and_signs",model) ~ outcome_names[3],
                                        grepl("suspected_vap_plus_signs_and_pathogen|suspected_vap_plus_abx_and_pathogen_and_signs",model) == FALSE ~ outcome_names[1])) %>%
  mutate(outcome_definition = factor(outcome_definition, levels = outcome_names)) %>%
  mutate(parameter = stringr::str_extract(string = model, pattern = "log_copy_16S_per_ml_sputum|shannon|max_read_prop")) %>%
  mutate(param = paste0(parameter,"_scaled")) %>%
  mutate(parameter = case_when(parameter == "shannon" ~ exposure_names[1],
                               parameter == "max_read_prop" ~ exposure_names[2],
                               parameter == "log_copy_16S_per_ml_sputum" ~ exposure_names[3])) %>%
  mutate(parameter = factor(parameter, levels = exposure_names)) %>%
  mutate(brms_model = map(.x = model, .f = ~ get(.x))) %>%
  mutate(fit = map(.x = brms_model, .f = ~ .x$fit)) %>%
  mutate(dat = map(.x = brms_model, .f = ~ .$dat %>% as_tibble())) %>%
  mutate(post_sum = map(.x = brms_model, .f = ~ brms::posterior_summary(.x) %>% as_tibble(rownames = "post_sum_param"))) %>%
  select(model,outcome_definition, parameter, param, post_sum) %>%
  unnest(cols = c(post_sum)) %>%
  filter(grepl("^b_",post_sum_param)) %>%
  filter(grepl("Intercept",post_sum_param) == FALSE) %>%
  mutate_at(.vars = vars(Estimate, Est.Error, Q2.5, Q97.5), .funs = list("OR" = ~ exp(.x))) %>%
  identity() -> crossectional_post_summary
crossectional_post_summary


crossectional_post_summary %>%
  select(outcome_definition, parameter, contains("_OR")) %>%
  mutate_at(.vars = vars(outcome_definition, parameter), .funs = ~ gsub("<br>","\n",.x)) %>%
  gt::gt() %>%
  gt::fmt_number(columns = 3:6, n_sigfig = 3)


#' contrasts
#' Shannon contrast
m_suspected_vap_plus_signs_and_pathogen_shannon_brms$data %>%
  as_tibble() %>%
  mutate(shannon = shannon_scaled * sd(d_vap2$shannon, na.rm = TRUE) + mean(d_vap2$shannon, na.rm = TRUE)) %>%
  expand(shannon = c(mean(d_vap2$shannon, na.rm = TRUE) + sd(d_vap2$shannon, na.rm = TRUE),
                     mean(d_vap2$shannon, na.rm = TRUE),
                     mean(d_vap2$shannon, na.rm = TRUE) - sd(d_vap2$shannon, na.rm = TRUE)),
         subject_id = unique(subject_id)) %>%
  replicate(n = 10, ., simplify = FALSE) %>%
  bind_rows() %>%
  mutate(shannon_scaled = (shannon - mean(d_vap2$shannon, na.rm = TRUE)) / sd(d_vap2$shannon, na.rm = TRUE)) %>%
  tidybayes::add_fitted_draws(model = m_suspected_vap_plus_signs_and_pathogen_shannon_brms, re_formula = NULL, n = 100, seed = 16) %>% # n = NULL for all draws
  identity() -> m_svp_shannon_contrast
m_svp_shannon_contrast


m_svp_shannon_contrast %>%
  ungroup() %>%
  mutate(shannon = paste0("shannon_",round(shannon,2))) %>%
  select(shannon, .value) %>%
  pivot_wider(names_from = shannon, values_from = .value) %>%
  unnest(cols = contains("shannon")) %>%
  # contrast as you move from higher Shannon to lower Shannon
  mutate(contrast = shannon_3.53 - shannon_0.85) %>%
  tidybayes::median_qi(.width = .95) %>%
  pivot_longer(cols = matches("shannon|contrast"), names_to = "posterior", values_to = "estimate") %>%
  gt::gt() %>%
  gt::fmt_number(columns = 5, n_sigfig = 3) %>%
  gt::tab_header(title = "Shannon Diversity")
  


#' contrasts
#' max_prop contrast
m_suspected_vap_plus_signs_and_pathogen_max_read_prop_brms$data %>%
  as_tibble() %>%
  mutate(max_read_prop = max_read_prop_scaled * sd(d_vap2$max_read_prop, na.rm = TRUE) + mean(d_vap2$max_read_prop, na.rm = TRUE)) %>%
  expand(max_read_prop = c(mean(d_vap2$max_read_prop, na.rm = TRUE) + sd(d_vap2$max_read_prop, na.rm = TRUE),
                     mean(d_vap2$max_read_prop, na.rm = TRUE),
                     mean(d_vap2$max_read_prop, na.rm = TRUE) - sd(d_vap2$max_read_prop, na.rm = TRUE)),
         subject_id = unique(subject_id)) %>%
  replicate(n = 10, ., simplify = FALSE) %>%
  bind_rows() %>%
  mutate(max_read_prop_scaled = (max_read_prop - mean(d_vap2$max_read_prop, na.rm = TRUE)) / sd(d_vap2$max_read_prop, na.rm = TRUE)) %>%
  tidybayes::add_fitted_draws(model = m_suspected_vap_plus_signs_and_pathogen_max_read_prop_brms, re_formula = NULL, n = 100, seed = 16) %>% # n = NULL for all draws
  identity() -> m_svp_max_read_prop_contrast
m_svp_max_read_prop_contrast


m_svp_max_read_prop_contrast %>%
  ungroup() %>%
  mutate(max_read_prop = paste0("max_read_prop_",round(max_read_prop,2))) %>%
  select(max_read_prop, .value) %>%
  pivot_wider(names_from = max_read_prop, values_from = .value) %>%
  unnest(cols = contains("max_read_prop")) %>%
  # contrast as you move from higher max_read_prop to lower max_read_prop
  mutate(contrast = max_read_prop_0.81 - max_read_prop_0.29) %>%
  tidybayes::median_qi(.width = .95) %>%
  pivot_longer(cols = matches("max_read_prop|contrast"), names_to = "posterior", values_to = "estimate") %>%
  gt::gt() %>%
  gt::fmt_number(columns = 5, n_sigfig = 3) %>%
  gt::tab_header(title = "Maximum ASV Proportional Abundance")





#' contrasts
#' max_prop contrast
m_suspected_vap_plus_signs_and_pathogen_log_copy_16S_per_ml_sputum_brms$data %>%
  as_tibble() %>%
  mutate(log_copy_16S_per_ml_sputum = log_copy_16S_per_ml_sputum_scaled * sd(d_vap2$log_copy_16S_per_ml_sputum, na.rm = TRUE) + mean(d_vap2$log_copy_16S_per_ml_sputum, na.rm = TRUE)) %>%
  expand(log_copy_16S_per_ml_sputum = c(mean(d_vap2$log_copy_16S_per_ml_sputum, na.rm = TRUE) + sd(d_vap2$log_copy_16S_per_ml_sputum, na.rm = TRUE),
                           mean(d_vap2$log_copy_16S_per_ml_sputum, na.rm = TRUE),
                           mean(d_vap2$log_copy_16S_per_ml_sputum, na.rm = TRUE) - sd(d_vap2$log_copy_16S_per_ml_sputum, na.rm = TRUE)),
         subject_id = unique(subject_id)) %>%
  replicate(n = 10, ., simplify = FALSE) %>%
  bind_rows() %>%
  mutate(log_copy_16S_per_ml_sputum_scaled = (log_copy_16S_per_ml_sputum - mean(d_vap2$log_copy_16S_per_ml_sputum, na.rm = TRUE)) / sd(d_vap2$log_copy_16S_per_ml_sputum, na.rm = TRUE)) %>%
  tidybayes::add_fitted_draws(model = m_suspected_vap_plus_signs_and_pathogen_log_copy_16S_per_ml_sputum_brms, re_formula = NULL, n = 100, seed = 16) %>% # n = NULL for all draws
  identity() -> m_svp_log_copy_16S_per_ml_sputum_contrast
m_svp_log_copy_16S_per_ml_sputum_contrast


m_svp_log_copy_16S_per_ml_sputum_contrast %>%
  ungroup() %>%
  mutate(log_copy_16S_per_ml_sputum = paste0("log_copy_16S_per_ml_sputum_",round(log_copy_16S_per_ml_sputum,2))) %>%
  select(log_copy_16S_per_ml_sputum, .value) %>%
  pivot_wider(names_from = log_copy_16S_per_ml_sputum, values_from = .value) %>%
  unnest(cols = contains("log_copy_16S_per_ml_sputum")) %>%
  # contrast as you move from higher log_copy_16S_per_ml_sputum to lower log_copy_16S_per_ml_sputum
  mutate(contrast = log_copy_16S_per_ml_sputum_9.59 - log_copy_16S_per_ml_sputum_7.54) %>%
  tidybayes::median_qi(.width = .95) %>%
  pivot_longer(cols = matches("log_copy_16S_per_ml_sputum|contrast"), names_to = "posterior", values_to = "estimate") %>%
  gt::gt() %>%
  gt::fmt_number(columns = 5, n_sigfig = 3) %>%
  gt::tab_header(title = "Total Bacterial Abundance\n(log copy 16S rRNA per mL sputum)")

