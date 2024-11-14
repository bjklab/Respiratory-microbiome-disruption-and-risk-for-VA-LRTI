#' ###########################################
#' MIXED EFFECT VAP ~ CROSS-SECTIONAL MDI models
#' ###########################################



#' load libraries and set seed / control variables
library(tidyverse)
library(tidybayes)
library(brms)
library(gt)
set.seed(16)
run_models = TRUE
save_models = FALSE
run_plots = TRUE



#' #######################################################################
#' SIMPLE CORRELATION ANALYSIS
#' #######################################################################

dat %>%
  # filter to include dates only up to VAP (or last date, if no VAP)
  group_by(subject_id) %>%
  mutate(first_suspected_vap = ifelse(is.na(first_suspected_vap), max(date), first_suspected_vap)) %>%
  filter(date <= first_suspected_vap) %>%
  ungroup() %>%
  select(subject_id, subject_day, shannon, max_read_prop, log_copy_16S_per_ml_sputum, dominant_genus) %>%
  mutate(dominant_genus = replace(dominant_genus, is.na(dominant_genus), "no dominant genus")) %>%
  mutate_at(.vars = vars(-subject_id, -subject_day, -dominant_genus), .funs = list("scaled" = ~ as.vector(scale(.x)))) %>%
  distinct() %>%
  filter(complete.cases(.)) %>%
  identity() -> d_vap

d_vap



#' correlation between Shannon diversity and total bacterial abundance
d_vap %>%
  brm(data = ., 
      family = gaussian,
      formula = bf(mvbind(shannon_scaled, log_copy_16S_per_ml_sputum_scaled) ~ 0) + set_rescor(TRUE),
      prior = c(prior(lkj(2), class = rescor)
      ),
      iter = 2000,
      warmup = 1000,
      chains = 4,
      cores = 4,
      control = list("adapt_delta" = 0.99, max_treedepth = 16),
      backend = "cmdstanr",
      seed = 16
      ) -> m_cor_shannon_log16S

m_cor_shannon_log16S %>%
  brms::posterior_summary() %>%
  as_tibble(rownames = "param") %>%
  gt::gt() %>%
  gt::fmt_number(columns = c(-param), n_sigfig = 3)
# correlation: 0.0191	(95%CI −0.0617 to 0.0989)

m_cor_shannon_log16S %>%
  write_rds("./models/correlation/m_cor_shannon_log16S.rds.gz")



#' correlation between max proportional abundance and total bacterial abundance
d_vap %>%
  brm(data = ., 
      family = gaussian,
      formula = bf(mvbind(max_read_prop_scaled, log_copy_16S_per_ml_sputum_scaled) ~ 0) + set_rescor(TRUE),
      prior = c(prior(lkj(2), class = rescor)
      ),
      iter = 2000,
      warmup = 1000,
      chains = 4,
      cores = 4,
      control = list("adapt_delta" = 0.99, max_treedepth = 16),
      backend = "cmdstanr",
      seed = 16
  ) -> m_cor_max_read_prop_log16S

m_cor_max_read_prop_log16S %>%
  brms::posterior_summary() %>%
  as_tibble(rownames = "param") %>%
  gt::gt() %>%
  gt::fmt_number(columns = c(-param), n_sigfig = 3)
# correlation: 0.0358 (95%CI −0.0438 to 0.117)

m_cor_max_read_prop_log16S %>%
  write_rds("./models/correlation/m_cor_max_read_prop_log16S.rds.gz")




#' correlation between max proportional abundance and Shannon diversity
d_vap %>%
  brm(data = ., 
      family = gaussian,
      formula = bf(mvbind(shannon, max_read_prop_scaled) ~ 0) + set_rescor(TRUE),
      prior = c(prior(lkj(2), class = rescor)
      ),
      iter = 2000,
      warmup = 1000,
      chains = 4,
      cores = 4,
      control = list("adapt_delta" = 0.99, max_treedepth = 16),
      backend = "cmdstanr",
      seed = 16
  ) -> m_cor_shannon_max_read_prop

m_cor_shannon_max_read_prop %>%
  brms::posterior_summary() %>%
  as_tibble(rownames = "param") %>%
  gt::gt() %>%
  gt::fmt_number(columns = c(-param), n_sigfig = 3)
# correlation: −0.483 (95%CI −0.539 to −0.422)

m_cor_shannon_max_read_prop %>%
  write_rds("./models/correlation/m_cor_shannon_max_read_prop.rds.gz")






#' correlation between diversity and total bacterial abundance -- stratified by genus
d_vap %>%
  add_count(dominant_genus) %>%
  filter(n > 5) %>%
  group_by(dominant_genus) %>%
  nest() %>%
  mutate(brmod = map(.x = data, .f = ~ brm(data = .x,
                                           family = gaussian,
                                           formula = bf(mvbind(shannon_scaled, log_copy_16S_per_ml_sputum_scaled) ~ 0) + set_rescor(TRUE),
                                           prior = c(prior(lkj(2), class = rescor)
                                                     ),
                                           iter = 2000,
                                           warmup = 1000,
                                           chains = 4,
                                           cores = 4,
                                           control = list("adapt_delta" = 0.99, max_treedepth = 16),
                                           backend = "cmdstanr",
                                           seed = 16
                                           )
                     )
  ) %>%
  mutate(postsum = map(.x = brmod, .f = ~ brms::posterior_summary(.x) %>%
                         as_tibble(rownames = "param"))) %>%
  identity() -> m_cor_shannon_log16S_by_genus


m_cor_shannon_log16S_by_genus %>%
  mutate(rescor_post_med = map_dbl(.x = postsum, .f = ~ filter(.x, grepl("rescor__", param)) %>%
                                     pull(Estimate)),
         rescor_Q2.5 = map_dbl(.x = postsum, .f = ~ filter(.x, grepl("rescor__", param)) %>%
                                     pull(Q2.5)),
         rescor_Q97.5 = map_dbl(.x = postsum, .f = ~ filter(.x, grepl("rescor__", param)) %>%
                                     pull(Q97.5))) %>%
  select(dominant_genus, contains("rescor")) %>%
  ungroup() %>%
  arrange(desc(abs(rescor_post_med))) %>%
  write_csv("./models/correlation/m_cor_shannon_log16S_by_genus.csv") %>%
  gt::gt() %>%
  gt::fmt_number(columns = c(-dominant_genus), n_sigfig = 3) %>%
  gt::tab_header(title = "Correlate Shannon & log 16S Copies") %>%
  tab_style(
    style = list(
      cell_fill(color = "lightcyan")#,
      #cell_text(weight = "bold")
    ),
    locations = cells_body(
      columns = gt::everything(),
      rows = (rescor_post_med > 0 & rescor_Q2.5 > 0 & rescor_Q97.5 > 0) | (rescor_post_med < 0 & rescor_Q2.5 < 0 & rescor_Q97.5 < 0)
    )
  )




#' correlation between max proportional abundance and total bacterial abundance -- stratified by genus
d_vap %>%
  add_count(dominant_genus) %>%
  filter(n > 5) %>%
  group_by(dominant_genus) %>%
  nest() %>%
  mutate(brmod = map(.x = data, .f = ~ brm(data = .x,
                                           family = gaussian,
                                           formula = bf(mvbind(max_read_prop_scaled, log_copy_16S_per_ml_sputum_scaled) ~ 0) + set_rescor(TRUE),
                                           prior = c(prior(lkj(2), class = rescor)
                                           ),
                                           iter = 2000,
                                           warmup = 1000,
                                           chains = 4,
                                           cores = 4,
                                           control = list("adapt_delta" = 0.99, max_treedepth = 16),
                                           backend = "cmdstanr",
                                           seed = 16
  )
  )
  ) %>%
  mutate(postsum = map(.x = brmod, .f = ~ brms::posterior_summary(.x) %>%
                         as_tibble(rownames = "param"))) %>%
  identity() -> m_cor_max_read_prop_log16S_by_genus


m_cor_max_read_prop_log16S_by_genus %>%
  mutate(rescor_post_med = map_dbl(.x = postsum, .f = ~ filter(.x, grepl("rescor__", param)) %>%
                                     pull(Estimate)),
         rescor_Q2.5 = map_dbl(.x = postsum, .f = ~ filter(.x, grepl("rescor__", param)) %>%
                                 pull(Q2.5)),
         rescor_Q97.5 = map_dbl(.x = postsum, .f = ~ filter(.x, grepl("rescor__", param)) %>%
                                  pull(Q97.5))) %>%
  select(dominant_genus, contains("rescor")) %>%
  ungroup() %>%
  arrange(desc(abs(rescor_post_med))) %>%
  write_csv("./models/correlation/m_cor_max_read_prop_log16S_by_genus.csv") %>%
  gt::gt() %>%
  gt::fmt_number(columns = c(-dominant_genus), n_sigfig = 3) %>%
  gt::tab_header(title = "Correlate Maximum ASV Abundance & log 16S Copies") %>%
  tab_style(
    style = list(
      cell_fill(color = "lightcyan")#,
      #cell_text(weight = "bold")
    ),
    locations = cells_body(
      columns = gt::everything(),
      rows = (rescor_post_med > 0 & rescor_Q2.5 > 0 & rescor_Q97.5 > 0) | (rescor_post_med < 0 & rescor_Q2.5 < 0 & rescor_Q97.5 < 0)
    )
  )








#' #######################################################################
#' #######################################################################
#' VAP as a function of lower respiratory tract diversity + total abundance
#' #######################################################################
#' #######################################################################


#' #######################################################################
#' SUSPECTED VAP + CX + SIGNS (same as SUSPECTED VAP + CX)
#' #######################################################################

m_ref <- read_rds("./models/longitudinal/m_minimpute_svp_minimpute_shannon_low_count_brms.rds.gz")
add_criterion(m_ref, criterion = c("waic","loo"))
loo(m_ref)
waic(m_ref)



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
  group_by(subject_id) %>%
  arrange(subject_day) %>%
  mutate(shannon_less_than_2 = shannon < 2,
         shannon_improved = shannon_less_than_2 == FALSE & lag(shannon_less_than_2, 1) == TRUE) %>%
  replace_na(replace = list("shannon_improved" = FALSE)) %>%
  mutate(shannon_flip_count = as.character(cumsum(shannon_improved))) %>%
  group_by(subject_id, shannon_flip_count) %>%
  arrange(subject_day) %>%
  mutate(shannon_low_count = cumsum(shannon_less_than_2)) %>%
  ungroup() %>%
  filter(complete.cases(.)) %>%
  identity() -> d_vap

d_vap



#' brms
d_vap %>%
  select(suspected_vap_plus_signs_and_pathogen, subject_id, shannon_low_count, shannon_scaled, log_copy_16S_per_ml_sputum_scaled) %>%
  filter(complete.cases(.)) %>%
  brm(data = ., family = bernoulli,
      #suspected_vap_plus_signs_and_pathogen ~ (1 | subject_id) + shannon_scaled,
      suspected_vap_plus_signs_and_pathogen ~ 1 + shannon_scaled,
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
        seed = 16) -> m1_suspected_vap_plus_signs_and_pathogen_shannon_only


d_vap %>%
  select(suspected_vap_plus_signs_and_pathogen, subject_id, shannon_low_count, shannon_scaled, log_copy_16S_per_ml_sputum_scaled) %>%
  filter(complete.cases(.)) %>%
  brm(data = ., family = bernoulli,
      #suspected_vap_plus_signs_and_pathogen ~ (1 | subject_id) + log_copy_16S_per_ml_sputum_scaled,
      suspected_vap_plus_signs_and_pathogen ~ 1 + log_copy_16S_per_ml_sputum_scaled,
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
      seed = 16) -> m2_suspected_vap_plus_signs_and_pathogen_log16S_only


d_vap %>%
  select(suspected_vap_plus_signs_and_pathogen, subject_id, shannon_low_count, shannon_scaled, log_copy_16S_per_ml_sputum_scaled) %>%
  filter(complete.cases(.)) %>%
  brm(data = ., family = bernoulli,
      #suspected_vap_plus_signs_and_pathogen ~ (1 | subject_id) + shannon_scaled + log_copy_16S_per_ml_sputum_scaled,
      suspected_vap_plus_signs_and_pathogen ~ 1 + shannon_scaled + log_copy_16S_per_ml_sputum_scaled,
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
      seed = 16) -> m3_suspected_vap_plus_signs_and_pathogen_shannon_log16S_nointerax


d_vap %>%
  select(suspected_vap_plus_signs_and_pathogen, subject_id, shannon_low_count, shannon_scaled, log_copy_16S_per_ml_sputum_scaled) %>%
  filter(complete.cases(.)) %>%
  brm(data = ., family = bernoulli,
      #suspected_vap_plus_signs_and_pathogen ~ (1 | subject_id) + shannon_scaled * log_copy_16S_per_ml_sputum_scaled,
      suspected_vap_plus_signs_and_pathogen ~ 1 + shannon_scaled * log_copy_16S_per_ml_sputum_scaled,
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
      seed = 16) -> m4_suspected_vap_plus_signs_and_pathogen_shannon_log16S_interax


d_vap %>%
  select(suspected_vap_plus_signs_and_pathogen, subject_id, shannon_low_count, shannon_scaled, log_copy_16S_per_ml_sputum_scaled) %>%
  filter(complete.cases(.)) %>%
  brm(data = ., family = bernoulli,
      #suspected_vap_plus_signs_and_pathogen ~ (1 | subject_id) + shannon_low_count,
      suspected_vap_plus_signs_and_pathogen ~ 1 + shannon_low_count,
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
      seed = 16) -> m5_ref



m1 <- add_criterion(m1_suspected_vap_plus_signs_and_pathogen_shannon_only, criterion = c("loo","waic"))
m2 <- add_criterion(m2_suspected_vap_plus_signs_and_pathogen_log16S_only, criterion = c("loo","waic"))
m3 <- add_criterion(m3_suspected_vap_plus_signs_and_pathogen_shannon_log16S_nointerax, criterion = c("loo","waic"))
m4 <- add_criterion(m4_suspected_vap_plus_signs_and_pathogen_shannon_log16S_interax, criterion = c("loo","waic"))
m5 <- add_criterion(m5_ref, criterion = c("loo","waic"))

loo(m1,m2,m3,m4,m5)
waic(m1,m2,m3,m4,m5)
loo_compare(m1, m2, m3, m4, m5) %>%
  as_tibble(rownames = "model")






