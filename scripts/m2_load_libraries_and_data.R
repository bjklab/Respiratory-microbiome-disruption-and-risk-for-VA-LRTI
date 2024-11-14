#' #############################################
#' load libraries and set seed
#' #############################################
library(tidyverse)
library(tidybayes)
library(brms)
set.seed(16)



#' #############################################
#' load data
#' #############################################
dat <- read_csv("./data/d_mdi_outcome.csv") %>%
  mutate(copy_number_per_ml_sputum = as.numeric(replace(copy_number_per_ml_sputum, copy_number_per_ml_sputum == "Undetermined", NA))) %>%
  mutate(log_copy_16S_per_ml_sputum = log10(as.numeric(copy_number_per_ml_sputum)))
dat


# dat_r <- read_csv("./data/d_mdi_outcome_restricted.csv") %>%
#   mutate(copy_number_per_ml_sputum = as.numeric(replace(copy_number_per_ml_sputum, copy_number_per_ml_sputum == "Undetermined", NA))) %>%
#   mutate(log_copy_16S_per_ml_sputum = log10(as.numeric(copy_number_per_ml_sputum)))
# dat_r


subj <- read_csv("./data/mira_subject_summary.csv") %>%
  left_join(distinct(select(filter(dat, subject_day == 0), subject_id, shannon)), by = "subject_id") %>%
  rename(admission_shannon = shannon) %>%
  mutate(admission_shannon_above_median = admission_shannon > median(admission_shannon, na.rm = TRUE)) %>%
  left_join(read_csv("./data/mira_subject_prior.csv"), by = "subject_id") %>%
  mutate(prior_pneumonia_noted = grepl("pneumonia", prior_icu_indication))
subj




#' #############################################
#' impute data
#' #############################################
dat %>%
  select(subject_id, date, subject_day, shannon, log_copy_16S_per_ml_sputum, max_read_prop, dominant_genus, contains("suspected")) %>%
  {. ->> d_vap_draft} %>%
  #
  # impute missing subject days
  # imputation = last recorded value
  #
  group_by(subject_id) %>%
  expand(subject_day = full_seq(subject_day, 1)) %>%
  ungroup() %>%
  left_join(d_vap_draft, by = c("subject_id", "subject_day")) %>%
  group_by(subject_id) %>%
  mutate_at(.vars = vars(contains("^suspected")), .funs = ~ replace(.x, is.na(.x), FALSE)) %>%
  mutate_at(.vars = vars(contains("^first")), .funs = ~ replace(.x, is.na(.x), unique(.x[!is.na(.x)]))) %>%
  arrange(subject_id, subject_day) %>%
  mutate_all(.funs = ~ replace(.x, is.nan(.x), NA)) %>%
  fill(shannon, log_copy_16S_per_ml_sputum, max_read_prop, dominant_genus, .direction = "down") %>%
  ungroup() %>%
  #
  # check imputation
  #View()
  #count(subject_id) %>% print(n=Inf)
  identity() -> dat_i
dat_i



# dat_r %>%
#   select(subject_id, date, subject_day, shannon, log_copy_16S_per_ml_sputum, max_read_prop, dominant_genus, contains("suspected")) %>%
#   {. ->> d_vap_draft} %>%
#   #
#   # impute missing subject days
#   # imputation = last recorded value
#   #
#   group_by(subject_id) %>%
#   expand(subject_day = full_seq(subject_day, 1)) %>%
#   ungroup() %>%
#   left_join(d_vap_draft, by = c("subject_id", "subject_day")) %>%
#   group_by(subject_id) %>%
#   mutate_at(.vars = vars(contains("^suspected")), .funs = ~ replace(.x, is.na(.x), FALSE)) %>%
#   mutate_at(.vars = vars(contains("^first")), .funs = ~ replace(.x, is.na(.x), unique(.x[!is.na(.x)]))) %>%
#   arrange(subject_id, subject_day) %>%
#   mutate_all(.funs = ~ replace(.x, is.nan(.x), NA)) %>%
#   fill(shannon, log_copy_16S_per_ml_sputum, max_read_prop, dominant_genus, .direction = "down") %>%
#   ungroup() %>%
#   #
#   # check imputation
#   #View()
#   #count(subject_id) %>% print(n=Inf)
#   identity() -> dat_ri
# dat_ri


# subj %>%
#   filter(subject_id %in% unique(dat_r$subject_id)) %>%
#   identity() -> subj_r
# subj_r





