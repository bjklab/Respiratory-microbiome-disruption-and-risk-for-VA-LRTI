#' ########################################
#' make Table 1 to display baseline (day 0) characteristics
#' ########################################
#' 
#' depends: micu_cohort
#' 

library(tidyverse)
library(gt)
library(gtsummary)
library(brms)
set.seed(16)


#' table 1 for MIRA cohort
subj %>%
  select(admission_shannon_above_median,
         admission_shannon,
         sputum_sampling_duration,
         age,
         gender,
         race,
         contains("admit"),
         days_mechan_vent_prior_enroll,
         days_trach_prior_enroll,
         prior_pneumonia_noted,
         copd,
         asthma,
         ild,
         lymphleuk,
         dm,
         chf,
         cirrhosis,
         contains('vanco_before0_7'),
         contains('metro_before0_7'),
         contains('linez_before0_7'),
         contains('dapto_before0_7'),
         contains('cefaz_before0_7'),
         contains('piptaz_before0_7'),
         contains('cefep_before0_7'),
         contains('mero_before0_7')) %>%
  mutate(race = stringr::str_to_title(race),
         gender = case_when(gender == "F" ~ "Female",
                            gender == "M" ~ "Male")) %>%
  mutate(admission_shannon_above_median = ifelse(admission_shannon_above_median, "High Admission Shannon Diversity", "Low Admission Shannon Diversity")) %>%
  rename_all(.funs = ~ stringr::str_to_title(stringr::str_replace_all(.x,"_"," "))) %>%
  rename_all(.funs = ~ stringr::str_remove(string = .x, pattern = "Admit | Before0 7")) %>%
  rename_at(.vars = vars(Wbc,Ast,Alt,Fio2,Peep,Copd,Ild,Dm,Chf), .funs = ~ stringr::str_to_upper(.x)) %>%
  rename_all(.funs = ~ stringr::str_replace_all(string = .x, pattern = "Pds", "PDS")) %>%
  rename_all(.funs = ~ stringr::str_replace_all(string = .x, pattern = "FIO2", "FiO2 (%)")) %>%
  rename_all(.funs = ~ stringr::str_replace_all(string = .x, pattern = "Sputum Sampling Duration", "Sampling Duration (days)")) %>%
  rename_all(.funs = ~ stringr::str_replace_all(string = .x, pattern = "Lymphleuk", "Lymphoma/Leukemia")) %>%
  rename_all(.funs = ~ stringr::str_replace_all(string = .x, pattern = "Admission Shannon$", "Admission Shannon Diversity (base e)")) %>%
  rename_all(.funs = ~ stringr::str_replace_all(string = .x, pattern = "Vanco", "Vancomycin (IV)")) %>%
  rename_all(.funs = ~ stringr::str_replace_all(string = .x, pattern = "Metro", "Metronidazole")) %>%
  rename_all(.funs = ~ stringr::str_replace_all(string = .x, pattern = "Linez", "Linezolid")) %>%
  rename_all(.funs = ~ stringr::str_replace_all(string = .x, pattern = "Dapto", "Daptomycin")) %>%
  rename_all(.funs = ~ stringr::str_replace_all(string = .x, pattern = "Cefaz", "Cefazolin")) %>%
  rename_all(.funs = ~ stringr::str_replace_all(string = .x, pattern = "Piptaz", "Piperacillin-tazobactam")) %>%
  rename_all(.funs = ~ stringr::str_replace_all(string = .x, pattern = "Cefep", "Cefepime")) %>%
  rename_all(.funs = ~ stringr::str_replace_all(string = .x, pattern = "Mero", "Meropenem")) %>%
  rename_all(.funs = ~ stringr::str_replace_all(string = .x, pattern = "Mechan Vent", "Mechanical Ventilation")) %>%
  rename_all(.funs = ~ stringr::str_replace_all(string = .x, pattern = "Trach", "Tracheostomy")) %>%
  rename_all(.funs = ~ stringr::str_replace_all(string = .x, pattern = "Prior Enroll", "Prior to LTACH Admission")) %>%
  rename_all(.funs = ~ stringr::str_replace_all(string = .x, pattern = "Prior Pneumonia Noted", "Pneumonia Noted Prior to LTACH Admission")) %>%
  tbl_summary(
    data = .,
    by = `Admission Shannon Above Median`,
    missing = "no",
    type = list("Admission Shannon Diversity (base e)" ~ "continuous", "Sampling Duration (days)" ~ "continuous", "FiO2 (%)" ~ "continuous", "PEEP" ~ "continuous")
  ) %>%
  #add_n() %>%
  add_overall() %>%
  add_p(pvalue_fun = function(x) style_pvalue(x, digits = 2)) %>%
  gtsummary::modify_header(update = list(label ~ "**Subject Characteristics**", stat_0 ~ "**All Subjects**, N = {N}")) %>%
  gtsummary::modify_spanning_header(c("stat_0", "stat_1", "stat_2") ~ "**Subject Category & Number**") %>%
  gtsummary::as_gt() %>%
  gt::tab_row_group(group = "Primary Exposure and Enrollment Duration", rows = c(1:2) ) %>%
  gt::tab_row_group(group = "Demographics and Medical Comorbidities", rows = c(3:12,19:28) ) %>%
  gt::tab_row_group(group = "Admission Laboratory Values and Ventilator Settings", rows = 13:18 ) %>%
  #gt::tab_row_group(group = "Medical Comorbidities", rows = 19:25 ) %>%
  gt::tab_row_group(group = "Antibiotics Within 7 Days Prior to Admission", rows = 29:36 ) %>%
  gt::row_group_order(groups = c("Primary Exposure and Enrollment Duration", "Demographics and Medical Comorbidities", "Admission Laboratory Values and Ventilator Settings", "Antibiotics Within 7 Days Prior to Admission")) %>%
  #gt::tab_header(title = "MIRA Cohort: High vs Low Admission Respiratory Microbiome Diversity", subtitle = "") %>%
  identity() -> t1_shannon

t1_shannon

# t1_shannon %>%
#   gt::as_raw_html() %>%
#   write_lines(file = "./tabs/t1_shannon.html")
# 
# t1_shannon %>%
#   gt::as_rtf() %>%
#   write_lines(file = "./tabs/t1_shannon.rtf")


subj %>%
  replace_na(replace = list("prior_icu_indication" = "prior records not available")) %>%
  mutate(prior_icu_indication = replace(prior_icu_indication, grepl("stroke|subarachnoid|neuromuscular|ALS",prior_icu_indication), "neurologic disorder")) %>%
  mutate(prior_icu_indication = replace(prior_icu_indication, grepl("cancer|carcinoma",prior_icu_indication), "cancer")) %>%
  count(prior_icu_indication) %>%
  arrange(desc(n), prior_icu_indication) %>%
  rename(`Indication for Critical Care Prior to LTACH Admission` = prior_icu_indication,
         `Subject Count` = n) %>%
  gt::gt() %>%
  identity() -> supp_t1

supp_t1

# supp_t1 %>%
#   gt::as_raw_html() %>%
#   write_lines(file = "./tabs/supp_t1.html")
# 
# supp_t1 %>%
#   gt::as_rtf() %>%
#   write_lines(file = "./tabs/supp_t1.rtf")


#' #' per reviewer request, interrogate admission shannon vs prior vent days
#' subj %>%
#'   select(admission_shannon, days_mechan_vent_prior_enroll, days_trach_prior_enroll) %>%
#'   filter(complete.cases(.)) %>%
#'   brms::brm(data = .,
#'             family = "gaussian",
#'             formula = admission_shannon ~ days_mechan_vent_prior_enroll + 1,
#'             seed = 16,
#'             cores = 4,
#'             iter = 2000,
#'             backend = "cmdstanr") %>%
#'   identity() -> m_admit_shannon_prior_vent_days
#' 
#' m_admit_shannon_prior_vent_days %>%
#'   brms::posterior_summary() %>%
#'   as_tibble(rownames = "param")
#' 
#' 
#' #' per reviewer request, interrogate admission shannon vs prior trach days
#' subj %>%
#'   select(admission_shannon, days_mechan_vent_prior_enroll, days_trach_prior_enroll) %>%
#'   filter(complete.cases(.)) %>%
#'   brms::brm(data = .,
#'             family = "gaussian",
#'             formula = admission_shannon ~ days_trach_prior_enroll + 1,
#'             seed = 16,
#'             cores = 4,
#'             iter = 2000,
#'             backend = "cmdstanr") %>%
#'   identity() -> m_admit_shannon_prior_trach_days
#' 
#' m_admit_shannon_prior_trach_days %>%
#'   brms::posterior_summary() %>%
#'   as_tibble(rownames = "param")


#
###
#####
###
#







