#' ###########################################
#' VA-LRTI OUTCOME COMPARISON
#' ###########################################


#' #######################################################################
#' #######################################################################
#' load libraries and set seed
#' #######################################################################
#' #######################################################################
library(tidyverse)
library(tidybayes)
library(brms)
set.seed(16)



#' #######################################################################
#' #######################################################################
#' VA-LRTI Definition Comparison
#' #######################################################################
#' #######################################################################

#' compare definitions
dat %>%
  select(subject_id, subject_day, contains("suspected")) %>%
  group_by(subject_id) %>%
  summarise(any_suspected_vap = sum(suspected_vap, na.rm = TRUE) > 0,
            suspected_vap_day = ifelse(any_suspected_vap, min(subject_day[suspected_vap == TRUE], na.rm = TRUE), max(subject_day, na.rm = TRUE)),
            any_vap_cx = sum(suspected_vap_plus_pathogen, na.rm = TRUE) > 0,
            vap_cx_day = ifelse(any_vap_cx, min(subject_day[suspected_vap_plus_pathogen == TRUE], na.rm = TRUE), max(subject_day, na.rm = TRUE)),
            any_vap_cx_signs = sum(suspected_vap_plus_signs_and_pathogen, na.rm = TRUE) > 0,
            vap_cx_signs_day = ifelse(any_vap_cx_signs, min(subject_day[suspected_vap_plus_signs_and_pathogen == TRUE], na.rm = TRUE), max(subject_day, na.rm = TRUE)),
            any_vap_cx_signs_abx = sum(suspected_vap_plus_abx_and_pathogen_and_signs, na.rm = TRUE) > 0,
            vap_cx_signs_abx_day = ifelse(any_vap_cx_signs_abx, min(subject_day[suspected_vap_plus_abx_and_pathogen_and_signs == TRUE], na.rm = TRUE), max(subject_day, na.rm = TRUE)),
            ) %>%
  ungroup() %>%
  identity() -> t_subject_vap
t_subject_vap


t_subject_vap %>%
  select(subject_id, contains("any")) %>%
  gather(key = "vap_definition", value = "vap_yn", -subject_id) %>%
  mutate(vap_definition = gsub("any_","",vap_definition)) %>%
  identity() -> t_vap_definition_yn
t_vap_definition_yn


#' count subjects 
t_subject_vap %>%
  select(subject_id, contains("any")) %>%
  gather(key = "vap_definition", value = "vap_yn", -subject_id) %>%
  mutate(vap_definition = gsub("any_","",vap_definition)) %>%
  group_by(vap_definition) %>%
  count(vap_yn) %>%
  rename(`VA-LRTI Occurred?` = vap_yn,
         `VA-LRTI Definition` = vap_definition) %>%
  ungroup() %>%
  gt::gt() %>%
  gt::tab_header(title = "Comparison of VA-LRTI Definitions")


#' time to VA-LRTI
t_subject_vap %>%
  select(subject_id, contains("day")) %>%
  gather(key = "vap_definition", value = "days", -subject_id) %>%
  mutate(vap_definition = gsub("_day","",vap_definition)) %>%
  left_join(t_vap_definition_yn, by = c("subject_id","vap_definition")) %>%
  group_by(vap_definition) %>%
  mutate(vap_tally = sum(vap_yn, na.rm = TRUE)) %>%
  mutate(vap_definition = paste0(vap_definition, " (",vap_tally," subjects)")) %>%
  ungroup() %>%
  group_by(vap_definition, vap_yn) %>%
  summarise_at(.vars = vars(days), .funs = list(mean = ~ mean(.x, na.rm = TRUE),
                                                          sd = ~ sd(.x, na.rm = TRUE),
                                                          median = ~ median(.x, na.rm = TRUE),
                                                          iqr = ~ IQR(.x, na.rm = TRUE))
               ) %>%
  rename(`VA-LRTI Occurred?` = vap_yn,
         `VA-LRTI Definition` = vap_definition,
         ) %>%
  ungroup() %>%
  gt::gt() %>%
  gt::tab_header(title = "Comparison of VA-LRTI Definitions") %>%
  gt::tab_spanner(columns = 3:6, label = "Days Until VA-LRTI (or study exit)") %>%
  gt::fmt_number(columns = 3:6, n_sigfig = 3)




#' #######################################################################
#' PATHOGEN PROFILES
#' #######################################################################

dat %>%
  # filter to include dates only up to VAP (or last date, if no VAP)
  group_by(subject_id) %>%
  mutate(first_suspected_vap_plus_pathogen = ifelse(is.na(first_suspected_vap_plus_pathogen), max(date), first_suspected_vap_plus_pathogen)) %>%
  filter(date == first_suspected_vap_plus_pathogen) %>%
  ungroup() %>%
  filter(suspected_vap_plus_pathogen == TRUE) %>%
  select(subject_id, date, respiratory_cx_pathogen) %>%
  distinct() %>%
  #add_count(subject_id) %>% arrange(desc(n)) %>% print(n = Inf)
  separate(col = respiratory_cx_pathogen, into = paste0("bacteria",seq(5)), sep = " ", remove = TRUE) %>%
  gather(key = "bacteria_count", value = "cultured_bacteria", -subject_id, -date) %>%
  filter(!is.na(cultured_bacteria)) %>%
  select(-bacteria_count) %>%
  mutate(cultured_bacteria = stringr::str_to_title(cultured_bacteria),
         cultured_bacteria = gsub("Staphylococcus$","Staphylococcus (non-aureus)",cultured_bacteria),
         cultured_bacteria = gsub("_"," ", cultured_bacteria)) %>%
  identity() -> t_respiratory_cx_pathogen
t_respiratory_cx_pathogen


t_respiratory_cx_pathogen %>%
  count(cultured_bacteria) %>%
  arrange(desc(n)) %>%
  rename(`Bacterial Genus` = cultured_bacteria,
         `Number of Isolates` = n) %>%
  gt::gt() %>%
  identity() -> supp_t2

supp_t2


# supp_t2 %>%
#   gt::as_raw_html() %>%
#   write_lines(file = "./tabs/supp_t2.html")
# 
# supp_t2 %>%
#   gt::as_rtf() %>%
#   write_lines(file = "./tabs/supp_t2.rtf")



t_respiratory_cx_pathogen %>%
  arrange(subject_id, date) %>%
  left_join(select(dat, subject_id, date, dominant_genus), by = c("subject_id","date")) %>%
  mutate(cultured_bacteria = gsub(" aureus","",cultured_bacteria),
         ) %>%
  group_by(subject_id, date) %>%
  # analysis restricted to positive cultures with same-day dominant ASV
  filter(!is.na(dominant_genus)) %>%
  summarise(any_match = sum(cultured_bacteria == dominant_genus, na.rm = TRUE) > 0) %>%
  ungroup()
# 6 of 9 dominant ASVs match same-day culture





