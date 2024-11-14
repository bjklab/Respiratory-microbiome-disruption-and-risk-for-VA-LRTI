#' ###########################################
#' MIXED EFFECT VAP ~ LONGITUDINAL MDI models
#' ###########################################
library(tidyverse)
library(tidybayes)
library(brms)
set.seed(16)



#' #######################################################################
#' #######################################################################
#' READ MIXED EFFECTS MODELS - OUTCOME = "suspected_vap_plus_pathogen"
#' #######################################################################
#' #######################################################################


#' #######################################################################
#' SHANNON DIVERSITY
#' #######################################################################

structure(.Data = list.files(path = "./models/", pattern = "shannon_brms_mixed_restricted.rds", full.names = TRUE),
          .Names = gsub("\\.rds\\.gz","", list.files(path = "./models/", pattern = "shannon_brms_mixed_restricted.rds"))) %>%
  map(.f = ~ read_rds(.x)) %>%
  identity() -> shannon_models

names(shannon_models)

shannon_models %>%
  setNames(object = ., nm = dplyr::case_when(grepl("suspected_vap_shannon", names(.)) ~ "LRTI = Respiratory Culture Ordered",
                                             grepl("vap_plus_pathogen_shannon", names(.)) ~ "LRTI = Positive Respiratory Culture",
                                             grepl("plus_abx_and_pathogen_shannon", names(.)) ~ "LRTI = Positive Culture + Antibiotic Change",
                                             grepl("vap_plus_abx_and_pathogen_and_signs", names(.)) ~ "LRTI = Positive Culture + Inflammation + Antibiotic Change",
  )) %>%
  identity() -> shannon_models
names(shannon_models)





#' #######################################################################
#' QPCR ABSOLUTE BACTERIAL ABUNDANCE
#' #######################################################################

structure(.Data = list.files(path = "./models/", pattern = "log_copy_16S_per_ml_sputum_brms_mixed_restricted.rds", full.names = TRUE),
          .Names = gsub("\\.rds\\.gz","", list.files(path = "./models/", pattern = "log_copy_16S_per_ml_sputum_brms_mixed_restricted.rds"))) %>%
  map(.f = ~ read_rds(.x)) %>%
  identity() -> qpcr_models

names(qpcr_models)

qpcr_models %>%
  setNames(object = ., nm = dplyr::case_when(grepl("suspected_vap_log_copy", names(.)) ~ "LRTI = Respiratory Culture Ordered",
                                             grepl("vap_plus_pathogen_log_copy", names(.)) ~ "LRTI = Positive Respiratory Culture",
                                             grepl("plus_abx_and_pathogen_log_copy", names(.)) ~ "LRTI = Positive Culture + Antibiotic Change",
                                             grepl("vap_plus_abx_and_pathogen_and_signs", names(.)) ~ "LRTI = Positive Culture + Inflammation + Antibiotic Change",
  )) %>%
  identity() -> qpcr_models
names(qpcr_models)




#' #######################################################################
#' MAX PROPORTIONAL BACTERIAL ABUNDANCE
#' #######################################################################

structure(.Data = list.files(path = "./models/", pattern = "max_read_prop_brms_mixed_restricted.rds", full.names = TRUE),
          .Names = gsub("\\.rds\\.gz","", list.files(path = "./models/", pattern = "max_read_prop_brms_mixed_restricted.rds"))) %>%
  map(.f = ~ read_rds(.x)) %>%
  identity() -> maxprop_models

names(maxprop_models)

maxprop_models %>%
  setNames(object = ., nm = dplyr::case_when(grepl("suspected_vap_max_read_prop", names(.)) ~ "LRTI = Respiratory Culture Ordered",
                                             grepl("vap_plus_pathogen_max_read_prop", names(.)) ~ "LRTI = Positive Respiratory Culture",
                                             grepl("plus_abx_and_pathogen_max_read_prop", names(.)) ~ "LRTI = Positive Culture + Antibiotic Change",
                                             grepl("vap_plus_abx_and_pathogen_and_signs", names(.)) ~ "LRTI = Positive Culture + Inflammation + Antibiotic Change",
  )) %>%
  identity() -> maxprop_models
names(maxprop_models)




#' #######################################################################
#' SCALE DATA
#' #######################################################################

dat_r %>%
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
  select(subject_id, subject_day, shannon, max_read_prop, log_copy_16S_per_ml_sputum, suspected_vap_plus_pathogen) %>%
  mutate_at(.vars = vars(-subject_id, -subject_day, -suspected_vap_plus_pathogen), .funs = list("scaled" = ~ as.vector(scale(.x)))) %>%
  distinct() %>%
  filter(complete.cases(.)) %>%
  identity() -> d_vap

d_vap


d_vap %>%
  select(shannon, log_copy_16S_per_ml_sputum, max_read_prop) %>%
  summarise_all(.funs = list(mean = ~ mean(.x, na.rm = TRUE), sd = ~ sd(.x, na.rm = TRUE))) %>%
  identity() -> d_scale
d_scale


#' #######################################################################
#' #######################################################################
#' COUNTERFACTUAL PLOTS
#' #######################################################################
#' #######################################################################

#' functions
counterfact <- function(brms_mod, exp_var_name = "shannon_low_count", draw_number = NULL) {
  brms_mod[["data"]] %>%
    expand(exp_var = modelr::seq_range(brms_mod[["data"]][[exp_var_name]], n = 100),
           subject_id = unique(subject_id)) %>%
    rename_with(.cols = exp_var, .fn = ~ exp_var_name) %>%
    tidybayes::add_fitted_draws(model = brms_mod, re_formula = NULL, n = draw_number, seed = 16) %>% # n = NULL for all draws
    identity() -> post_fitted
  return(post_fitted)
}



#' shannon

map_dfr(shannon_models,
        .f = ~ counterfact(brms_mod = .x,
                           exp_var_name = "shannon_scaled",
                           draw_number = 100),
        .id = "model") %>%
  ungroup() %>%
  filter(model == "LRTI = Positive Respiratory Culture") %>%
  mutate(exp_var_name = "shannon_scaled") %>%
  rename(exp_var_scaled = shannon_scaled) %>%
  mutate(exp_var_raw = exp_var_scaled * d_scale$shannon_sd + d_scale$shannon_mean) %>%
  identity() -> shannon_cf
shannon_cf




#' qpcr

map_dfr(qpcr_models,
        .f = ~ counterfact(brms_mod = .x,
                           exp_var_name = "log_copy_16S_per_ml_sputum_scaled",
                           draw_number = 100),
        .id = "model") %>%
  ungroup() %>%
  filter(model == "LRTI = Positive Respiratory Culture") %>%
  mutate(exp_var_name = "log_copy_16S_per_ml_sputum_scaled") %>%
  rename(exp_var_scaled = log_copy_16S_per_ml_sputum_scaled) %>%
  mutate(exp_var_raw = exp_var_scaled * d_scale$log_copy_16S_per_ml_sputum_sd + d_scale$log_copy_16S_per_ml_sputum_mean) %>%
  identity() -> qpcr_cf
qpcr_cf



#' maxprop

map_dfr(maxprop_models,
        .f = ~ counterfact(brms_mod = .x,
                           exp_var_name = "max_read_prop_scaled",
                           draw_number = 100),
        .id = "model") %>%
  ungroup() %>%
  filter(model == "LRTI = Positive Respiratory Culture") %>%
  mutate(exp_var_name = "max_read_prop_scaled") %>%
  rename(exp_var_scaled = max_read_prop_scaled) %>%
  mutate(exp_var_raw = exp_var_scaled * d_scale$max_read_prop_sd + d_scale$max_read_prop_mean) %>%
  identity() -> maxprop_cf
maxprop_cf



#' combine counterfactual data
bind_rows(shannon_cf,
          qpcr_cf,
          maxprop_cf) %>%
  mutate(exp_var_label = case_when(grepl("shannon",exp_var_name) ~ "Shannon Diversity (base e)",
                                   grepl("log_copy", exp_var_name) ~ "Total Bacterial Abundance<br>by 16S rRNA Gene qPCR<br>(log copies per mL sputum)",
                                   grepl("max_read_prop", exp_var_name) ~ "Maximum ASV<br>Proportional Abundance")) %>%
  mutate(exp_var_label = factor(exp_var_label, levels = c("Shannon Diversity (base e)",
                                                          "Maximum ASV<br>Proportional Abundance",
                                                          "Total Bacterial Abundance<br>by 16S rRNA Gene qPCR<br>(log copies per mL sputum)"))) %>%
  identity() -> combined_mdi_cf
combined_mdi_cf



#' plot combined data
cf_plot <- function(brms_counterfact_output,
                    exp_var_name = "exp_var_raw",
                    exp_lab = "",
                    out_lab = "Probability of VA-LRTI",
                    facet_var_name = "exp_var_label",
                    facet_scale_string = "free") {
  brms_counterfact_output %>%
    ggplot(data = ., aes_string(x = exp_var_name, y = ".value")) +
    tidybayes::stat_lineribbon(.width = c(0.5,0.8,0.95),
                               alpha = 0.7,
                               color = colorspace::sequential_hcl(5, palette = "Blues 3")[1]) +
    colorspace::scale_fill_discrete_sequential(palette = "Blues 3", nmax = 5, order = 2:4) +
    facet_wrap(eval(expr( ~ !!ensym(facet_var_name))), scales = facet_scale_string) +
    theme_bw() +
    theme(legend.position = c(0.25,0.82),
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
    labs(x = exp_lab, y = out_lab, fill = "Posterior<br>Credible<br>Interval") -> p
  return(p)
}



combined_mdi_cf %>%
  cf_plot(brms_counterfact_output = .,
          exp_var_name = "exp_var_raw",
          exp_lab = "Microbiome Disruption Index",
          out_lab = "Probability of VA-LRTI",
          facet_var_name = "exp_var_label",
          facet_scale_string = "free_x") -> p_combined_restricted_mdi_crosssection_counterfact
p_combined_restricted_mdi_crosssection_counterfact



p_combined_restricted_mdi_crosssection_counterfact +
  tidybayes::stat_lineribbon(.width = c(0.5,0.8,0.95),
                             alpha = 0.7,
                             color = colorspace::sequential_hcl(5, palette = "Light Grays")[1]) +
  colorspace::scale_fill_discrete_sequential(palette = "Light Grays", nmax = 5, order = 2:4)



# p_combined_restricted_mdi_crosssection_counterfact %>%
#   ggsave(plot = ., filename = "./figs/p_combined_restricted_mdi_crosssection_counterfact.pdf", height = 5, width = 6, units = "in")
# p_combined_restricted_mdi_crosssection_counterfact %>%
#   ggsave(plot = ., filename = "./figs/p_combined_restricted_mdi_crosssection_counterfact.svg", height = 5, width = 6, units = "in")
# p_combined_restricted_mdi_crosssection_counterfact %>%
#   ggsave(plot = ., filename = "./figs/p_combined_restricted_mdi_crosssection_counterfact.png", height = 5, width = 6, units = "in", dpi = 600)






