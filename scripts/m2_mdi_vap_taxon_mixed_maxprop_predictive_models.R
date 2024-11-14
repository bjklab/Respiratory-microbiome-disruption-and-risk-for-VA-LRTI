#' ###########################################
#' MIXED EFFECT VAP ~ MDI models
#' - depends on libraries and data
#' ###########################################


#' #######################################################################
#' #######################################################################
#' VAP as a function of lower respiratory tract max prop abundance
#' #######################################################################
#' #######################################################################
library(tidyverse)
library(tidybayes)
library(brms)
set.seed(16)
run_models = TRUE
save_models = TRUE
run_plots = TRUE




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
  gt::gt()





#' #######################################################################
#' SUSPECTED VAP + SIGNS + PATHOGEN
#' #######################################################################

if(run_models) {
dat %>%
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
  select(subject_id, subject_day, shannon, max_read_prop, log_copy_16S_per_ml_sputum, suspected_vap_plus_pathogen, dominant_genus, dominant_asv) %>%
  mutate_at(.vars = vars(-subject_id, -subject_day, -suspected_vap_plus_pathogen, -contains("dominant")), .funs = list("scaled" = ~ as.vector(scale(.x)))) %>%
  distinct() %>%
  filter(complete.cases(.)) %>%
  identity() -> d_vap

d_vap




d_vap %>%
  get_prior(data = ., family = bernoulli,
            suspected_vap_plus_pathogen ~ (1 | subject_id) + (1 + max_read_prop_scaled | dominant_asv)
  )


d_vap %>%
  brm(data = ., family = bernoulli,
      suspected_vap_plus_pathogen ~ (1 | subject_id) + (1 + max_read_prop_scaled | dominant_genus),
      prior = c(prior(student_t(3,0,2.5), class = Intercept)
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
#' genus
d_vap %>%
  brm(data = ., family = bernoulli,
      suspected_vap_plus_pathogen ~ (1 | subject_id) + (1 + max_read_prop_scaled | dominant_genus),
      prior = c(prior(student_t(3,0,2.5), class = Intercept)
      ),
      #sample_prior = "only", # prior predictive
      iter = 2000,
      warmup = 1000,
      chains = 4,
      cores = 4,
      control = list("adapt_delta" = 0.999, max_treedepth = 16),
      backend = "cmdstanr",
      seed = 16) -> m_suspected_vap_plus_pathogen_max_read_prop_genus_brms
}


if(save_models) {
m_suspected_vap_plus_pathogen_max_read_prop_genus_brms %>% write_rds(file = "./models/m_suspected_vap_plus_pathogen_max_read_prop_genus_brms_mixed_taxon.rds.gz", compress = "gz")
m_suspected_vap_plus_pathogen_max_read_prop_genus_brms$fit %>% write_rds(file = "./models/m_suspected_vap_plus_pathogen_max_read_prop_genus_brms_stanfit_mixed_taxon.rds.gz", compress = "gz")
}

m_suspected_vap_plus_pathogen_max_read_prop_genus_brms <- read_rds(file = "./models/m_suspected_vap_plus_pathogen_max_read_prop_genus_brms_mixed_taxon.rds.gz")

m_suspected_vap_plus_pathogen_max_read_prop_genus_brms$formula

pp_check(m_suspected_vap_plus_pathogen_max_read_prop_genus_brms)
rstan::check_hmc_diagnostics(m_suspected_vap_plus_pathogen_max_read_prop_genus_brms$fit)

m_suspected_vap_plus_pathogen_max_read_prop_genus_brms %>%
  tidybayes::get_variables()



#' #' asv
if(run_models) {
d_vap %>%
  brm(data = ., family = bernoulli,
      suspected_vap_plus_pathogen ~ (1 | subject_id) + (1 + max_read_prop_scaled | dominant_asv),
      prior = c(prior(student_t(3,0,2.5), class = Intercept)
      ),
      #sample_prior = "only", # prior predictive
      iter = 2000,
      warmup = 1000,
      chains = 4,
      cores = 4,
      control = list("adapt_delta" = 0.999, max_treedepth = 16),
      backend = "cmdstanr",
      seed = 16) -> m_suspected_vap_plus_pathogen_max_read_prop_asv_brms
}


if(save_models) {
m_suspected_vap_plus_pathogen_max_read_prop_asv_brms %>% write_rds(file = "./models/m_suspected_vap_plus_pathogen_max_read_prop_asv_brms_mixed_taxon.rds.gz", compress = "gz")
m_suspected_vap_plus_pathogen_max_read_prop_asv_brms$fit %>% write_rds(file = "./models/m_suspected_vap_plus_pathogen_max_read_prop_asv_brms_stanfit_mixed_taxon.rds.gz", compress = "gz")
}

m_suspected_vap_plus_pathogen_max_read_prop_asv_brms <- read_rds(file = "./models/m_suspected_vap_plus_pathogen_max_read_prop_asv_brms_mixed_taxon.rds.gz")

m_suspected_vap_plus_pathogen_max_read_prop_asv_brms$formula

pp_check(m_suspected_vap_plus_pathogen_max_read_prop_asv_brms)
rstan::check_hmc_diagnostics(m_suspected_vap_plus_pathogen_max_read_prop_asv_brms$fit)

m_suspected_vap_plus_pathogen_max_read_prop_asv_brms %>%
  tidybayes::get_variables()






#' functions
counterfact_asv <- function(brms_mod, exp_var_name = "shannon_low_count", draw_number = NULL, asv_vector) {
  brms_mod[["data"]] %>%
    expand(exp_var = modelr::seq_range(brms_mod[["data"]][[exp_var_name]], n = 100),
           subject_id = unique(subject_id),
           dominant_asv = asv_vector) %>%
    rename_with(.cols = exp_var, .fn = ~ exp_var_name) %>%
    tidybayes::add_fitted_draws(model = brms_mod, re_formula = NULL, n = draw_number, seed = 16) %>% # n = NULL for all draws
    identity() -> post_fitted
  return(post_fitted)
}

counterfact_genus <- function(brms_mod, exp_var_name = "shannon_low_count", draw_number = NULL, genus_vector) {
  brms_mod[["data"]] %>%
    expand(exp_var = modelr::seq_range(brms_mod[["data"]][[exp_var_name]], n = 100),
           subject_id = unique(subject_id),
           dominant_genus = genus_vector) %>%
    rename_with(.cols = exp_var, .fn = ~ exp_var_name) %>%
    tidybayes::add_fitted_draws(model = brms_mod, re_formula = NULL, n = draw_number, seed = 16) %>% # n = NULL for all draws
    identity() -> post_fitted
  return(post_fitted)
}





#' counterfactual

#' most commonly observed genera
d_vap %>%
  count(dominant_genus) %>%
  arrange(desc(n)) %>%
  slice(1:10) %>%
  pull(dominant_genus) -> top10_genera
top10_genera


#' most commonly observed ASVs
d_vap %>%
  count(dominant_asv) %>%
  arrange(desc(n)) %>%
  slice(1:10) %>%
  pull(dominant_asv) -> top10_asv
top10_asv



#' genera
top10_genera %>%
  map(.f = ~ counterfact_genus(brms_mod = m_suspected_vap_plus_pathogen_max_read_prop_genus_brms,
                               exp_var_name = "max_read_prop_scaled",
                               draw_number = 100,
                               genus_vector = .x)) %>%
  bind_rows() %>%
  identity() -> m_suspected_vap_plus_pathogen_max_read_prop_genus_cf
m_suspected_vap_plus_pathogen_max_read_prop_genus_cf


m_suspected_vap_plus_pathogen_max_read_prop_genus_cf %>%
  group_by(max_read_prop_scaled, dominant_genus) %>%
  tidybayes::median_qi(.value) %>%
  ggplot(data = .) +
  geom_line(aes(x = max_read_prop_scaled, y = .value, group = dominant_genus, color = dominant_genus))



#' asv
top10_asv %>%
  map(.f = ~ counterfact_asv(brms_mod = m_suspected_vap_plus_pathogen_max_read_prop_asv_brms,
                             exp_var_name = "max_read_prop_scaled",
                             draw_number = 100,
                             asv_vector = .x)) %>%
  bind_rows() %>%
  identity() -> m_suspected_vap_plus_pathogen_max_read_prop_asv_cf
m_suspected_vap_plus_pathogen_max_read_prop_asv_cf

if(run_plots) {
m_suspected_vap_plus_pathogen_max_read_prop_asv_cf %>%
  group_by(max_read_prop_scaled, dominant_asv) %>%
  tidybayes::median_qi(.value) %>%
  ggplot(data = .) +
  geom_line(aes(x = max_read_prop_scaled, y = .value, group = dominant_asv, color = dominant_asv))
}




#' #######################################################################
#' FORMAT ASV COUNTERFACTUAL PLOT - LINES
#' #######################################################################

dat %>%
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
  select(subject_id, subject_day, shannon, max_read_prop, log_copy_16S_per_ml_sputum, suspected_vap_plus_signs_and_pathogen, dominant_genus, dominant_asv) %>%
  mutate_at(.vars = vars(-subject_id, -subject_day, -suspected_vap_plus_signs_and_pathogen, -contains("dominant")), .funs = list("scaled" = ~ as.vector(scale(.x)))) %>%
  lm(max_read_prop ~ max_read_prop_scaled, data = .) # scale param: 0.5524, 0.2580



dat %>%
  select(dominant_asv, dominant_genus, dominant_taxon) %>%
  distinct() %>%
  #View()
  mutate(dominant_genus = replace(dominant_genus, grepl("Corynebacterium", dominant_genus), "Corynebacterium"),
         dominant_lab = glue::glue("{dominant_genus}<br>({dominant_asv})")) %>%
  identity() -> asv_genus_key
asv_genus_key



m_suspected_vap_plus_pathogen_max_read_prop_asv_cf %>%
  mutate(max_read_prop = (max_read_prop_scaled * 0.2580) + 0.5524) %>%
  left_join(asv_genus_key, by = "dominant_asv") %>%
  group_by(max_read_prop, dominant_lab) %>%
  tidybayes::median_qi(.value) %>%
  ungroup() %>%
  ggplot(data = .) +
  geom_line(aes(x = max_read_prop, y = .value, group = dominant_lab, color = dominant_lab)) +
  colorspace::scale_color_discrete_qualitative(palette = "Dark 3") +
  theme_bw() +
  theme(#legend.position = "bottom",
        axis.text.x = ggtext::element_markdown(color = "black"),
        axis.text.y = ggtext::element_markdown(color = "black"),
        axis.title.y = ggtext::element_markdown(),
        strip.text.y = ggtext::element_markdown(),
        strip.text.x = ggtext::element_markdown(),
        legend.title = ggtext::element_markdown(size = 10),
        legend.text = ggtext::element_markdown(),
        #legend.background = element_rect(color = "black", fill = "white", size = 0.25),
        strip.background = element_blank()) +
  labs(x = "Maximum ASV Proportional Abundance",
       y = "Probability of VA-LRTI",
       color = "") -> p_read_prop_mdi_by_asv
p_read_prop_mdi_by_asv


p_read_prop_mdi_by_asv +
  gghighlight::gghighlight(grepl("Pseudomonas|Staph", dominant_lab), label_params = list(fill = NA, color = "black", size = 4))


asv_genus_key %>%
  filter(grepl("f011d|9db1d",dominant_asv)) %>%
  pull(dominant_taxon)


p_read_prop_mdi_by_asv +
  annotate(geom = "text", x = 0.65, y = 0.0275, label = "Staphylococcus", size = 2, angle = 0) +
  annotate(geom = "text", x = 0.9, y = 0.036, label = "Pseudomonas", size = 2, angle = 45) -> p_read_prop_mdi_by_asv_labelled
p_read_prop_mdi_by_asv_labelled


# p_read_prop_mdi_by_asv %>%
#   ggsave(plot = ., filename = "./figs/p_read_prop_mdi_by_asv.pdf", height = 4, width = 7, units = "in")
# p_read_prop_mdi_by_asv %>%
#   ggsave(plot = ., filename = "./figs/p_read_prop_mdi_by_asv.svg", height = 4, width = 7, units = "in")
# p_read_prop_mdi_by_asv %>%
#   ggsave(plot = ., filename = "./figs/p_read_prop_mdi_by_asv.png", height = 4, width = 7, units = "in", dpi = 600)


# p_read_prop_mdi_by_asv_labelled %>%
#   ggsave(plot = ., filename = "./figs/p_read_prop_mdi_by_asv_labelled.pdf", height = 4, width = 7, units = "in")
# p_read_prop_mdi_by_asv_labelled %>%
#   ggsave(plot = ., filename = "./figs/p_read_prop_mdi_by_asv_labelled.svg", height = 4, width = 7, units = "in")
# p_read_prop_mdi_by_asv_labelled %>%
#   ggsave(plot = ., filename = "./figs/p_read_prop_mdi_by_asv_labelled.png", height = 4, width = 7, units = "in", dpi = 600)



#' #######################################################################
#' FORMAT ASV COUNTERFACTUAL PLOT - FUZZY LINES
#' #######################################################################

dat %>%
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
  select(subject_id, subject_day, shannon, max_read_prop, log_copy_16S_per_ml_sputum, suspected_vap_plus_signs_and_pathogen, dominant_genus, dominant_asv) %>%
  mutate_at(.vars = vars(-subject_id, -subject_day, -suspected_vap_plus_signs_and_pathogen, -contains("dominant")), .funs = list("scaled" = ~ as.vector(scale(.x)))) %>%
  lm(max_read_prop ~ max_read_prop_scaled, data = .) # scale param: 0.5524, 0.2580



dat %>%
  select(dominant_asv, dominant_genus, dominant_taxon) %>%
  distinct() %>%
  #View()
  mutate(dominant_genus = replace(dominant_genus, grepl("Corynebacterium", dominant_genus), "Corynebacterium"),
         dominant_lab = glue::glue("{dominant_genus}<br>({dominant_asv})")) %>%
  identity() -> asv_genus_key
asv_genus_key



m_suspected_vap_plus_pathogen_max_read_prop_asv_cf %>%
  mutate(max_read_prop = (max_read_prop_scaled * 0.2580) + 0.5524) %>%
  left_join(asv_genus_key, by = "dominant_asv") %>%
  # group_by(max_read_prop, dominant_lab) %>%
  # tidybayes::median_qi(.value) %>%
  ungroup() %>%
  ggplot(data = ., aes(x = max_read_prop, y = .value, fill = dominant_lab, color = dominant_lab)) +
  stat_lineribbon(alpha = 1/41, .width = ppoints(40)) +
  colorspace::scale_color_discrete_qualitative(palette = "Dark 3") +
  theme_bw() +
  theme(#legend.position = "bottom",
    axis.text.x = ggtext::element_markdown(color = "black"),
    axis.text.y = ggtext::element_markdown(color = "black"),
    axis.title.y = ggtext::element_markdown(),
    strip.text.y = ggtext::element_markdown(),
    strip.text.x = ggtext::element_markdown(),
    legend.title = ggtext::element_markdown(size = 10),
    legend.text = ggtext::element_markdown(),
    #legend.background = element_rect(color = "black", fill = "white", size = 0.25),
    strip.background = element_blank()) +
  coord_cartesian(ylim = c(0,0.2)) +
  labs(x = "Maximum ASV Proportional Abundance",
       y = "Probability of VA-LRTI",
       color = "",
       fill = "") -> p_read_prop_mdi_by_asv_fuzzy
p_read_prop_mdi_by_asv_fuzzy




m_suspected_vap_plus_pathogen_max_read_prop_asv_cf %>%
  mutate(max_read_prop = (max_read_prop_scaled * 0.2580) + 0.5524) %>%
  left_join(asv_genus_key, by = "dominant_asv") %>%
  # group_by(max_read_prop, dominant_lab) %>%
  # tidybayes::median_qi(.value) %>%
  ungroup() %>%
  ggplot(data = ., aes(x = max_read_prop, y = .value, fill = dominant_lab, color = dominant_lab)) +
  stat_lineribbon(alpha = 1/3, .width = .5) +
  colorspace::scale_color_discrete_qualitative(palette = "Dark 3") +
  theme_bw() +
  theme(#legend.position = "bottom",
    axis.text.x = ggtext::element_markdown(color = "black"),
    axis.text.y = ggtext::element_markdown(color = "black"),
    axis.title.y = ggtext::element_markdown(),
    strip.text.y = ggtext::element_markdown(),
    strip.text.x = ggtext::element_markdown(),
    legend.title = ggtext::element_markdown(size = 10),
    legend.text = ggtext::element_markdown(),
    #legend.background = element_rect(color = "black", fill = "white", size = 0.25),
    strip.background = element_blank()) +
  coord_cartesian(ylim = c(0,NA)) +
  labs(x = "Maximum ASV Proportional Abundance",
       y = "Probability of VA-LRTI",
       color = "",
       fill = "") -> p_read_prop_mdi_by_asv_crisp
p_read_prop_mdi_by_asv_crisp


# p_read_prop_mdi_by_asv_crisp %>%
#   ggsave(plot = ., filename = "./figs/p_read_prop_mdi_by_asv_crisp.pdf", height = 4, width = 7, units = "in")
# p_read_prop_mdi_by_asv_crisp %>%
#   ggsave(plot = ., filename = "./figs/p_read_prop_mdi_by_asv_crisp.svg", height = 4, width = 7, units = "in")
# p_read_prop_mdi_by_asv_crisp %>%
#   ggsave(plot = ., filename = "./figs/p_read_prop_mdi_by_asv_crisp.png", height = 4, width = 7, units = "in", dpi = 600)





#' #######################################################################
#' FORMAT ASV COUNTERFACTUAL PLOT - CROSS SECTIONAL SUMMARIES
#' #######################################################################

expand(exp_var = modelr::seq_range(brms_mod[["data"]][[exp_var_name]], n = 100),
       subject_id = unique(subject_id),
       dominant_asv = asv_vector) %>%
  rename_with(.cols = exp_var, .fn = ~ exp_var_name) %>%
  tidybayes::add_fitted_draws(model = brms_mod, re_formula = NULL, n = draw_number, seed = 16) %>% # n = NULL for all draws
  


m_suspected_vap_plus_pathogen_max_read_prop_asv_brms$data %>%
  expand(dominant_asv = unique(top10_asv),
         subject_id = unique(subject_id),
         max_read_prop_scaled = modelr::seq_range(x = max_read_prop_scaled, n = 10)) %>%
  tidybayes::add_fitted_draws(model = m_suspected_vap_plus_pathogen_max_read_prop_asv_brms, re_formula = NULL, n = 100, seed = 16) %>%
  mutate(max_read_prop = (max_read_prop_scaled * 0.2580) + 0.5524) %>%
  ungroup() %>%
  left_join(asv_genus_key, by = "dominant_asv") %>%
  identity() %>%
  ggplot(data = ., aes(y = dominant_lab, x = .value)) +
  tidybayes::stat_interval(.width = c(0.5,0.8,0.95)) +
  tidybayes::stat_pointinterval(.width = c(0)) +
  facet_wrap(~ factor(max_read_prop)) +
  #geom_vline(xintercept = 0, linetype = "dashed") +
  colorspace::scale_color_discrete_sequential(palette = "Blues 3", nmax = 5, order = 2:4) +
  theme_bw() +
  theme(axis.text.x = ggtext::element_markdown(color = "black"),
        axis.text.y = ggtext::element_markdown(color = "black"),
        axis.title.x = ggtext::element_markdown(color = "black"),
        legend.title = ggtext::element_markdown(color = "black"),
        legend.background = element_rect(fill = "white", color = "black", size = 0.25),
        legend.position = c(0.88,0.22)) +
  labs(x = "Regression Coefficient Relating Microbiome Disruption<br>to Duration of Mechanical Ventilation (days)",
       y = "",
       color = "Posterior<br>Credible<br>Interval")



