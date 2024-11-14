#' #############################################
#' load libraries and set seed
#' #############################################
library(tidyverse)
library(tidybayes)
library(brms)
set.seed(16)
run_models = FALSE
save_models = FALSE
run_plots = FALSE




#' #############################################
#' mixed effects model for MDI ~ ventilation days
#' #############################################

if(run_models) {
dat %>%
  select(subject_id, subject_day, shannon, max_read_prop, log_copy_16S_per_ml_sputum) %>%
  gather(key = "mdi", value = "mdi_value", -subject_id, -subject_day) %>%
  filter(!is.na(mdi_value)) %>%
  group_by(mdi) %>%
  nest() %>%
  mutate(model = map(.x = data, .f = ~ brms::brm(data = .x,
                                                 family = "gaussian",
                                                 formula = mdi_value ~ subject_day + 0 + (1 | subject_id),
                                                 seed = 16,
                                                 cores = 4,
                                                 backend = "cmdstanr"))) %>%
  identity() -> d_mdi_vs_ventdays_model
d_mdi_vs_ventdays_model
}


d_mdi_vs_ventdays_model %>%
  mutate(postsum = map(.x = model, .f = ~ brms::posterior_summary(.x) %>%
                         as_tibble(rownames = "param"))) %>%
  select(mdi, postsum) %>%
  unnest(col = "postsum") %>%
  ungroup() %>%
  filter(param == "b_subject_day") %>%
  mutate(mdi = case_when(mdi == "shannon" ~ "Shannon Diversity (base e)",
                         mdi == "max_read_prop" ~ "Max ASV Proportional Abundance",
                         mdi == "log_copy_16S_per_ml_sputum" ~ "Log Copies 16S rRNA per mL Sputum")) %>%
  identity() %>%
  gt::gt()
  




d_mdi_vs_ventdays_model %>%
  mutate(postsum = map(.x = model, .f = ~ brms::posterior_summary(.x) %>%
                         as_tibble(rownames = "param"))) %>%
  select(mdi, postsum) %>%
  unnest(col = "postsum") %>%
  ungroup() %>%
  filter(param == "b_subject_day") %>%
  mutate(mdi = case_when(mdi == "shannon" ~ "Shannon Diversity<br>(base e)",
                         mdi == "max_read_prop" ~ "Max ASV<br>Proportional<br>Abundance",
                         mdi == "log_copy_16S_per_ml_sputum" ~ "Log Copies<br>16S rRNA<br>per mL Sputum")) %>%
  ggplot(data = .) +
  geom_segment(aes(x = Q2.5, xend = Q97.5, y = mdi, yend = mdi)) +
  geom_point(aes(x = Estimate, y = mdi)) +
  geom_vline(xintercept = 0, linetype = 2) +
  theme_bw() +
  theme(axis.text.x = ggtext::element_markdown(color = "black"),
        axis.text.y = ggtext::element_markdown(color = "black"),
        axis.title.x = ggtext::element_markdown()) +
  labs(x = "Regression Coefficient Relating Microbiome Disruption<br>to Duration of Mechanical Ventilation (days)",
       y = "MDI")



d_mdi_vs_ventdays_model %>%
  mutate(postsamp = map(.x = model, .f = ~ posterior::as_draws_df(.x$fit) %>%
                          as_tibble() %>%
                          select(b_subject_day)
                        )
         ) %>%
  select(mdi, postsamp) %>%
  unnest(col = "postsamp") %>%
  ungroup() %>%
  mutate(mdi = case_when(mdi == "shannon" ~ "Shannon Diversity",
                         mdi == "max_read_prop" ~ "Maximum ASV<br>Proportional<br>Abundance",
                         mdi == "log_copy_16S_per_ml_sputum" ~ "Total<br>Bacterial<br>Abundance<br>(log copies<br>16S rRNA per<br>mL sputum)")) %>%
  identity() %>%
  ggplot(data = ., aes(x = b_subject_day, y = mdi)) +
  tidybayes::stat_interval(.width = c(0.5,0.8,0.95)) +
  tidybayes::stat_pointinterval(.width = c(0)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  colorspace::scale_color_discrete_sequential(palette = "Blues 3", nmax = 5, order = 2:4) +
  theme_bw() +
  theme(axis.text.x = ggtext::element_markdown(color = "black"),
        axis.text.y = ggtext::element_markdown(color = "black"),
        axis.title.x = ggtext::element_markdown(color = "black"),
        legend.title = ggtext::element_markdown(color = "black"),
        legend.background = element_rect(fill = "white", color = "black", size = 0.25),
        legend.position = c(0.88,0.22)) +
  labs(x = "Linear Regression Coefficient Relating MDI<br>to Duration of Mechanical Ventilation (days)",
       y = "",
       color = "Posterior<br>Credible<br>Interval"
       ) -> p_mdi_vs_ventdays
p_mdi_vs_ventdays




# p_mdi_vs_ventdays %>%
#   ggsave(filename = "./figs/p_mdi_vs_ventdays.pdf", height = 4, width = 5, units = "in")
# p_mdi_vs_ventdays %>%
#   ggsave(filename = "./figs/p_mdi_vs_ventdays.svg", height = 4, width = 5, units = "in")
# p_mdi_vs_ventdays %>%
#   ggsave(filename = "./figs/p_mdi_vs_ventdays.png", height = 4, width = 5, units = "in", dpi = 600)
