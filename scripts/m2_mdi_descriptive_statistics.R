#' ###########################################
#' MDI descriptive statistics & plots
#' - depends on libraries and data
#' ###########################################

#' how many subjects/specimens to inform MDI estimates

dat %>%
  select(subject_id, date, subject_day, shannon, max_read_prop, log_copy_16S_per_ml_sputum) %>%
  gather(key = "mdi", value = "mdi_value", -subject_id, -date, -subject_day) %>%
  filter(!is.na(mdi_value)) %>%
  mutate(time_period = case_when(subject_day >= 0 & subject_day < 7 ~ "Week 1",
                                 subject_day >= 7 & subject_day < 14 ~ "Week 2",
                                 subject_day >= 14 ~ "Beyond Week 2")) %>%
  count(mdi) %>%
  mutate(mdi = case_when(mdi == "shannon" ~ "Shannon Diversity (base e)",
                         mdi == "max_read_prop" ~ "Max ASV Proportional Abundance",
                         mdi == "log_copy_16S_per_ml_sputum" ~ "Log Copies 16S rRNA per mL Sputum")) %>%
  identity() -> t_specimen_totals
t_specimen_totals


dat %>%
  select(subject_id, date, subject_day, shannon, max_read_prop, log_copy_16S_per_ml_sputum) %>%
  gather(key = "mdi", value = "mdi_value", -subject_id, -date, -subject_day) %>%
  filter(!is.na(mdi_value)) %>%
  mutate(time_period = case_when(subject_day >= 0 & subject_day < 7 ~ "Week 1",
                                 subject_day >= 7 & subject_day < 14 ~ "Week 2",
                                 subject_day >= 14 ~ "Beyond Week 2")) %>%
  mutate(time_period = factor(time_period, levels = c("Week 1","Week 2", "Beyond Week 2"))) %>%
  group_by(mdi,time_period) %>%
  summarise(subject_n = n_distinct(subject_id),
            specimen_n = n()) %>%
  mutate(mdi = case_when(mdi == "shannon" ~ "Shannon Diversity (base e)",
                         mdi == "max_read_prop" ~ "Max ASV Proportional Abundance",
                         mdi == "log_copy_16S_per_ml_sputum" ~ "Log Copies 16S rRNA per mL Sputum")) %>%
  identity() -> t_weekly_subject_specimen_totals
t_weekly_subject_specimen_totals



t_weekly_subject_specimen_totals %>%
  rename(MDI = "mdi",
         `LTACH Course` = time_period,
         `Active Subject Number` = subject_n,
         `Specimens Collected` = specimen_n) %>%
  gt::gt()


t_weekly_subject_specimen_totals %>%
  rename(MDI = "mdi",
         `LTACH Course` = time_period,
         `Active Subject Number` = subject_n,
         `Specimens Collected` = specimen_n) %>%
  knitr::kable(format = "html") %>% 
  kableExtra::kable_styling(full_width = FALSE) %>% 
  kableExtra::column_spec(column = c(1,2,3,4), color = "black") %>%
  kableExtra::collapse_rows(columns = c(1))



#' plot MDIs over categorical time post enrollment
dat %>%
  select(subject_id, subject_day, shannon, max_read_prop, copy_number_per_ml_sputum, log_copy_16S_per_ml_sputum) %>%
  distinct() %>%
  gather(key = "mdi", value = "mdi_value", -subject_id, -subject_day) %>%
  mutate(mdi_name = stringr::str_to_title(gsub("_"," ", mdi)),
         mdi_name = stringr::str_replace(mdi_name, "Ml", "mL"),
         mdi_name = stringr::str_replace(mdi_name, "Log Copy 16s", "Log<sub>10</sub> Copies 16S rRNA"),
         mdi_name = stringr::str_replace(mdi_name, "Shannon", "Shannon Diversity (base e)"),
         mdi_name = stringr::str_replace(mdi_name, "Max Read Prop", "Max ASV Proportional Abundance")) %>%
  mutate(subject_day_category = ifelse(subject_day == 0, "Enrollment", ifelse(subject_day > 0 & subject_day < 7, "Week 1", ifelse(subject_day >= 7 & subject_day < 14, "Week 2", ifelse(subject_day >= 14, "Beyond<br>Week 2", NA)))),
         subject_day_category = factor(subject_day_category, levels = c("Enrollment", "Week 1", "Week 2", "Beyond<br>Week 2"))) %>%
  filter(mdi != "copy_number_per_ml_sputum") %>%
  ggplot(data = ., aes(x = subject_day_category, y = mdi_value)) +
  geom_boxplot(fill = "white", color = "black") +
  #scale_fill_viridis_d(option = "plasma") +
  facet_wrap(facets = ~ mdi_name, scales = "free") +
  theme_bw() +
  theme(axis.text = ggtext::element_markdown(color = "black"),
        axis.text.x = ggtext::element_markdown(angle = 325, hjust = 0.5, vjust = 0.5),
        strip.background = element_blank(),
        strip.text = ggtext::element_markdown(),
        legend.position = "none",
        plot.title.position = "plot") +
  labs(x = "", y = "MDI Value", title = "") -> p_mdi_longitudinal

p_mdi_longitudinal

# p_mdi_longitudinal %>%
#   ggsave(plot = ., filename = "./figs/p_mdi_longitudinal.pdf", height = 4, width = 8, units = "in")
# p_mdi_longitudinal %>%
#   ggsave(plot = ., filename = "./figs/p_mdi_longitudinal.svg", height = 4, width = 8, units = "in")
# p_mdi_longitudinal %>%
#   ggsave(plot = ., filename = "./figs/p_mdi_longitudinal.png", height = 4, width = 8, units = "in", dpi = 600)












