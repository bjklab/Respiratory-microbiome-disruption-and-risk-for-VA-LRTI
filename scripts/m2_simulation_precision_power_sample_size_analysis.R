#' #############################################
#' load libraries and set seed
#' #############################################
library(tidyverse)
library(tidybayes)
library(brms)
set.seed(16)






#' #############################################################
#' #############################################################
#' 
#' DETECTABLE DIFFERENCE IN MDI
#' 
#' #############################################################
#' #############################################################

#' #############################################################
#' simulation for power/precision analysis
#' #############################################################

#' simulate data
tibble(va_lrti = rbinom(n = 10000, size = 1, prob = 0.1)) %>%
  mutate(mdi_0.1 = map_dbl(.x = va_lrti, .f = ~ ifelse(.x == 1, rnorm(1, mean = 0.1, sd = 1), rnorm(1, mean = 0, sd = 1))),
         mdi_0.2 = map_dbl(.x = va_lrti, .f = ~ ifelse(.x == 1, rnorm(1, mean = 0.2, sd = 1), rnorm(1, mean = 0, sd = 1))),
         mdi_0.3 = map_dbl(.x = va_lrti, .f = ~ ifelse(.x == 1, rnorm(1, mean = 0.3, sd = 1), rnorm(1, mean = 0, sd = 1))),
         mdi_0.4 = map_dbl(.x = va_lrti, .f = ~ ifelse(.x == 1, rnorm(1, mean = 0.4, sd = 1), rnorm(1, mean = 0, sd = 1))),
         mdi_0.5 = map_dbl(.x = va_lrti, .f = ~ ifelse(.x == 1, rnorm(1, mean = 0.5, sd = 1), rnorm(1, mean = 0, sd = 1))),
         mdi_0.6 = map_dbl(.x = va_lrti, .f = ~ ifelse(.x == 1, rnorm(1, mean = 0.6, sd = 1), rnorm(1, mean = 0, sd = 1))),
         mdi_0.7 = map_dbl(.x = va_lrti, .f = ~ ifelse(.x == 1, rnorm(1, mean = 0.7, sd = 1), rnorm(1, mean = 0, sd = 1))),
         mdi_0.8 = map_dbl(.x = va_lrti, .f = ~ ifelse(.x == 1, rnorm(1, mean = 0.8, sd = 1), rnorm(1, mean = 0, sd = 1))),
         mdi_0.9 = map_dbl(.x = va_lrti, .f = ~ ifelse(.x == 1, rnorm(1, mean = 0.9, sd = 1), rnorm(1, mean = 0, sd = 1))),
         mdi_1.0 = map_dbl(.x = va_lrti, .f = ~ ifelse(.x == 1, rnorm(1, mean = 1.0, sd = 1), rnorm(1, mean = 0, sd = 1))),
         ) %>%
  gather(key = "mdi", value = "mdi_measure", -va_lrti) %>%
  mutate(actual_microbiome_disruption = as.numeric(gsub("mdi_","",mdi))) %>%
  identity() -> sim_data
sim_data



# run models with simulated data
sim_data %>%
  group_by(mdi) %>%
  nest() %>%
  mutate(model = map(.x = data,
                     .f = ~ brm(data = .x,
                                family = "bernoulli",
                                formula = va_lrti ~ mdi_measure,
                                backend = "cmdstanr",
                                seed = 16,
                                cores = 4))) %>%
  mutate(model_100 = map(.x = model, .f = ~ update()))
  identity() -> sim_model
sim_model



#' update model with different sample sizes
sim_model %>%
  mutate(data_100 = map(.x = data, .f = ~ .x %>% slice_sample(n = 100)),
         data_200 = map(.x = data, .f = ~ .x %>% slice_sample(n = 200)),
         data_300 = map(.x = data, .f = ~ .x %>% slice_sample(n = 300)),
         data_400 = map(.x = data, .f = ~ .x %>% slice_sample(n = 400)),
         data_500 = map(.x = data, .f = ~ .x %>% slice_sample(n = 500)),
         data_600 = map(.x = data, .f = ~ .x %>% slice_sample(n = 600)),) %>%
  mutate(mod_100 = map2(.x = data_100, .y = model, .f = ~ update(.y, newdata = .x)),
         mod_200 = map2(.x = data_200, .y = model, .f = ~ update(.y, newdata = .x)),
         mod_300 = map2(.x = data_300, .y = model, .f = ~ update(.y, newdata = .x)),
         mod_400 = map2(.x = data_400, .y = model, .f = ~ update(.y, newdata = .x)),
         mod_500 = map2(.x = data_500, .y = model, .f = ~ update(.y, newdata = .x)),
         mod_600 = map2(.x = data_600, .y = model, .f = ~ update(.y, newdata = .x)),) %>%
  identity() -> sim_power
sim_power



#' posterior samples of parameters
sim_power %>%
  mutate(ps_100 = map(.x = mod_100,
                      .f = ~ .x %>%
                        pluck("fit") %>%
                        posterior::as_draws_df() %>%
                        as_tibble() %>%
                        select(b_mdi_measure)
  ),
  ps_200 = map(.x = mod_200,
               .f = ~ .x %>%
                 pluck("fit") %>%
                 posterior::as_draws_df() %>%
                 as_tibble() %>%
                 select(b_mdi_measure)
  ),
  ps_300 = map(.x = mod_300,
               .f = ~ .x %>%
                 pluck("fit") %>%
                 posterior::as_draws_df() %>%
                 as_tibble() %>%
                 select(b_mdi_measure)
  ),
  ps_400 = map(.x = mod_400,
               .f = ~ .x %>%
                 pluck("fit") %>%
                 posterior::as_draws_df() %>%
                 as_tibble() %>%
                 select(b_mdi_measure)
  ),
  ps_500 = map(.x = mod_500,
               .f = ~ .x %>%
                 pluck("fit") %>%
                 posterior::as_draws_df() %>%
                 as_tibble() %>%
                 select(b_mdi_measure)
  ),
  ps_600 = map(.x = mod_600,
               .f = ~ .x %>%
                 pluck("fit") %>%
                 posterior::as_draws_df() %>%
                 as_tibble() %>%
                 select(b_mdi_measure)
  ),
  ) %>%
  identity() -> sim_post
sim_post



sim_post %>%
  select(mdi, contains("ps")) %>%
  gather(key = "sample_size", value = "posterior_samples", -mdi) %>%
  mutate(sample_size = as.numeric(gsub("ps_","",sample_size)),
         effect_size = as.numeric(gsub("mdi_","",mdi))) %>%
  unnest(cols = contains("posterior")) %>%
  ungroup() %>%
  mutate(b_mdi_OR = exp(b_mdi_measure)) %>%
  identity()









#' summarise parameter estimates
sim_power %>%
  mutate(posterior_100 = map(.x = mod_100,
                             .f = ~ brms::posterior_summary(.x) %>%
                               as_tibble(rownames = "param") %>%
                               filter(param == "b_mdi_measure")
                             ),
         posterior_200 = map(.x = mod_200,
                             .f = ~ brms::posterior_summary(.x) %>%
                               as_tibble(rownames = "param") %>%
                               filter(param == "b_mdi_measure")
         ),
         posterior_300 = map(.x = mod_300,
                             .f = ~ brms::posterior_summary(.x) %>%
                               as_tibble(rownames = "param") %>%
                               filter(param == "b_mdi_measure")
         ),
         posterior_400 = map(.x = mod_400,
                             .f = ~ brms::posterior_summary(.x) %>%
                               as_tibble(rownames = "param") %>%
                               filter(param == "b_mdi_measure")
         ),
         posterior_500 = map(.x = mod_500,
                             .f = ~ brms::posterior_summary(.x) %>%
                               as_tibble(rownames = "param") %>%
                               filter(param == "b_mdi_measure")
         ),
         posterior_600 = map(.x = mod_600,
                             .f = ~ brms::posterior_summary(.x) %>%
                               as_tibble(rownames = "param") %>%
                               filter(param == "b_mdi_measure")
         ),
         ) %>%
  identity() -> sim_summ
sim_summ





sim_summ %>%
  select(mdi, contains("posterior")) %>%
  gather(key = "sample_size", value = "posterior", -mdi) %>%
  mutate(sample_size = as.numeric(gsub("posterior_","",sample_size)),
         effect_size = as.numeric(gsub("mdi_","",mdi))) %>%
  unnest(cols = contains("posterior")) %>%
  mutate_at(.vars = vars(Estimate, Q2.5, Q97.5), .funs = list("OR" = ~ exp(.x))) %>%
  identity() %>%
  ggplot(data = .) +
  geom_segment(aes(x = Q2.5_OR, xend = Q97.5_OR, y = effect_size, yend = effect_size)) +
  geom_point(aes(x = Estimate_OR, y = effect_size)) +
  facet_wrap(facets = ~ factor(paste0(sample_size," specimens"))) +
  geom_vline(xintercept = 1, linetype = 2)





