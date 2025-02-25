#### Script Info ####

# Authors: Dario San Segundo Molina, Ignacio Morales Castilla
# 
# GitHub repo: github.com/dario-ssm/PhenoBrassicaPests

# Aim: apply phyisiological models from phenology studies of Brassica pests 
#      to spatio-temporal contexts. Specifically to obtain long-term temporal trends
#      of phenological cues of a Brassica pest species ("Pieris rapae") in Spain with different methodologies

# Description: we use Spain02 climatic database and thermal traits (i.e. DDs, LDTs) of insect pest to predict
#              dates of emergence and voltinism across Spain and years. We analyse temporal trends for yearly variation
#              of thermal traits obtained by linear degree-days modeling.
#
# Species: Pieris rapae
# Climate dataset: Spain02 v5  (Herrera et al. 2016 and Kotlarsky et al. 2017; see http://www.meteo.unican.es/datasets/spain02 )
## Aknowledgements: The authors thank AEMET and UC by the data provided for this work (Spain02v5 gridded temperature data set). 

# 0. Load ----
library(terra)
library(tidyverse)
#library(lubridate)
library(here)
library(rTPC)
library(nls.multstart)
library(ggthemes)
library(chillR)
library(sf)
library(viridis)
source(here("R/1-functions_phenobraspests.R"))
library(mappestRisk)

pieris_data <- read_delim(here("data/gilbert_pupa.csv"),
                          delim = ";") |> 
  rename(stage = ...3,
         parasite = ...4, 
         reference = interact)

schmalensee_pieris_data <- read_delim("~/Data and scripts von Schmalensee et al. 2023/Data/full_data.txt", 
                                      delim = "\t") |> 
  filter(species == "rapae",
         life.stage == "pupa") |> 
  select(temp, dev.rate) |> 
  drop_na() 
  


# 1. Model fitting ------------------------------------------------------

#####  a) von Schmalensee 2021 ---------------------------------------------------
fit_models_pieris_rapae <- mappestRisk::fit_devmodels(temp = schmalensee_pieris_data$temp,
                                                      dev_rate = schmalensee_pieris_data$dev.rate,
                                                      model_name = "all")
my_species_name <- expression(~italic("Pieris rapae"))

plot_devmodels(temp = schmalensee_pieris_data$temp,
               dev_rate = schmalensee_pieris_data$dev.rate,
               fitted_parameters = fit_models_pieris_rapae,
               species = "Pieris rapae",
               life_stage = "Pupa")+
  labs(title = my_species_name)+
  khroma::scale_color_batlow(discrete = TRUE)+
  khroma::scale_fill_batlow(discrete = TRUE)+
  theme_few()+
  theme(legend.position = "none")
ggsave(here("figures/schmalensee_tpcs.png"),
       height = 2600,
       width = 2600,
       units = "px")
  

boots_pieris_rapae <- mappestRisk::predict_curves(temp = schmalensee_pieris_data$temp,
                                                  dev_rate = schmalensee_pieris_data$dev.rate,
                                                  fitted_parameters = fit_models_pieris_rapae,
                                                  model_name_2boot = c("schoolfield", "wang", "mod_polynomial",
                                                                       "ratkowsky", "thomas","lrf", "boatman",
                                                                       "lactin2", "joehnk", "briere2", "beta", 
                                                                       "mod_weibull",                                                                       "oneill", "kamykowski", "lactin1", "pawar"),
                                                  propagate_uncertainty = TRUE,
                                                  n_boots_samples = 100) 
                                                  
mappestRisk::plot_uncertainties(bootstrap_uncertainties_tpcs = boots_pieris_rapae,
                                temp = schmalensee_pieris_data$temp,
                                dev_rate = schmalensee_pieris_data$dev.rate,
                                species = "Pieris rapae",
                                life_stage = "Pupa")
# exclude unrealistic TPC shapes or unreasonable uncertainties 
selected_models_pieris_rapae <- fit_models_pieris_rapae #|> 
 # filter(!model_name %in% c("mod_polynomial", "oneill", "mod_weibull",
  #                          "ratkowsky", "schoolfield", "wang"))
plot_devmodels(temp = schmalensee_pieris_data$temp,
               dev_rate = schmalensee_pieris_data$dev.rate,
               fitted_parameters = selected_models_pieris_rapae,
               species = "Pieris rapae",
               life_stage = "Pupa")+
  labs(title = my_species_name)+
  khroma::scale_color_batlow(discrete = TRUE)+
  khroma::scale_fill_batlow(discrete = TRUE)

#####  b) Gilbert, 1987 ---------------------------------------------------

gilbert_pieris_data <- read_delim(here("data/gilbert_pupa.csv"),
                                  delim = ";") |> 
  rename(temp = temperature,
         dev.rate = devrate,
         life_stage = ...3,
         notes = ...4,
         reference = interact) |> 
  filter(life_stage == "pupa") |> 
  dplyr::select(temp, dev.rate)
  
  
fit_models_pieris_rapae <- mappestRisk::fit_devmodels(temp = gilbert_pieris_data$temp,
                                                      dev_rate = gilbert_pieris_data$dev.rate,
                                                      model_name = "all")
my_species_name <- expression(~italic("Pieris rapae"))

plot_devmodels(temp = gilbert_pieris_data$temp,
               dev_rate = gilbert_pieris_data$dev.rate,
               fitted_parameters = fit_models_pieris_rapae,
               species = "Pieris rapae",
               life_stage = "Pupa")+
  labs(title = my_species_name)+
  khroma::scale_color_batlow(discrete = TRUE)+
  khroma::scale_fill_batlow(discrete = TRUE)+
  theme_few()+
  theme(legend.position = "none")
ggsave(here("figures/gilbert_tpcs.png"),
       width = 2600,
       height = 2600,
       units = "px")



boots_pieris_rapae <- mappestRisk::predict_curves(temp = gilbert_pieris_data$temp,
                                                  dev_rate = gilbert_pieris_data$dev.rate,
                                                  fitted_parameters = fit_models_pieris_rapae,
                                                  model_name_2boot = c("oneill", "mod_weibull", "ratkowsky",
                                                                       "lrf", "thomas", "briere2", "beta", "boatman",
                                                                       "wang", "joehnk", "lactin1", "kamykowski", "flextpc"),
                                                  propagate_uncertainty = TRUE,
                                                  n_boots_samples = 100) 

mappestRisk::plot_uncertainties(bootstrap_uncertainties_tpcs = boots_pieris_rapae,
                                temp = gilbert_pieris_data$temp,
                                dev_rate = gilbert_pieris_data$dev.rate,
                                species = "Pieris rapae",
                                life_stage = "Pupa")
# exclude unrealistic TPC shapes or unreasonable uncertainties 
selected_models_pieris_rapae <- fit_models_pieris_rapae |> 
  filter(!model_name %in% c("mod_polynomial"))
plot_devmodels(temp = gilbert_pieris_data$temp,
               dev_rate = gilbert_pieris_data$dev.rate,
               fitted_parameters = selected_models_pieris_rapae,
               species = "Pieris rapae",
               life_stage = "Pupa")+
  labs(title = my_species_name)+
  khroma::scale_color_batlow(discrete = TRUE)+
  khroma::scale_fill_batlow(discrete = TRUE)
ggsave(here("figures/gilbert_tpcs.png"),
       width = 2600,
       height = 2600,
       units = "px")# 2. Predict phenology ------------------------------------------------------
##### 2.1. Daily res. ----------------------------------------------------------
###### a) nonlinear rate summation ----
load(here("Data/daily_tmax_df.RData"))
print(daily_tmax_df)
load(here("Data/daily_tmin_df.RData"))
print(daily_tmin_df)
daily_temperatures <- inner_join(daily_tmin_df, daily_tmax_df) |> 
  mutate(year = year(date),
         daily_tavg = map2_dbl(.x = daily_tmin,
                               .y = daily_tmax,
                               .f = ~mean(c(.x, .y))
         )
  )
rate_pred <- function(model2ratesum, temperature, fitted_params, res) {
  fitted_param_tbl_i <- fitted_params |> 
    filter(model_name == model2ratesum)
  params_i <- coef(fitted_param_tbl_i[1, ]$model_fit[[1]])
  available_models_i <- available_models |> 
    filter(model_name == model2ratesum)
  temperature_i <- temperature
  working_formula_i <- available_models_i$working_formula
  if(res == "daily") {
    rate_preds <- temperature_i |> 
      mutate(rate_pred = map_dbl(.x = daily_tavg,
                                 .f = reformulate(working_formula_i)),
             model_name = model2ratesum)
  } else if (res == "hourly") {
    rate_preds <- temperature_i |> 
      mutate(rate_pred = map_dbl(.x = temperature,
                                 .f = reformulate(working_formula_i)),
             model_name = model2ratesum)
  } else {
    stop("`res` must be either `daily` or `hourly`")
    }
  
  return(rate_preds)
}

daily_rate_preds <- tibble()

for(model_i in unique(selected_models_pieris_rapae$model_name)) {
  rate_daily_i <- rate_pred(model2ratesum = model_i,
                          temperature = daily_temperatures,
                          fitted_params = selected_models_pieris_rapae,
                          res = "daily")
  daily_rate_preds <- bind_rows(daily_rate_preds, rate_daily_i)
}

rate_sum_preds_daily <- daily_rate_preds |> 
  mutate(rate_pred = case_when(rate_pred < 0 ~ 0, # <- avoid negative summation
                               is.na(rate_pred) ~0, # <- zeros where model yields NA
                               .default = rate_pred), 
         rate_pred = 100*rate_pred) |>  # as rate_dev equals 1/days_dev
  group_by(year, model_name) |> 
  mutate(rate_sum = cumsum(rate_pred), #heat summation
         emergence_lgl = logic_ratesum(rate_sum), #threshold point
         doy = yday(date)) |> 
  summarise(day_of_emergence = sum(emergence_lgl))  #number of days until threshold is reach
   
rate_sum_preds_daily

###### b) linear DDs ----
lm_ratedev <- lm(dev.rate ~ temp,
                 data = schmalensee_pieris_data) #or gilbert pieris data

lm_pred <- function(lm_fit, tavg) {
  alpha <- coef(lm_fit)[1]
  beta <- coef(lm_fit)[2]
  preds_lm <-  alpha + beta*tavg
  return(preds_lm)
}

calc_dd <- function(lm_fit, tavg) {
  ldt <- -coef(lm_fit)[1]/coef(lm_fit)[2]
  time_id_dds <- case_when(tavg >= ldt ~ tavg-ldt,
                          .default = 0)
  return(time_id_dds)
}


dd_preds_daily <- daily_temperatures |> 
  mutate(model_name = "linear") |> 
  group_by(year, model_name) |> 
  mutate(rate_pred = map_dbl(.x = daily_tavg,
                             .f = ~lm_pred(lm_fit = lm_ratedev,
                                           .x))) |> 
  mutate(rate_pred = if_else(rate_pred <0,
                             0,
                             rate_pred)) |>
  mutate(rate_sum = map_dbl(.x = daily_tavg,
                            .f = ~calc_dd(lm_fit = lm_ratedev,
                                          tavg = .x)), #heat summation
         rate_sum = cumsum(rate_sum),
         emergence_lgl = map_dbl(.x = rate_sum,
                                 .f = ~logic_dd(heat_units = 1/coef(lm_ratedev)[2],
                                                rate_cumsum = .x)), #threshold point
         doy = yday(date)) |> 
  summarise(day_of_emergence = sum(emergence_lgl))

###### c) combined  ----
doe_preds_daily <- rate_sum_preds_daily |> 
  full_join(dd_preds_daily) |> 
  mutate(time_slot = "daily")

##### 2.2. Hourly res. ----------------------------------------------------------
##  we use chillR package to simulate temperature variation across hours within days.
###### a) nonlinear rate summation ----

daily_temps_doy <- daily_temperatures  |> 
  mutate(JDay = yday(date)) |>  
  rename(Tmin = daily_tmin,
         Tmax = daily_tmax) 

hourly_temp_doy <- chillR::make_hourly_temps(latitude = 40.5,
                                             year_file = daily_temps_doy,
                                             keep_sunrise_sunset = TRUE )  |> 
  pivot_longer(cols = 10:33,
               names_to = "hour_of_day",
               values_to = "temperature") |>  
  rename(tmax = Tmax,
         tmin = Tmin,
         doy = JDay,
         sunrise = Sunrise,
         sunset = Sunset,
         daylength = Daylength) 


hourly_rate_preds <- tibble()

for(model_i in unique(selected_models_pieris_rapae$model_name)) {
  rate_hourly_i <- rate_pred(model2ratesum = model_i,
                            temperature = hourly_temp_doy,
                            fitted_params = selected_models_pieris_rapae,
                            res = "hourly")
  hourly_rate_preds <- bind_rows(hourly_rate_preds, rate_hourly_i)
}
hourly_rate_preds

rate_sum_preds_hourly <- hourly_rate_preds |> 
  mutate(rate_pred = case_when(rate_pred < 0 ~ 0, # <- avoid negative summation
                               is.na(rate_pred) ~0, # <- zeros where model yields NA
                               .default = rate_pred), 
         rate_pred = rate_pred*100/24) |>  # 24h a day
  group_by(year, model_name) |> 
  mutate(rate_sum = cumsum(rate_pred), #heat summation
         emergence_lgl = logic_ratesum(rate_sum), #threshold point
         doy = yday(date)) |> 
  summarise(day_of_emergence = round(sum(emergence_lgl)/24))  #number of days until threshold is reach

rate_sum_preds_hourly

###### b) linear DDs ----
dd_preds_hourly <- hourly_temp_doy |> 
  mutate(model_name = "linear") |> 
  group_by(year, model_name) |> 
  mutate(rate_pred = map_dbl(.x = temperature,
                             .f = ~lm_pred(lm_fit = lm_ratedev,
                                           .x))) |> 
  mutate(rate_pred = if_else(rate_pred <0,
                             0,
                             rate_pred)) |> #fraction of 1/24 hours
  mutate(rate_sum = map_dbl(.x = temperature,
                            .f = ~calc_dd(lm_fit = lm_ratedev,
                                          tavg = .x)), #heat summation
         rate_sum = cumsum(rate_sum),
         emergence_lgl = map_dbl(.x = rate_sum,
                                 .f = ~logic_dd(heat_units = 24*1/coef(lm_ratedev)[2],
                                                rate_cumsum = .x)), #threshold point
         doy = yday(date)) |> 
  summarise(day_of_emergence = round(sum(emergence_lgl)/24))

###### c) combined  ----
doe_preds_hourly <- rate_sum_preds_hourly |> 
  full_join(dd_preds_hourly) |> 
  mutate(time_slot = "hourly")

# 3. Comparisons (long-term trends) ---------------------------------------
##### a) Visual comparison ----

# incorporate observational data from Gordo & Sanz (2006)
pieris_obs_trends <- read_csv(here("data/gordo_sanz_2006_doe.csv")) |> 
  rename(day_of_emergence = doe) |> 
  relocate(year, day_of_emergence) |> 
  group_by(year) |> 
  slice(1)

rate_sum_validation <- doe_preds_daily |> 
  bind_rows(doe_preds_hourly)
doe_nonlinear_trends <- ggplot(rate_sum_validation,
                               aes(x = year, 
                                   y = day_of_emergence))+
  geom_point(aes(color = time_slot),
             alpha = .6)+
  geom_line(aes(color = time_slot),
            linetype = "dashed",
            alpha = .6)+
  geom_smooth(aes(color = time_slot,
                  fill = time_slot))+
  facet_wrap(~model_name)+
  theme_clean()+
  scale_color_manual(values = c("#3A848B", "#9C4E64"))+
  scale_fill_manual(values = c("#3A848B", "#9C4E64"))+
  geom_point(data = pieris_obs_trends,
             aes(x = year,
                 y = day_of_emergence),
             color = "#E9C86B",
             alpha = .6)+
  geom_line(data = pieris_obs_trends,
            aes(x = year,
                y = day_of_emergence),
            color = "#E9C86B",
            linetype = "dashed",
            alpha = .6)+
  geom_smooth(data = pieris_obs_trends,
              aes(x = year, y  = day_of_emergence),
              color = "#E9C86B",
              fill = "#E9C86B")+
  labs(x = "Year",
       y = "Day of first adult emergence")
doe_nonlinear_trends

#gilbert
ggsave(filename = here("figures/doe_validation_gilbert.png"),
       width = 2600, height = 2600, units = "px")

ggsave(filename = here("figures/doe_validation_gilbert.svg"))

#schmalensee

ggsave(filename = here("figures/doe_validation_schmalensee2023.png"),
       width = 2600, height = 2600, units = "px")

ggsave(filename = here("figures/doe_validation_schmalensee.svg"))


##### b) RMSE ----
preds_vs_obs_trends <- rate_sum_validation  |>
  filter(year > 1951 & year < 2005) |> 
  group_by(model_name, time_slot) |> 
  mutate(observed_emergence = pieris_obs_trends$day_of_emergence) |> 
  summarise(rmsep = chillR::RMSEP(predicted = day_of_emergence,
                                  observed = observed_emergence)) |>  
  arrange(rmsep) %>% 
  print() #gaussian, oneill and weibull are better models, followed by ssilow, hourly resolution does not improve accuracy
#for worse models, such as linear, briere, ssi and lactin, hourly resolution clearly improve accuracy.

transfer_rank_resolution <- preds_vs_obs_trends %>% 
  group_by(time_slot) %>% 
  summarise(rmsep = mean(rmsep)) %>% 
  arrange(rmsep) %>% 
  print()# overall hourly rmsep is lower

transfer_rank_model <- preds_vs_obs_trends %>% 
  group_by(model_name) %>% 
  summarise(rmsep = mean(rmsep)) %>% 
  arrange(rmsep) %>% 
  print() # overall mod_weibull, beta, flextpc, schoolfield and oneill

#re-plot the graph (manually) annotating rmsep

preds_obs_compare_rmsep <- rate_sum_validation |> 
  full_join(preds_vs_obs_trends) |> 
  group_by(year, model_name, time_slot) |> 
  summarise(rmsep = min(rmsep)) |>  
  arrange(rmsep) |> 
  print()
rmsep_text <-  preds_obs_compare_rmsep |>  
  group_by(model_name, time_slot) |>  
  summarise(rmsep = mean(rmsep)) |> 
  arrange(rmsep)
rmsep_order <- rmsep_text |> 
  pull(model_name)
rmsep_order_min <- rmsep_text |> 
  slice_min(rmsep) |>
  arrange(rmsep) |> 
  pull(model_name)
rmsep_values <- rmsep_text |> 
  mutate(rmsep =  case_when(time_slot == "daily" ~paste("d-RMSEP =",round(rmsep,2)),
                             time_slot == "hourly" ~paste("h-RMSEP =",round(rmsep,2))),
         year = 2004,
         day_of_emergence = case_when(time_slot == "daily" ~165,
                                      time_slot == "hourly" ~145))
preds_obs_compare_plot <- ggplot()+
  geom_point(data = rate_sum_validation |> 
               filter(year %in% pieris_obs_trends$year),
             aes(x = year, y = day_of_emergence, color = time_slot),
             alpha = .66)+
  geom_line(data = rate_sum_validation |> 
              filter(year %in% pieris_obs_trends$year),
            aes(x = year, y = day_of_emergence, color = time_slot),
            linetype = "longdash",
            alpha = .66)+
  geom_smooth(data =  rate_sum_validation |> 
                filter(year %in% pieris_obs_trends$year),
              aes(x = year, y = day_of_emergence, 
                  color = time_slot, 
                  fill = time_slot))+
  scale_color_manual(values = c("#3A848B", "#E85038"))+
  scale_fill_manual(values = c("#3A848B", "#E85038"))+
  geom_point(data = pieris_obs_trends, 
             aes(x = year, y = day_of_emergence),
             color = "#E9C86B",
             fill = "#E9C86B",
             alpha = .66)+
  geom_line(data = pieris_obs_trends, 
            aes(x = year, y = day_of_emergence),
            color = "#E9C86B",
            linetype = "longdash",
            alpha = .66)+
  geom_smooth(data = pieris_obs_trends, 
              aes(x = year, y = day_of_emergence),
              color = "#E9C86B",
              fill = "#E9C86B")+
  labs(x = "Year",
       y = "Day of first adult emergence",
       color = "Resolution (model)",
       fill = "Resolution (model)")+
  facet_wrap(~factor(model_name, levels = rmsep_order_min))+
  theme_few()+
  theme(strip.text.x = element_text(size = 10, 
                                    color = "grey42", 
                                    face = "bold"))+
  geom_text(data = rmsep_values, 
            aes(label = rmsep,
                x = year,
                y = day_of_emergence), 
            size = 3,
            fontface = "italic",hjust = "right")
preds_obs_compare_plot  

## save for schmalensee
ggsave(filename = here("figures/schmalensee_modelpreds_validate_rmsep.png"),
       width = 25,
       height = 25,
       units = "cm")
ggsave(filename = here("figures/schmalensee_modelpreds_validate_rmsep.svg"),
       width = 25,
       height = 25,
       units = "cm")

## save for gilbert
ggsave(filename = here("figures/gilbert_modelpreds_validate_rmsep.png"),
       width = 25,
       height = 25,
       units = "cm")
ggsave(filename = here("figures/gilbert_modelpreds_validate_rmsep.svg"),
       width = 25,
       height = 25,
       units = "cm")


# 4. Thermal regimes ---------------------------------------

daily_temperatures_doy <- daily_temperatures |> 
  filter(year %in% 1952:2004) |> 
  mutate(doy = yday(date)) |> 
  rename(tmin = daily_tmin,
         tmax = daily_tmax,
         tavg = daily_tavg) |> 
  pivot_longer(cols = c(2, 3, 5),
               names_to = "temp_var",
               values_to = "daily_temperature")
  
ggplot(daily_temperatures_doy, aes(x = doy, y = daily_temperature)) +
  geom_point(aes(color = temp_var), alpha = .1)+
  scale_color_manual(values = c("#F1CB56", "#E36944", "#7FD3C3"))+
  ggdark::dark_theme_minimal()+
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "none")+
  labs(x = "Day of Year",
       y = "Daily temperature",
       color = "Temperature variables")

ggsave(here("figures/daily_temperature_regimes.png"),
       width = 2400,
       height = 2400,
       units = "px")
ggsave(here("figures/daily_temperature_regimes.svg"),
       width = 2400,
       height = 2400,
       units = "px") 

daily_temperatures_doy_decades <- daily_temperatures_doy |> 
  mutate(decade = case_when(year %in% c(1950:1959) ~ "1950's",
                            year %in% c(1960:1969) ~ "1960's",
                            year %in% c(1970:1979) ~ "1970's",
                            year %in% c(1980:1989) ~ "1980's",
                            year %in% c(1990:1999) ~ "1990's",
                            year %in% c(2000:2009) ~ "2000's",
                            year %in% c(2010:2015) ~ "2010's"
                            )
         ) |> 
  group_by(decade) |> 
  mutate(decade_count_day = rep(1:(n()/3), each = 3)) |>  # Assign the same value to all rows for the same doy
  ungroup()

ggplot(daily_temperatures_doy_decades, aes(x = decade_count_day, y = daily_temperature)) +
  #geom_point(aes(color = temp_var), alpha =.3)+
  geom_line(aes(color = temp_var),
            linewidth = .5,
            alpha = 1)+
  scale_color_manual(values = c("#F1CB56", "#E36944", "#7FD3C3"))+
  scale_x_reverse()+
  ggdark::dark_theme_minimal()+
  labs(x = "Days",
       y = "Daily temperature",
       color = "Temperature variables")+
  coord_flip()+
  facet_wrap(~decade,
             scales = "free_x",
             nrow = 1
             )+
  theme(legend.position = "none")

ggsave(here("figures/daily_temperature_regimes_decades.png"),
       width = 2500,
       height = 4000,
       units = "px")

ggsave(here("figures/daily_temperature_regimes_decades.svg"),
       width = 2500,
       height = 4000,
       units = "px")

## and plot oneill

plot_devmodels(temp = gilbert_pieris_data$temp,
               dev_rate = gilbert_pieris_data$dev.rate,
               fitted_parameters = selected_models_pieris_rapae |> 
                 filter(model_name == "oneill"),
               species = "Pieris rapae",
               life_stage = "Pupa")+
  labs(title = NULL,
       subtitle = NULL,
       y = NULL,
       x = NULL)+
  scale_fill_manual(values = "#f39e71ff")+
  scale_color_manual(values = "#f39e71ff")+
  theme_minimal()+
  theme(legend.position = "none",
        axis.text = element_blank())+
  facet_null()
  
ggsave(here("figures/oneill.svg"),
       width = 600,
       height = 600,
       units = "px")

# 5. Data from CBMS -------------------------------------------------------


