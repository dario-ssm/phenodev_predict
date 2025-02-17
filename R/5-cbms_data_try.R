# 1. CBMS data ----
library(tidyverse)
library(here)
library(lubridate)
library(leaflet)
library(leaflet.extras2)
library(RColorBrewer)
library(terra)
library(sf)
library(easyclimate)

## load pieris rapae only data counting
piepra_raw_count_data_csv <- read_delim("~/GitHub/dev_rapae/dev_rapae/data_source/peticio_sansegundo.csv") |> 
  mutate(Nindiv = if_else(Nindiv == 0.5, # <- misread with R
                          1,
                          Nindiv))

## a) Regular data set -----------------------------------------------------

piepra_raw_count_data <- piepra_raw_count_data_csv |> 
  rename(id_transect = IDitin,
         id_species = IDesp,
         year = Any,
         date = Date,
         n_indiv = Nindiv,
         length = longitud,
         lat = LAT,
         lon = LNG) |> 
  select(-year) |> 
  mutate(date = lubridate::dmy(date),
         year = lubridate::year(date),
         id_transect = as_factor(id_transect)) |> 
  relocate(id_transect, id_species, year, date, n_indiv, length, lat, lon) |> 
  glimpse()

## a) Combined data set -----------------------------------------------------
cbms_all_transects <- read_delim("~/GitHub/dev_rapae/dev_rapae/data_source/m_visit_sub.csv") |> 
 rename(id_transect = SITE_ID,
        date = DATE)|> 
  select(id_transect, date) |> 
  mutate(id_transect = as_factor(id_transect),
         year = year(date))

pieris_pres_abs <- full_join(piepra_raw_count_data,
                             cbms_all_transects) |>
  group_by(id_transect) |> 
  mutate(pres_abs = map_chr(.x = n_indiv,
                            ~if_else(is.na(.x),
                                     "absent",
                                     "present"))) |> 
  mutate(n_indiv = map2_dbl(.x = pres_abs,
                            .y = n_indiv,
                            .f = ~if_else(.x == "absent",
                                          0,
                                          .y))) |> 
  mutate(doy = yday(date))

## c) Data wrangling -----------------------------------------------------
## examine trends

## now we will apply exclusion criteria:
##  
## (1) we'll exclude complete transects with >50% data absences than presences.
## since populations at these places may be transient or population dynamics may
## reflect other processes rather than heat accumulation.

## (2) we'll exclude anomalous years for each transect, i.e., those
##  with late first appearances (i.e., > doy150) with no previous or very early absences
##  followed by no data for months that may indicate lack of sampling rather than lack of presence.
##  Similarly, we'll exclude years within tranects having exceptionally low abundances 
##  when more  than ten year data have been collected (i.e., years with extraordinary sitations that may 
##  indicate noise, rather than a phenological response to climate (see eg. transect 68 year 2012)

## filter by presence/absence ratio
pieris_cbms_pres <- pieris_pres_abs |> 
  group_by(id_transect) |> 
  count(pres_abs) |>
  mutate(total_counts_transect = sum(n),
         prop_pres = n/total_counts_transect) |> 
  filter(pres_abs == "present" & prop_pres > 0.5) |> 
  select(id_transect, prop_pres)


sites_selection <- unique(pieris_pres_abs$id_transect)

for(transect_i in sites_selection) {
  pieris_selected_sites_i <- pieris_selected_sites |> 
    filter(id_transect == transect_i)
  ggplot_pieris_i <- ggplot(data = pieris_selected_sites_i, 
                            aes(x = doy,
                                y = n_indiv,
                                color = pres_abs)
  )+
    geom_point()+
    geom_line()+
    labs(x = "Day of Year",
         y = "N individuals observed",
         color = NULL,
         title = paste("Transect", transect_i))+
    facet_wrap(.~year)+
    theme_bw()
  ggsave(plot = ggplot_pieris_i, 
         filename = here(paste0("figures/cbms_transects_pres_abs/cbms_pieris_count_", transect_i,".png")),
         width = 2600,
         height = 2600,
         units = "px")
  
}

## filtering the data set:
pieris_cbms_selection <- pieris_pres_abs |> 
  filter(id_transect %in% sites_selection) |> 
  filter(!(id_transect == 17 & year == 2000),
         !(id_transect == 29 & year == 1997),    
         !(id_transect == 47 & year == 2020),
         !(id_transect == 68 & year == 2012),
         !(id_transect == 89 & year == 2020), #likely COVID-19 restrictions
         !(id_transect == 90 & year %in% c(2007, 2021)),
         !(id_transect == 94 & year == 2020), #likely COVID-19 restrictions
         !(id_transect == 106 & year == 2020), #likely COVID-19 restrictions
         !(id_transect == 107 & year == 2020), #likely COVID-19 restrictions
         !(id_transect == 110 & year == 2020), #likely COVID-19 restrictions
         !(id_transect == 114 & year == 2020), #likely COVID-19 restrictions
         !(id_transect == 116 & year == 2020), #likely COVID-19 restrictions
         !(id_transect == 117 & year == 2020)) #likely COVID-19 restrictions
  


coords_transects <- pieris_cbms_selection |> 
  select(lat, lon) |> 
  group_by(id_transect) |> 
  slice(1) |> 
  ungroup() |> 
  select(lat, lon)

sites <- pieris_cbms_selection |> 
  rename(ID_coords = id_transect)

first_emergence_adult_sites<- sites |> 
  group_by(year) |> 
  slice_min(date) |> 
  mutate(doy = yday(date))

years_sites <- sort(unique(sites$year))

sites_coords <- sites |> 
  group_by(ID_coords) |> 
  slice(1) |> 
  ungroup() |> 
  dplyr::select(lat, lon)

easyclimate_temp_transects <- easyclimate::get_daily_climate(coords = sites_coords,
                                                             climatic_var = c("Tmin", "Tmax"),
                                                             period = 1988:2022,
                                                             output = "df",
                                                             version = 4) |> 
  as_tibble() |> 
  mutate(date = ymd(date))

easyclimate_temp_transects

daily_temperatures <- easyclimate_temp_transects |> 
  mutate(daily_tavg = map2_dbl(.x = Tmin,
                               .y = Tmax,
                               .f = ~mean(c(.x, .y), na.rm = TRUE))) |> 
  mutate(year = year(date))

##daily
daily_rate_preds_cbms <- tibble()

  for(model_i in unique(selected_models_pieris_rapae$model_name)) {
  
  rate_daily_i <- rate_pred(model2ratesum = model_i,
                            temperature = daily_temperatures,
                            fitted_params = selected_models_pieris_rapae,
                            res = "daily")
  daily_rate_preds_cbms <- bind_rows(daily_rate_preds_cbms, rate_daily_i)
}

rate_sum_preds_daily <- daily_rate_preds_cbms |> 
  mutate(rate_pred = case_when(rate_pred < 0 ~ 0, # <- avoid negative summation
                               is.na(rate_pred) ~0, # <- zeros where model yields NA
                               .default = rate_pred), 
         rate_pred = 100*rate_pred) |>  # as rate_dev equals 1/days_dev
  group_by(year, model_name, ID_coords) |> 
  mutate(rate_sum = cumsum(rate_pred), #heat summation
         emergence_lgl = logic_ratesum(rate_sum), #threshold point
         doy = yday(date)) |> 
  summarise(day_of_emergence = sum(emergence_lgl))  #number of days until threshold is reach

daily_preds_ratesum_cbms <- rate_sum_preds_daily |> 
  mutate(time_slot = "daily")

##hourly
daily_temps_doy <- daily_temperatures  |> 
  mutate(JDay = yday(date)) 

  hourly_temps_cbms <- tibble()

for(transect_i in unique(daily_temps_doy$ID_coords)) {
  print(paste0(transect_i, "/", length(unique(daily_temps_doy$ID_coords))))
  daily_temps_doy_i <- daily_temps_doy |> 
    filter(ID_coords == transect_i)
  lat_i <- unique(daily_temps_doy_i$lat)
  hourly_temps_i <- chillR::make_hourly_temps(latitude = lat_i,
                                              year_file = daily_temps_doy_i,
                                              keep_sunrise_sunset = TRUE) |> 
    pivot_longer(cols = 13:36,
                 names_to = "hour_of_day",
                 values_to = "temperature") |>  
    rename(tmax = Tmax,
           tmin = Tmin,
           sunrise = Sunrise,
           sunset = Sunset,
           daylength = Daylength,
           doy = JDay
           ) 
  
  hourly_temps_cbms <- bind_rows(hourly_temps_cbms, hourly_temps_i) 
}
  
hourly_temps_cbms

hourly_rate_preds <- tibble()

for(model_i in unique(selected_models_pieris_rapae$model_name)) {
  print(paste0("Predicting phenology from TPC model ", model_i))
  rate_hourly_i <- rate_pred(model2ratesum = model_i,
                             temperature = hourly_temps_cbms,
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
  group_by(year, model_name, ID_coords) |> 
  mutate(rate_sum = cumsum(rate_pred), #heat summation
         emergence_lgl = logic_ratesum(rate_sum), #threshold point
         doy = yday(date)) |> 
  summarise(day_of_emergence = round(sum(emergence_lgl)/24))  #number of days until threshold is reach

hourly_preds_ratesum_cbms <- rate_sum_preds_hourly |> 
  mutate(time_slot = "hourly")

## combine them:
preds_ratesum_cbms <- bind_rows(daily_preds_ratesum_cbms,
                                hourly_preds_ratesum_cbms)
write_rds(preds_ratesum_cbms, file = here("data/preds_ratesum_cbms.rds"))

preds_ratesum_cbms <- readRDS(here("data/preds_ratesum_cbms.rds"))
## and observations:
first_emergence_adult_sites <- sites |> 
  group_by(year, ID_coords) |> 
  slice_min(date) |> 
  mutate(doy = yday(date))

doe_nonlinear_cbms <- ggplot(preds_ratesum_cbms,
                             aes(x = year, 
                                 y = day_of_emergence,
                                 ))+
  geom_point(aes(color = time_slot), 
             alpha = .2)+
  geom_smooth(aes(color = time_slot,
                  fill = time_slot))+
  scale_fill_manual(values = c("#3A848B", "#E85038"))+
  scale_color_manual(values = c("#3A848B", "#E85038"))+
  facet_wrap(~model_name)+
  theme_clean()+
  geom_point(data = first_emergence_adult_sites,
             aes(x = year,
                 y = doy),
             color = "#E9C86B",
             alpha = .2)+

  geom_smooth(data = first_emergence_adult_sites,
              aes(x = year, y  = doy),
              color = "#E9C86B",
              fill = "#E9C86B")+
  labs(x = "Year",
       y = "Day of first adult emergence")
doe_nonlinear_cbms
#gilbert
ggsave(filename = here("figures/doe_nonlinear_cbms.png"),
       width = 2600, height = 2600, units = "px")

ggsave(filename = here("figures/doe_nonlinear_cbms.svg"))


###### by transect -------------------------------------------------------------

for(transect_i in unique(preds_ratesum_cbms$ID_coords)){
  first_emergence_adult_sites_i <- first_emergence_adult_sites |>
    filter(ID_coords == transect_i)
  preds_ratesum_cbms_i <- preds_ratesum_cbms |> 
    filter(ID_coords == transect_i)
  doe_nonlinear_cbms_transect_i <- ggplot(preds_ratesum_cbms_i,
                               aes(x = year, 
                                   y = day_of_emergence,
                               ))+
    geom_point(aes(color = time_slot), 
               alpha = .2)+
    geom_smooth(aes(color = time_slot,
                    fill = time_slot))+
    scale_fill_manual(values = c("#3A848B", "#E85038"))+
    scale_color_manual(values = c("#3A848B", "#E85038"))+
    facet_wrap(~model_name)+
    theme_clean()+
    geom_point(data = first_emergence_adult_sites_i,
               aes(x = year,
                   y = doy),
               color = "#E9C86B",
               alpha = .2)+
    
    geom_smooth(data = first_emergence_adult_sites_i,
                aes(x = year, y  = doy),
                color = "#E9C86B",
                fill = "#E9C86B")+
    labs(x = "Year",
         y = "Day of first adult emergence",
         title = paste("Transect", transect_i))
  
  ggsave(filename = here(paste0("figures/cbms_transects_preds/doe_nonlinear_cbms_", transect_i, ".png")),
         width = 2600, height = 2600, units = "px")
  }
