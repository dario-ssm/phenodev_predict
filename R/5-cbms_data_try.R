# 0. Load ----
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
piepra_raw_count_data_csv <- read_delim(("~/GitHub/dev_rapae/dev_rapae/data_source/peticio_sansegundo.csv")) |> 
  glimpse()


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

##transects with at least ten years
piepra_selected_sites <- piepra_raw_count_data |> 
  filter(id_transect %in% c(1, 5, 8, 9, 10, 11, 12, 13, 18, 19, 23, 24, 29, 33, 
                            34, 40, 42, 43, 51, 52, 58, 59, 60, 64, 67, 68, 69, 
                            75, 76, 77, 79, 80, 85, 88, 89, 90, 92, 94, 96, 98, 
                            106, 107, 109, 110, 113, 114, 115, 116, 120, 121, 
                            123, 125, 127, 128, 155, 168)) |> 
  #and exclude 2020 due to COVID-19 related lack of sampling
  filter(year != 2020) |> 
  group_by(id_transect) |> 
  filter(
    !(id_transect == 5 & year %in% c(2004, 2019)) &
      !(id_transect == 9 & year %in% c(2004, 2009, 2022)) &
      !(id_transect == 11 & year %in% c(1997, 1998, 2003, 2004, 2012)) &
      !(id_transect == 12 & year %in% c(1994, 1995, 2004, 2012, 2018)) &
      !(id_transect == 13 & year %in% c(2003, 2004)) &
      !(id_transect == 18 & year %in% 2021) &
      !(id_transect == 19 & year %in% c(2003, 2004)) &
      !(id_transect == 24 & year %in% c(1998, 2001, 2012)) &
      !(id_transect == 29 & year %in% 1997) &
      !(id_transect == 33 & year %in% c(2003, 2004)) &
      !(id_transect == 40 & year %in% c(2012, 2017 ,2018)) &
      !(id_transect == 68 & year %in% 2012) &
      !(id_transect == 69 & year %in% 2021) &
      !(id_transect == 90 & year %in% 2007) &
      !(id_transect == 96 & year %in% 2021) &
      !(id_transect == 107 & year %in% c(2012, 2017)) &
      !(id_transect == 109 & year %in% 2022) &
      !(id_transect == 110 & year %in% 2012) &
      !(id_transect == 115 & year %in% 2010) &
      !(id_transect == 116 & year %in% 2018) &
      !(id_transect == 125 & year %in% 2021) &
      !(id_transect == 155 & year %in% c(2016, 2022)) &
      !(id_transect == 168 & year %in% 2018)
    ) |> 
  filter(n_indiv >= 1)
   
## and now let's extract climate
load(here("data/daily_tmax_df.RData"))


coords_transects <- piepra_selected_sites |> 
  select(lat, lon) |> 
  group_by(id_transect) |> 
  slice(1) |> 
  ungroup() |> 
  select(lat, lon)
site_1 <- piepra_selected_sites |> 
  filter(id_transect %in% (1:10)) |> 
  rename(ID_coords = id_transect)

first_emergence_adult_site_1 <- site_1 |> 
  group_by(year) |> 
  slice_min(date) |> 
  mutate(doy = yday(date))

years_site1 <- sort(unique(site_1$year))

sites_coords <- site_1 |> 
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
  group_by(year, model_name) |> 
  mutate(rate_sum = cumsum(rate_pred), #heat summation
         emergence_lgl = logic_ratesum(rate_sum), #threshold point
         doy = yday(date)) |> 
  summarise(day_of_emergence = sum(emergence_lgl))  #number of days until threshold is reach

rate_sum_preds_daily

first_emergence_adult_sites <- site_1 |> 
  group_by(year, ID_coords) |> 
  slice_min(date) |> 
  mutate(doy = yday(date))

doe_nonlinear_cbms <- ggplot(rate_sum_preds_daily,
                             aes(x = year, 
                                 y = day_of_emergence))+
  geom_point(color = "#3A848B", 
             alpha = .6)+
  geom_line(color = "#3A848B", 
            alpha = .6,
            linetype = "dashed")+
  geom_smooth(color = "#3A848B",
              fill = "#3A848B")+
  facet_wrap(~model_name)+
  theme_clean()+
  geom_point(data = first_emergence_adult_sites,
             aes(x = year,
                 y = doy),
             color = "#E9C86B",
             alpha = .6)+
  geom_line(data = first_emergence_adult_sites,
            aes(x = year,
                y = doy),
            color = "#E9C86B",
            linetype = "dashed",
            alpha = .6)+
  geom_smooth(data = first_emergence_adult_sites,
              aes(x = year, y  = doy),
              color = "#E9C86B",
              fill = "#E9C86B")+
  labs(x = "Year",
       y = "Day of first adult emergence")
doe_nonlinear_cbms
