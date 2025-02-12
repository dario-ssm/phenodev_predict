library(leaflet)
## load pieris rapae only data counting
piepra_raw_count_data_csv <- read_delim("~/GitHub/dev_rapae/dev_rapae/data_source/peticio_sansegundo.csv") 
 
piepra_counts <- piepra_raw_count_data_csv |> 
  rename(id_transect = IDitin,
         species = IDesp,
         year = Any,
         date_obs = Date,
         n_individuals = Nindiv,
         transect_m = longitud,
         lat = LAT,
         lon = LNG) |> 
  mutate(species = "Pieris rapae",
         date_obs = dmy(date_obs),
         doy = yday(date_obs),
         id_transect = as_factor(id_transect))# |> 
  #group_by(id_transect, lon, lat, year) |> 
  #slice_min(doy)


all_transects_csv <- read_delim("~/GitHub/dev_rapae/dev_rapae/data_source/m_visit_sub.csv") |> 
  glimpse()
all_transects_data <- all_transects_csv |> 
  rename(id_obs = ...1,
         id_transect = SITE_ID,
         date_obs = DATE) |> 
  mutate(id_transect = as_factor(id_transect)) |> 
  filter(id_transect != 20001) # <- rare

#merge it with piepra data to obtain true absences
piepra_all_count_data <- right_join(piepra_counts, all_transects_data) |> 
  group_by(id_transect) |> 
  mutate(id_species = "Pieris rapae",
         n_individuals = replace_na(n_individuals, 0),
         lat = na.omit(unique(lat)),
         lon = na.omit(unique(lon)),
         pres_or_abs = if_else(n_individuals > 0,
                               "presence",
                               "absence"))
piepra_all_count_data_lgl <- piepra_all_count_data |> 
  group_by(id_transect, lon, lat, year) |> 
  mutate(change_value = map_dbl(.x = pres_or_abs,
                                .f = ~case_when(.x == "presence" &
                                                  lag(.x, default = "absence") == "absence" ~ 1,
                                                .default = 0)))

piepra_first_no_absent <- piepra_all_count_data_lgl |> 
  filter(change_value == 1) |> 
  group_by(id_transect, lon, lat, year) |> 
  slice_min(doy)
