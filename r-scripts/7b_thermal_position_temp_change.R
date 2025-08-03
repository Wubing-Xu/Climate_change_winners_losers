###### calculate changes in thermal position and temperatures between the last and first year of time series, and between 
## the late and early period 

rm(list = ls())

# load packages
needed_libs <- c("tidyverse")

usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1])) {   
    if(p == "mobr") {install_github('MoBiodiv/mobr')}  
    install.packages(p, dep = TRUE)
  }
  require(p, character.only = TRUE)
}

sapply(needed_libs, usePackage)
rm(usePackage)

# Set user dependent working directories
user <- Sys.info()["nodename"]
path2wd <- switch(user,
                  "IDIVNB341" = "C:/Dropbox/iDiv/Populations_warming/analyses_202406",
                  "IDIVTS01" = "H:/wubing/iDiv/Populations_warming/analyses_202406")
setwd(path2wd)

load("data/Assemblages_species_thermal_position.RDATA")
load("data/Species_thermal_limits.RDATA")
load("intermediate_results/Assemblage_occupancy_abundance_change.RDATA")


#### calculate thermal position in the first year and changes of position between the last and first year

# thermal position in the first and last year
tempos_2yr <- dat_spe_year_tempos %>%
  inner_join(occupancy_2yr %>% distinct(study, year, specieskey, period))

# calculate changes in thermal position
tempos_change_2yr <- tempos_2yr %>%
  dplyr::select(study, period, species:nloc) %>%
  pivot_longer(cols = tmeanpos:tmeanpos_sampyr, names_to = "temp_variable", values_to = "position") %>%
  pivot_wider(names_from = period, values_from = position) %>%
  mutate(change = last - first) %>% 
  dplyr::select(-c(nloc, first, last)) %>%
  pivot_wider(names_from = temp_variable, names_prefix = "change_", values_from = change)

# combined thermal position in the first year and changes in thermal position
tempos_first_change_2yr <- tempos_2yr %>%
  filter(period == "first") %>%
  dplyr::select(c(study, species, specieskey, tmeanpos:tmeanpos_sampyr)) %>%
  left_join(tempos_change_2yr)


#### calculate mean thermal position in each period and changes of position between the later and early periods

# mean thermal position in early and late period
tempos_period <- dat_spe_year_tempos %>% 
  left_join(years_period %>% dplyr::select(study, year, period)) %>%
  filter(period %in% c("early", "late")) %>%
  group_by(study, species, specieskey, period) %>%
  summarise(tmeanpos = mean(tmeanpos),
            tmaxpos = mean(tmaxpos),
            tminpos = mean(tminpos),
            tmeanpos98 = mean(tmeanpos98),
            tmeanpos_sampyr = mean(tmeanpos_sampyr)) %>%
  ungroup() 

tempos_period <- tempos_period %>% 
  inner_join(occupancy_period %>% distinct(study, specieskey, period))

# calculate changes in thermal position
tempos_change_period <- tempos_period %>%
  pivot_longer(cols = tmeanpos:tmeanpos_sampyr, names_to = "temp_variable", values_to = "position") %>%
  pivot_wider(names_from = period, values_from = position) %>%
  mutate(change = late - early) %>% 
  dplyr::select(-c(early, late)) %>%
  pivot_wider(names_from = temp_variable, names_prefix = "change_", values_from = change)

# combined thermal position in the early period and changes in thermal position
tempos_early_change_period <- tempos_period %>%
  filter(period == "early") %>%
  dplyr::select(-period) %>%
  left_join(tempos_change_period)


#### calculate changes of temperature between the last and first years
# temperature in the first and last years
temp_2yr <- dat_region_year_temp %>%
  dplyr::select(-c(tempmax_sampyr, tempmin_sampyr)) %>% 
  inner_join(occupancy_2yr %>% distinct(study, year, period))

# changes of temperature
temp_change_2yr <- temp_2yr %>% 
  dplyr::select(-year) %>%
  pivot_longer(cols = tempmean:tempmean_sampyr, names_to = "temp_variable", values_to = "tempe_value") %>%
  pivot_wider(names_from = period, values_from = tempe_value) %>% 
  mutate(temp_change = last - first) %>%
  dplyr::select(-c(first, last)) %>%
  pivot_wider(names_from = temp_variable, names_prefix = "change_", values_from = temp_change)

# combine temperature in the first years and changes of temperature
temp_first_change_2yr <- temp_2yr %>%
  filter(period == "first") %>%
  dplyr::select(-c(year, period)) %>%
  left_join(temp_change_2yr)


#### calculate changes of temperature between the later and early periods
# mean temperature in early and late period
temp_period <- dat_region_year_temp %>%
  left_join(years_period %>% dplyr::select(study, year, period)) %>%
  filter(period %in% c("early", "late")) %>%
  group_by(study, realm, period) %>%
  summarise(tempmean = mean(tempmean),
            tempmax = mean(tempmax),
            tempmin = mean(tempmin),
            tempmean_sampyr = mean(tempmean_sampyr)) %>%
  ungroup()

# temperature change between the late and early periods
temp_change_period <- temp_period %>%
  pivot_longer(cols = tempmean:tempmean_sampyr, names_to = "temp_variable", values_to = "tempe_value") %>%
  pivot_wider(names_from = period, values_from = tempe_value) %>% 
  mutate(temp_change = late - early) %>%
  dplyr::select(-c(early, late)) %>%
  pivot_wider(names_from = temp_variable, names_prefix = "change_", values_from = temp_change)

# combine temperature in the early period and changes of temperature
temp_early_change_period <- temp_period %>%
  filter(period == "early") %>%
  dplyr::select(-period) %>%
  left_join(temp_change_period)


ggplot(temp_change_period) + 
  geom_point(aes(change_tempmean, change_tempmax)) +
  geom_abline(intercept = 0, slope = 1)

ggplot(temp_change_period) + 
  geom_point(aes(change_tempmean, change_tempmean_sampyr)) +
  geom_abline(intercept = 0, slope = 1)

ggplot(temp_change_period) + 
  geom_density(aes(change_tempmean, fill = realm), alpha = 0.5) +
  geom_vline(data = temp_change_period %>% group_by(realm) %>% summarise(change_tempmean = mean(change_tempmean)), 
             aes(xintercept = change_tempmean, col = realm), linewidth =1) +
  geom_vline(xintercept = 0, lty =2)


#### species-level thermal preference and observed thermal range
spe_thermal_range <- splimit %>% 
  dplyr::select(c(species, specieskey, nocc_climate,
                  tempmean_median, tempmax_median, tempmin_median, 
                  tempmean_range, tempmax_range, tempmin_range))

# save data
save(tempos_early_change_period, tempos_first_change_2yr, 
     temp_early_change_period, temp_first_change_2yr, spe_thermal_range, 
     file = "intermediate_results/Assemblage_thermal_position_tempChange.RDATA")
