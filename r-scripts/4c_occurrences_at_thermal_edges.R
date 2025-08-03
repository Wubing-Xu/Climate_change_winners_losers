################## to test whether lower data availability toward the equator and thus bias species’ warm limits,
# we compare the number of occurrence records in the warmest and coldest 10% of each species’ observed distribution


rm(list = ls())

# load packages
packages <- c("tidyverse")

for(x in packages){
  if(!require(x, character.only=TRUE)){
    install.packages(x, repos = "http://cran.us.r-project.org", dependencies = TRUE)
    require(x, character.only = TRUE)
  }
}

# Set user dependent working directories
user <- Sys.info()["nodename"]
path2wd <- switch(user,
                  "IDIVNB341" = "C:/Dropbox/iDiv/Populations_warming/analyses_202406",
                  "IDIVTS01" = "H:/wubing/iDiv/Populations_warming/analyses_202406")
setwd(path2wd)

# load temperatures across distribution occurrences
load("data/Species_distribution_temperature.RDATA")

# a self-defined function to calculate percentile of each value within all values 
value_to_percent <- function(x){ 
  fn <- ecdf(x)
  p_x <- fn(x)
  return(p_x)
}

# calculate the percentile (position) of each values from the lowest and the highest,
# and choose the occurrences with temperatures in the lowest 10% and the highest 90%
occ_rasid10_temp_egdes <- occ_rasid10_temp %>% 
  select(-c(ras10id, keep, tempmax, tempmin)) %>% 
  group_by(specieskey) %>% 
  mutate(tempmean_p_low = value_to_percent(tempmean),
         tempmean_p_high = value_to_percent(-tempmean)) %>% 
  ungroup() %>% 
  filter(tempmean_p_low <= 0.1 | tempmean_p_high <= 0.1) %>% 
  mutate(edges = ifelse(tempmean_p_low <= 0.1, "cold", "warm"))

# the number of occurrences in the observed cold and warm distribution edges
nocc_thermal_egdes <- occ_rasid10_temp_egdes %>% 
  group_by(species, specieskey) %>% 
  mutate(nocc_10km = length(nocc)) %>% 
  group_by(species, specieskey, edges, nocc_10km) %>% 
  summarise(nocc = sum(nocc)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = edges, names_prefix = "nocc_", values_from = nocc) 


ggplot(nocc_thermal_egdes) +
  geom_point(aes(x = nocc_cold, y = nocc_warm), alpha = 0.5) +
  scale_x_log10() +
  scale_y_log10()

save(nocc_thermal_egdes, file = "data/species_occurrences_at_thermal_edges.RDATA")
