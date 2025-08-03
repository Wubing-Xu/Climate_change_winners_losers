## find the continents or oceans that each region (meta-community) located

rm(list = ls())

# Set user dependent working directories
user <- Sys.info()["nodename"]
path2wd <- switch(user,
                  "IDIVNB341" = "C:/Dropbox/iDiv/Populations_warming/analyses_202406",
                  "IDIVTS01" = "H:/wubing/iDiv/Populations_warming/analyses_202406")
setwd(path2wd)

# load packages
needed_libs <- c("tidyverse", "terra", "rworldmap", "rnaturalearth", "oceanmap")

usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1])) {   
    install.packages(p, dep = TRUE)
  }
  require(p, character.only = TRUE)
}

sapply(needed_libs, usePackage)
rm(usePackage)


## load coordinates of studies
load("data/Combined_assemblages.RDATA")

# download and prepare spatial data
ocean <- rnaturalearth::ne_download(scale=10, type="geography_marine_polys", category='physical', load=TRUE)
country <- vect(getMap(resolution='low'))
study_points <- terra::vect(x = data.frame(study = dat_meta$study, lon = dat_meta$cent_long, lat = dat_meta$cent_lat), geom = c("lon", "lat"), crs=crs(country))  

# extract to get names of polygons containing each point
study_continent <- terra::extract(country, study_points) %>% 
  tibble() %>% 
  dplyr::select(ID = id.y, continent = REGION) %>% 
  distinct()

study_ocean <- terra::extract(vect(ocean), study_points) %>% 
  tibble() %>% 
  dplyr::select(ID = id.y, ocean = name_en) %>% 
  distinct()

study_ocean %>% filter(duplicated(.[, "ID"])) # check study points with multiple matching
study_ocean %>% filter(ID %in% c(42, 142, 173, 174, 186, 369))
study_ocean <- study_ocean %>% distinct(ID, .keep_all = TRUE)

region <- dat_meta %>% 
  mutate(continent = as.character(study_continent$continent), # return the continent (7 continent model)
         ocean = study_ocean$ocean) # return the name of ocean or sea

# Correct the continents and oceans manually
bind_rows(filter(region, realm == "Marine" & is.na(ocean)),
          filter(region, realm != "Marine" & is.na(continent))) %>%
  dplyr::select(study:climate, continent, ocean, cent_long, cent_lat) %>% 
  pull(study)


# visual check distributions of marine and non-marine studies 
world = map_data('world')

# oceans
ggplot() +
  geom_polygon(data=world, aes(long, lat, group = group), colour=NA, fill= "darkgray") +
  geom_point(data = filter(region, realm == "Marine"),
             aes(x = cent_long, y = cent_lat, colour = ocean), size = 1, stroke = 1.2)

ggplot() +
  geom_polygon(data=world, aes(long, lat, group = group), colour=NA, fill= "darkgray") +
  geom_point(data = filter(region, realm == "Marine" & is.na(ocean)),
             aes(x = cent_long, y = cent_lat), size = 0.5, stroke = 1.2)

# classify Marine studies into Australian Ocean, East Atlantic, West Atlantic, and Pacific Ocean based on coordinates
region %>% filter(realm == "Marine" & cent_long < 0 & cent_lat < 0)
region %>% filter(realm == "Marine" & cent_long > 100 & cent_lat > 0)
region %>% filter(realm == "Marine" & cent_long < -100 & cent_lat > 25)
region %>% filter(realm == "Marine" & cent_long < -50 & cent_lat < 10)
region %>% filter(realm == "Marine" & cent_long < -50 & cent_lat > 10)

region <- region %>% 
  mutate(ocean = ifelse(cent_long > 100 & cent_lat < 0, "Australian Ocean", ocean),
         ocean = ifelse(cent_long > -40 & cent_long < 50 & cent_lat > 0, "East Atlantic", ocean),
         ocean = ifelse(cent_long < -50 & cent_long  > -90 & cent_lat > 10, "West Atlantic", ocean),
         ocean = ifelse(cent_long < -45 & cent_long  > -55 & cent_lat < 0, "West Atlantic", ocean),
         ocean = ifelse(!ocean %in% c("Australian Ocean", "East Atlantic", "West Atlantic"), "Pacific Ocean", ocean))

table(region$ocean)

ggplot() +
  geom_polygon(data=world, aes(long, lat, group = group), colour=NA, fill= "darkgray") +
  geom_point(data = filter(region, realm == "Marine"),
             aes(x = cent_long, y = cent_lat, colour = ocean), size = 1, stroke = 1.2)


# continents
ggplot() +
  geom_polygon(data=world, aes(long, lat, group = group), colour=NA, fill= "darkgray") +
  geom_point(data = filter(region, realm != "Marine"),
             aes(x = cent_long, y = cent_lat, colour = continent), size = 1, stroke = 1.2, alpha = 0.9) 

ggplot() +
  geom_polygon(data=world, aes(long, lat, group = group), colour=NA, fill= "darkgray") +
  geom_point(data = filter(region, realm != "Marine" & is.na(continent)),
             aes(x = cent_long, y = cent_lat), size = 0.5, stroke = 1.2) 

# terrestrial or freshwater studies to be checked. Correct them manually
region %>% 
  filter(realm != "Marine" & is.na(continent)) %>%
  dplyr::select(study:climate, continent, ocean, cent_long, cent_lat)

region <- region %>% 
  mutate(continent = ifelse(study == "bt_525", "Asia", continent),
         continent = ifelse(study == "mr_1512", "North America", continent),
         continent = ifelse(study == "mr_1513", "North America", continent))

table(region$continent)

ggplot() +
  geom_polygon(data=world, aes(long, lat, group = group), colour=NA, fill= "darkgray") +
  geom_point(data = filter(region, realm != "Marine"),
             aes(x = cent_long, y = cent_lat, colour = continent), size = 1, stroke = 1.2, alpha = 0.9) 


## determine regions and check manually
region <- region %>% 
  mutate(region =  ifelse(realm == "Marine", ocean, continent))

table(region$region, useNA = "always")
region %>% dplyr::select(study, realm, region, ocean, continent, cent_lat, cent_long)
region %>% distinct(realm, region) %>% print(n = 16)

ggplot() +
  geom_polygon(data=world, aes(long, lat, group = group), colour=NA, fill= "darkgray") +
  geom_point(data = region, aes(x = cent_long, y = cent_lat, colour = region), size = 1, stroke = 1.2)


save(region, file = "data/Assemblages_regions.RDATA")
