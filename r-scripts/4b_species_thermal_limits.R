## link species occurrences with climate data and calculate observed species thermal limits within its current ranges

rm(list = ls())

# load packages
packages <- c("tidyverse", "terra")

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

# load rasterized occurrences and summary
load("data/GBIF/Species_distributions_gridcells.RDATA")
# load("data/GBIF/rasterized_records/occ_10km_0127999-210914110416597.RDATA")
load("data/GBIF/Species_distribution_summary.RDATA")


spsuma <- spsuma %>%
  # remove species with few records
  filter(!is.na(aoo10) & aoo10 >= 20) 

occ_rasid <- occ_rasid %>%
  filter(specieskey %in% spsuma$specieskey) %>% 
  filter(!is.na(ras10id))

# read land and marine climatic data
dir_land_temp <-paste0("data/climate/wc2.1_2.5m/wc2.1_2.5m_", c("bio_1", "templmax", "templmin"),".tif", sep="")
land_temp <- terra::rast(dir_land_temp)

dir_marine_temp <-paste0("data/climate/Bio_ORACLE/Present.Surface.Temperature.",c("Mean", "Lt.max", "Lt.min"),".tif",sep="")
marine_temp <- terra::rast(dir_marine_temp)

## project and resample climatic data to raster in the resolution of 10 km that was used to calculate range size
behr <- "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=km +no_defs" # behrmann equal area projection
ras <- terra::rast(xmin=-180, xmax=180, ymin=-90, ymax=90, crs="+proj=longlat +datum=WGS84", resolution = 0.01)
ras10 <- terra::project(ras, behr)
res(ras10) <- 10

land_temp_behr <- terra::project(land_temp, y = behr)
land_temp_ras10 <- terra::resample(land_temp_behr, ras10, method = 'bilinear')
# land_temp_ras10 <- readAll(land_temp_ras10)

marine_temp_behr <-terra::project(marine_temp, y = behr)
marine_temp_ras10 <- raster::resample(marine_temp_behr, ras10, method='bilinear')
# marine_temp_ras10 <- readAll(marine_temp_ras10)


# the temperature values for grid cells that has distribution records
ras10id <- pull(occ_rasid, ras10id) %>% unique()

temp_ras10 <- tibble(ras10id = ras10id) %>%
  mutate(land_tempmean = land_temp_ras10[[1]][ras10id][, 1],
         land_tempmax = land_temp_ras10[[2]][ras10id][, 1],
         land_tempmin = land_temp_ras10[[3]][ras10id][, 1],
         marine_tempmean = marine_temp_ras10[[1]][ras10id][, 1],
         marine_tempmax = marine_temp_ras10[[2]][ras10id][, 1],
         marine_tempmin = marine_temp_ras10[[3]][ras10id][, 1])

# link distribution records with climate
occ_rasid10_temp <- occ_rasid %>% 
  inner_join(temp_ras10, by = "ras10id")

# update which habitat (land or ocean) to keep according to the availability of climate data 
# # If 80%> records in one habitat and more than in the other habitat, we assume it as the main habitat
n_climate_records <- occ_rasid10_temp %>% 
  group_by(species, specieskey) %>% 
  summarise(n_climate_total =  sum(!is.na(land_tempmean) | !is.na(marine_tempmean)),
            n_climate_land = sum(!is.na(land_tempmean)),
            n_climate_marine = sum(!is.na(marine_tempmean))) %>% 
  mutate(keep_new = ifelse(n_climate_land/n_climate_total > 0.8  & n_climate_land > n_climate_marine,"land", NA)) %>%
  mutate(keep_new = ifelse(n_climate_marine/n_climate_total > 0.8 & n_climate_marine > n_climate_land, "ocean", keep_new)) 

spsuma <- spsuma %>% 
  left_join(n_climate_records) %>% 
  mutate(keep_new = ifelse(is.na(keep_new), keep, keep_new)) %>%
  relocate(keep_new, .after = keep)

# compare the original and updated keep. output the inconsistent species to check manually
spsuma %>% filter(keep != keep_new)
sum(spsuma$keep != spsuma$keep_new & !is.na(spsuma$keep)) # 136 species
sum(is.na(spsuma$keep))  # 142 species with keep as NA
sum(is.na(spsuma$keep_new))  # 30 species with keep_new as NA
spsuma_check <- spsuma %>% filter(keep != keep_new | is.na(spsuma$keep) | is.na(spsuma$keep_new))
write_csv(spsuma_check, file = "data/GBIF/Species_distribution_summary_check_habitat.csv")
spsuma_check <- read_csv("data/GBIF/Species_distribution_summary_check_habitat_filled.csv")
spsuma <- spsuma %>% 
  left_join(spsuma_check) %>% 
  mutate(keep_new = ifelse(!is.na(keep_new_checked), keep_new_checked, keep_new)) %>% 
  dplyr::select(-keep_new_checked)
table(spsuma$keep_new, useNA = "always")

# choose the land or marine climate based on main habitats used by species
occ_rasid10_temp <- occ_rasid10_temp %>% 
  inner_join(distinct(spsuma,  species, specieskey, keep, keep_new), by = c("species", "specieskey"))

occ_rasid10_temp <- occ_rasid10_temp %>%
  filter(!is.na(keep_new) & keep_new != "both" & keep_new != "uncertain") %>% 
  mutate(tempmean = ifelse(keep_new == "land", land_tempmean, marine_tempmean),
         tempmax = ifelse(keep_new == "land", land_tempmax, marine_tempmax),
         tempmin = ifelse(keep_new == "land", land_tempmin, marine_tempmin)) %>%
  dplyr::select(-(land_tempmean:marine_tempmin))

# remove rows with temperature as NA
sum(is.na(occ_rasid10_temp$tempmean))/nrow(occ_rasid10_temp) # about 3.6% records have no climate data. remove them
occ_rasid10_temp <- occ_rasid10_temp %>% filter(!is.na(tempmean))

save(occ_rasid10_temp, file = "data/Species_distribution_temperature.RDATA")


# calculate the lower and higher limits of realized temperature and also medium values for each species.
#  measured as the 1% and 99% quantile; 
# for alternative, use the 2.5% and 97.5% quantile
splimit <- occ_rasid10_temp %>% 
  group_by(species, specieskey, keep, keep_new) %>%
  summarise(nocc_climate = n(),
            # based on annual mean temperature
            tempmean_low = quantile(tempmean, 0.01, na.rm=TRUE),
            tempmean_high = quantile(tempmean, 0.99, na.rm=TRUE),
            tempmean_median = quantile(tempmean, 0.5, na.rm=TRUE),
            # based on temperature of warmest month
            tempmax_low = quantile(tempmax, 0.001, na.rm=TRUE),
            tempmax_high = quantile(tempmax, 0.99, na.rm=TRUE),
            tempmax_median = quantile(tempmax, 0.5, na.rm=TRUE),
            # based on temperature of coldest month
            tempmin_low = quantile(tempmin, 0.01, na.rm=TRUE),
            tempmin_high = quantile(tempmin, 0.99, na.rm=TRUE),
            tempmin_median = quantile(tempmin, 0.5, na.rm=TRUE),
            # use the quantile 2.5% and 97.5% to calculate low and high limit of annual mean temperature
            tempmean_low2 = quantile(tempmean, 0.025, na.rm=TRUE),
            tempmean_high98 = quantile(tempmean, 0.975, na.rm=TRUE)
            ) %>%
  ungroup() 

# add thermal range, 
splimit <- splimit %>%
  mutate(tempmean_range = tempmean_high - tempmean_low,
         tempmax_range = tempmax_high - tempmax_low,
         tempmin_range = tempmin_high - tempmin_low,
         tempmean_range98 = tempmean_high98 - tempmean_low2)

# save data
save(splimit, file = "data/Species_thermal_limits.RDATA")


####################
## comparing thermal limits based on 2.5%/97.5% quantile and 1%/99% quantile

ggplot(splimit[sample(1:nrow(splimit), 500),]) +
  geom_point(aes(x = tempmean_high, y = tempmean_high98, col = log(nocc_climate)), alpha = 0.5) +
  geom_abline(col = "red") + 
  scale_color_gradient(low = "blue", high = "red")

ggplot(splimit[sample(1:nrow(splimit), 500),]) +
  geom_point(aes(x = tempmean_low, y = tempmean_low2, col = log(nocc_climate)), alpha = 0.5) +
  geom_abline(col = "red") + 
  scale_color_gradient(low = "blue", high = "red")

summary(splimit$tempmean_range) # min = 0.5293, Q1=8.74, Q2=12.05, mean=12.88, Q3=16.01,Q4=44.32 
summary(splimit$tempmean_range98) # min = 0.4915, Q1=7.1576, Q2=10.0518, mean=10.8911, Q3=13.4437,Q4=39.9238 

splimit <- splimit %>%
  mutate(diff_trange = tempmean_range - tempmean_range98,
         diff_tlow = tempmean_low2 - tempmean_low,
         diff_thigh = tempmean_high - tempmean_high98)

summary(splimit$diff_trange)
summary(splimit$diff_tlow)
summary(splimit$diff_thigh)

ggplot(splimit) + geom_bin2d(aes(x = nocc_climate, y = diff_trange)) + scale_x_log10()
ggplot(splimit) + geom_bin2d(aes(x = nocc_climate, y = diff_tlow)) + scale_x_log10() 
ggplot(splimit) + geom_bin2d(aes(x = nocc_climate, y = diff_thigh)) + scale_x_log10() 
ggplot(splimit) + geom_bin2d(aes(x = nocc_climate, y = tempmean_range)) + scale_x_log10() # positive relations

ggplot(splimit) + geom_bin2d(aes(x = diff_tlow, y = diff_thigh)) # negative relations
ggplot(splimit) + geom_bin2d(aes(x = tempmean_median, y = diff_tlow))
ggplot(splimit) + geom_bin2d(aes(x = tempmean_median, y = diff_thigh)) # species with intermediate temperature have high differences
