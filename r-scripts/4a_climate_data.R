## download and organize spatial climate data and climate times series for the land and ocean  

rm(list = ls())

# load packages
packages <- c("tidyverse", "terra", "sdmpredictors", "cruts", "ncdf4", "lubridate", "R.utils")

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


#################
## WorldClim data for land

dir.create("data/climate", showWarnings = FALSE)
dir.create("data/climate/wc2.1_2.5m", showWarnings = FALSE)
dir.create("data/climate/wc2.1_2.5m_tavg", showWarnings = FALSE)

## download worldclim data: use bio1 - long term (1970-2000) average annula mean temperature
download.file("https://geodata.ucdavis.edu/climate/worldclim/2_1/base/wc2.1_2.5m_bio.zip",
              destfile = "data/climate/wc2.1_2.5m_bio.zip")

unzip(zipfile = "data/climate/wc2.1_2.5m_bio.zip", 
      files = c("wc2.1_2.5m_bio_1.tif"), exdir = "data/climate/wc2.1_2.5m", list = FALSE, overwrite = TRUE)

file.remove("data/climate/wc2.1_2.5m_bio.zip")


## monthly temperature
download.file("https://geodata.ucdavis.edu/climate/worldclim/2_1/base/wc2.1_2.5m_tavg.zip",
              destfile = "data/climate/wc2.1_2.5m_tavg.zip")

unzip(zipfile = "data/climate/wc2.1_2.5m_tavg.zip", 
      exdir = "data/climate/wc2.1_2.5m_tavg", list = FALSE, overwrite = TRUE)

file.remove("data/climate/wc2.1_2.5m_tavg.zip")


# calculate long-term average temperature of warmest and coldest month
dir.temp <- paste0("data/climate/wc2.1_2.5m_tavg/wc2.1_2.5m_tavg_", sprintf("%02d", 1:12), ".tif")
tavg <- terra::rast(dir.temp)
wc_templmax <-  max(tavg)
wc_templmin <-  min(tavg)

terra::writeRaster(wc_templmax, file = "data/climate/wc2.1_2.5m/wc2.1_2.5m_templmax.tif")
terra::writeRaster(wc_templmin, file = "data/climate/wc2.1_2.5m/wc2.1_2.5m_templmin.tif")

# remove the unzipped monthely temperature folder
unlink("data/climate/wc2.1_2.5m_tavg", recursive = TRUE)

# read climatic raster: long-term temperature in the land
dir.temp <-paste0("data/climate/wc2.1_2.5m/wc2.1_2.5m_", c("bio_1", "templmin", "templmax"), ".tif")
land_temp <- terra::rast(dir.temp)

# check correlation among mean, maximum and minimum temperature
land_temp_values <- data.frame("tmean" = terra::values(land_temp[[1]]), 
                        "tlmin" = terra::values(land_temp[[2]]), 
                        "tlmax" = terra::values(land_temp[[3]]))
id <- !is.na(land_temp_values[,1])
cor(land_temp_values[id, 1:3], method ="pearson") # r-tmean-tmax = 0.971; R-tmean-tmin = 0.977; r-tmax-tmin = 0.900
cor(land_temp_values[id, 1:3], method ="spearman") # r-tmean-tmax = 0.968; R-tmean-tmin = 0.984; r-tmax-tmin = 0.919



################
## Bio-ORACLE for ocean

# find climate layers for marine data
marine_layers <- list_layers(c("Bio-ORACLE"), marine = TRUE, monthly = FALSE, version = 22)
dim(marine_layers)
head(marine_layers)
id <- grepl(pattern = "temperature", marine_layers$name)
marine_layers <- marine_layers[id, ]
dim(marine_layers)

# download and unzip climate layers: long-term (2000-2014) mena, maximum, minimum temperature at the sea surface
layers_load <- c("BO22_tempmean_ss", "BO22_templtmax_ss", "BO22_templtmin_ss")
load_layers(layers_load, datadir = "data/climate/Bio_ORACLE")

unzip(zipfile = "data/climate/Bio_ORACLE/BO22_tempmean_ss_lonlat.zip", 
      exdir = "data/climate/Bio_ORACLE", list = FALSE, overwrite = TRUE)

unzip(zipfile = "data/climate/Bio_ORACLE/BO22_templtmax_ss_lonlat.zip", 
      exdir = "data/climate/Bio_ORACLE", list = FALSE, overwrite = TRUE)

unzip(zipfile = "data/climate/Bio_ORACLE/BO22_templtmin_ss_lonlat.zip", 
      exdir = "data/climate/Bio_ORACLE", list = FALSE, overwrite = TRUE)

file.remove("data/climate/Bio_ORACLE/BO22_tempmean_ss_lonlat.zip")
file.remove("data/climate/Bio_ORACLE/BO22_templtmax_ss_lonlat.zip")
file.remove("data/climate/Bio_ORACLE/BO22_templtmin_ss_lonlat.zip")


# read climatic raster: long-term temperature in the ocean
dir.temp <-paste0("data/climate/Bio_ORACLE/Present.Surface.Temperature.", c("Mean", "Lt.max", "Lt.min"), ".tif")
ocean_temp <- terra::rast(dir.temp)

# check correlation among mean, maximum and minimum temperature
ocean_temp_values <- data.frame("tmean" = terra::values(ocean_temp[[1]]), 
                               "tlmax" = terra::values(ocean_temp[[2]]), 
                               "tlmin" = terra::values(ocean_temp[[3]]))
id <- !is.na(ocean_temp_values[,1])
cor(ocean_temp_values[id, 1:3], method ="pearson") # # r-tmean-tmax = 0.989; r-tmean-tmin = 0.993; r-tmax-tmin = 0.964
cor(ocean_temp_values[id, 1:3], method ="spearman") # r-tmean-tmax = 0.991; r-tmean-tmin = 0.992; r-tmax-tmin = 0.972


################
# CRU TS 4.07:  monthly high-resolution gridded climates

# download the monthly temperature dataset. resolution = 0.5 degree
dir.create("data/climate/CRUTS", showWarnings = FALSE)
download.file("https://crudata.uea.ac.uk/cru/data/hrg/cru_ts_4.07/cruts.2304141047.v4.07/tmp/cru_ts4.07.1901.2022.tmp.dat.nc.gz",
              destfile = "data/climate/CRUTS/cru_ts4.07.1901.2022.tmp.dat.nc.gz")

# unzipping the dataset
gunzip("data/climate/CRUTS/cru_ts4.07.1901.2022.tmp.dat.nc.gz",
       remove = TRUE, overwrite = TRUE)

# calculate annual mean temperature, temperature of warmest and coldest months
cruts <- terra::rast("data/climate/CRUTS/cru_ts4.07.1901.2022.tmp.dat.nc", subds = "tmp")
years <- 1901:2022

cruts_tempmean <- NULL
cruts_tempmax <- NULL
cruts_tempmin <- NULL
for(i in 1:length(years)){
  cruts_yeari <- cruts[[(12*i - 11):(12*i)]]
  # annual mean
  tmean_yeari <- mean(cruts_yeari, na.rm = TRUE)
  cruts_tempmean <- c(cruts_tempmean, tmean_yeari)
  # temperatfeure of warmest month
  tmax_yeari <- max(cruts_yeari, na.rm = TRUE)
  cruts_tempmax <- c(cruts_tempmax, tmax_yeari)
  # temperature of coldest month
  tmin_yeari <- min(cruts_yeari, na.rm = TRUE)
  cruts_tempmin <- c(cruts_tempmin, tmin_yeari)
}

names(cruts_tempmean) <- as.character(years)
names(cruts_tempmax) <- as.character(years)
names(cruts_tempmin) <- as.character(years)

# file.remove("data/climate/CRUTS/cru_ts4.07.1901.2022.tmp.dat.nc")

# calculate long term average of temperature of annul, warmest and coldest month
cruts_templmean <- subset(cruts_tempmean, years <= 2000 & years >= 1970) %>% 
  terra::rast() %>% mean()

cruts_templmax <- subset(cruts_tempmax, years <= 2000 & years >= 1970) %>% 
  terra::rast() %>% mean()

cruts_templmin <- subset(cruts_tempmin, years <= 2000 & years >= 1970) %>% 
  terra::rast() %>% mean()

# save long term average as raster
terra::writeRaster(cruts_templmean, file = "data/climate/CRUTS/cruts_4.07_1970_2000_templmean.tif", overwrite = TRUE)
terra::writeRaster(cruts_templmax, file = "data/climate/CRUTS/cruts_4.07_1970_2000_templmax.tif", overwrite = TRUE)
terra::writeRaster(cruts_templmin, file = "data/climate/CRUTS/cruts_4.07_1970_2000_templmin.tif", overwrite = TRUE)


# save annual values
cruts_tempmean <- lapply(cruts_tempmean, terra::wrap)
cruts_tempmax <- lapply(cruts_tempmax, terra::wrap)
cruts_tempmin <- lapply(cruts_tempmin, terra::wrap)
save(cruts_tempmean, cruts_tempmax, cruts_tempmin, file = "data/climate/CRUTS/cruts_4.07_annual_mean_max_min.RDATA")


################
# HadISST1: monthly sea surface temperature from 1870 - 2022 in 1 degree resolution

# download the monthly temperature dataset
dir.create("data/climate/HadISST1", showWarnings = FALSE)
download.file("https://www.metoffice.gov.uk/hadobs/hadisst/data/HadISST_sst.nc.gz",
              destfile = "data/climate/HadISST1/HadISST_sst.nc.gz")

# unzipping the dataset
gunzip("data/climate/HadISST1/HadISST_sst.nc.gz",
       remove = TRUE, overwrite = TRUE)

# calculate annual mean temperature, temperature of warmest and coldest months
hadisst <- terra::rast("data/climate/HadISST1/HadISST_sst.nc", subds = "sst")
years <- 1870:2023

hadisst_tempmean <- NULL
hadisst_tempmax <- NULL
hadisst_tempmin <- NULL
for(i in 1:length(years)){
  hadisst_yeari <- hadisst[[(12*i - 11):(12*i)]]
  hadisst_yeari[hadisst_yeari == -1000] <- NA # few outlier values as -1000 degree --> NA
  # annual mean
  tmean_yeari <- mean(hadisst_yeari, na.rm = TRUE)
  hadisst_tempmean <- c(hadisst_tempmean, tmean_yeari)
  # temperatfeure of warmest month
  tmax_yeari <- max(hadisst_yeari, na.rm = TRUE)
  hadisst_tempmax <- c(hadisst_tempmax, tmax_yeari)
  # temperature of coldest month
  tmin_yeari <- min(hadisst_yeari, na.rm = TRUE)
  hadisst_tempmin <- c(hadisst_tempmin, tmin_yeari)
}

names(hadisst_tempmean) <- as.character(years)
names(hadisst_tempmax) <- as.character(years)
names(hadisst_tempmin) <- as.character(years)

# file.remove("data/climate/HadISST1/HadISST_sst.nc")

# calculate long term average of temperature of annul, warmest and coldest month
hadisst_templmean <- subset(hadisst_tempmean, years <= 2014 & years >= 2000) %>% 
  terra::rast() %>% mean()

hadisst_templmax <- subset(hadisst_tempmax, years <= 2014 & years >= 2000) %>% 
  terra::rast() %>% mean()

hadisst_templmin <- subset(hadisst_tempmin, years <= 2014 & years >= 2000) %>% 
  terra::rast() %>% mean()

# save long term average as raster
terra::writeRaster(hadisst_templmean, file = "data/climate/HadISST1/HadISST1_2000_2014_templmean.tif", overwrite=TRUE)
terra::writeRaster(hadisst_templmax, file = "data/climate/HadISST1/HadISST1_2000_2014_templmax.tif", overwrite=TRUE)
terra::writeRaster(hadisst_templmin, file = "data/climate/HadISST1/HadISST1_2000_2014_templmin.tif", overwrite=TRUE)

# save annual values
hadisst_tempmean <- lapply(hadisst_tempmean, terra::wrap)
hadisst_tempmax <- lapply(hadisst_tempmax, terra::wrap)
hadisst_tempmin <- lapply(hadisst_tempmin, terra::wrap)
save(hadisst_tempmean, hadisst_tempmax, hadisst_tempmin, file = "data/climate/HadISST1/HadISST1_annual_mean_max_min.RDATA")


##########################
# compare temperature from different database: worldclim vs. cruts; Bio_ORACLE vs. HadISST

# compare CRUTS and WorldClim data
cruts_templmean <- terra::rast("data/climate/CRUTS/cruts_4.07_1970_2000_templmean.tif") # o.5 degree
cruts_templmax <- terra::rast("data/climate/CRUTS/cruts_4.07_1970_2000_templmax.tif")
cruts_templmin <- terra::rast("data/climate/CRUTS/cruts_4.07_1970_2000_templmin.tif")

wc_templmean <- terra::rast("data/climate/wc2.1_2.5m/wc2.1_2.5m_bio_1.tif") # 0.04166667 degree, or 2.5 minute
wc_templmax <- terra::rast("data/climate/wc2.1_2.5m/wc2.1_2.5m_templmax.tif")
wc_templmin <- terra::rast("data/climate/wc2.1_2.5m/wc2.1_2.5m_templmin.tif")

# resample to same raster tempelates
wc_templmean <- terra::aggregate(wc_templmean, fact = 10)
wc_templmean <- terra::resample(wc_templmean, cruts_templmean, method='bilinear')

wc_templmax <- terra::aggregate(wc_templmax, fact = 10)
wc_templmax <- terra::resample(wc_templmax, cruts_templmax, method='bilinear')

wc_templmin <- terra::aggregate(wc_templmin, fact = 10)
wc_templmin <- terra::resample(wc_templmin, cruts_templmin, method='bilinear')

# compare long term annual mean temperature
land_templmean_compare <- data.frame(wc_templmean = terra::values(wc_templmean)[,1], cruts_templmean = values(cruts_templmean)[,1]) %>%
  drop_na()
cor(land_templmean_compare) # 0.997

ggplot(data = land_templmean_compare[sample(1:nrow(land_templmean_compare), 1000), ]) +
  geom_point(aes(x = wc_templmean, y = cruts_templmean)) +
  geom_abline(intercept = 0, slope = 1, col = "red")

# compare long term annual max temperature
land_templmax_compare <- data.frame(wc_templmax = values(wc_templmax)[,1], cruts_templmax = values(cruts_templmax)[,1]) %>%
  drop_na()
cor(land_templmax_compare) # 0.994

ggplot(data = land_templmax_compare[sample(1:nrow(land_templmax_compare), 1000), ]) +
  geom_point(aes(x = wc_templmax, y = cruts_templmax)) +
  geom_abline(intercept = 0, slope = 1, col = "red")

# compare long term annual min temperature
land_templmin_compare <- data.frame(wc_templmin = values(wc_templmin)[,1], cruts_templmin = values(cruts_templmin)[,1]) %>%
  drop_na()
cor(land_templmin_compare) # 0.997

ggplot(data = land_templmin_compare[sample(1:nrow(land_templmin_compare), 1000), ]) +
  geom_point(aes(x = wc_templmin, y = cruts_templmin)) +
  geom_abline(intercept = 0, slope = 1, col = "red")


# compare Bio_ORACLE and HadISST
hadisst_templmean <- terra::rast("data/climate/HadISST1/HadISST1_2000_2014_templmean.tif") # 1 degree
hadisst_templmax <- terra::rast("data/climate/HadISST1/HadISST1_2000_2014_templmax.tif")
hadisst_templmin <- terra::rast("data/climate/HadISST1/HadISST1_2000_2014_templmin.tif")

bo_templmean <- terra::rast("data/climate/Bio_ORACLE/Present.Surface.Temperature.Mean.tif") # 0.083 degree, or 5 minute
bo_templmax <- terra::rast("data/climate/Bio_ORACLE/Present.Surface.Temperature.Lt.max.tif")
bo_templmin <- terra::rast("data/climate/Bio_ORACLE/Present.Surface.Temperature.Lt.min.tif")

# resample to same raster tempelates
bo_templmean <- terra::aggregate(bo_templmean, fact = 10)
bo_templmean <- terra::resample(bo_templmean, hadisst_templmean, method='bilinear')

bo_templmax <- terra::aggregate(bo_templmax, fact = 10)
bo_templmax <- terra::resample(bo_templmax, hadisst_templmax, method='bilinear')

bo_templmin <- terra::aggregate(bo_templmin, fact = 10)
bo_templmin <- terra::resample(bo_templmin, hadisst_templmin, method='bilinear')

# compare long term annual mean temperature
ocean_templmean_compare <- data.frame(bo_templmean = values(bo_templmean)[,1], hadisst_templmean = values(hadisst_templmean)[,1]) %>%
  drop_na()
cor(ocean_templmean_compare) #0.9996

ggplot(data = ocean_templmean_compare[sample(1:nrow(ocean_templmean_compare), 1000), ]) +
  geom_point(aes(x = bo_templmean, y = hadisst_templmean)) +
  geom_abline(intercept = 0, slope = 1, col = "red")

# compare long term annual max temperature
ocean_templmax_compare <- data.frame(bo_templmax = values(bo_templmax)[,1], hadisst_templmax = values(hadisst_templmax)[,1]) %>%
  drop_na()
cor(ocean_templmax_compare) #0.9990

ggplot(data = ocean_templmax_compare[sample(1:nrow(ocean_templmax_compare), 1000), ]) +
  geom_point(aes(x = bo_templmax, y = hadisst_templmax)) +
  geom_abline(intercept = 0, slope = 1, col = "red")

# compare long term annual min temperature
ocean_templmin_compare <- data.frame(bo_templmin = values(bo_templmin)[,1], hadisst_templmin = values(hadisst_templmin)[,1]) %>%
  drop_na()
cor(ocean_templmin_compare) #0.9994

ggplot(data = ocean_templmin_compare[sample(1:nrow(ocean_templmin_compare), 1000), ]) +
  geom_point(aes(x = bo_templmin, y = hadisst_templmin)) +
  geom_abline(intercept = 0, slope = 1, col = "red")

