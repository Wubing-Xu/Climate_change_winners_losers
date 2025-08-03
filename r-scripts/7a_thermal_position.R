## extract temperature time series for the community data, 
# and link the temperature with species' observed thermal limits for calculating thermal positions

rm(list = ls())

# load packages
needed_libs <- c("tidyverse", "terra", "sf")

usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1])) {   
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

# load assemblage data and climate time series data
load("data/Combined_assemblages.RDATA")
load("data/climate/CRUTS/cruts_4.07_annual_mean_max_min.RDATA")
load("data/climate/HadISST1/HadISST1_annual_mean_max_min.RDATA")
load("data/Species_thermal_limits.RDATA")

# unwrap climate spatRaster
cruts_tempmean <- lapply(cruts_tempmean, rast)
cruts_tempmax <- lapply(cruts_tempmax, rast)
cruts_tempmin <- lapply(cruts_tempmin, rast)
hadisst_tempmean <- lapply(hadisst_tempmean, rast)
hadisst_tempmax <- lapply(hadisst_tempmax, rast)
hadisst_tempmin <- lapply(hadisst_tempmin, rast)


# locations in the dataset "myers-smith_2019_Herschel" are located in the margin of land and ocean. 
# there are no terrestrial climate data in the described locations. we modify the latitude to that in the neighboring land from 69.58 to 69.50 to extract climate data
id <- which(dat$study_name  %in% c("myers-smith_2019_Herschel", "myers-smith_2019_Komakuk"))
dat$latitude[id] <- 69.50

# locations of each region
dat_loc <- dat %>% 
  distinct(study, sample, latitude, longitude) %>%
  left_join(dat_meta %>% dplyr::select(study, realm))

# add the cells of climate rasters where locations located 
dat_loc_xy <- as.matrix(dat_loc[, c("longitude","latitude")]) 
dat_loc <- dat_loc %>%
  mutate(cell_marine = dat_loc_xy %>% terra::cellFromXY(hadisst_tempmean[[1]], .),
         cell_land = dat_loc_xy %>% terra::cellFromXY(cruts_tempmean[[1]], .)) %>%
  mutate(cell = ifelse(realm == "Marine", cell_marine, cell_land))

# summarize locations in each cell across regions
dat_cell <- dat_loc %>%
  group_by(study, realm, cell) %>%
  mutate(nloc = n_distinct(sample)) %>%
  ungroup() %>%
  dplyr::select(study, realm, cell, latitude, longitude, nloc) %>%
  distinct(study, realm, cell, .keep_all = TRUE)

# add sampled years of each study
dat_cell_yr <- dat_cell %>%
  left_join(dat %>% distinct(study, year), relationship = "many-to-many")

# combinations of realm, year; use it for extracting temperature from climate raster
dat_cell_yr <- dat_cell_yr %>%
  unite("realm_year", "realm", "year", remove = FALSE) %>%
  group_by(realm_year) %>%
  nest(data = c(study, realm, year, cell, latitude, longitude, nloc))


#  a function to get temperature for each cell from a climate raster (depend on realm and year) 
# previous_year: whether to get the temperature at one year before the sampling 
get_temp <- function(x, previous_year = FALSE) {
  if(previous_year) x$year <- x$year - 1
  if(x$realm[1] == "Marine"){
    id <- which(names(hadisst_tempmean) == as.character(x$year[1]))  
    tempmean_ras_obsyr <- hadisst_tempmean[[id]]
    tempmax_ras_obsyr <- hadisst_tempmax[[id]]
    tempmin_ras_obsyr <- hadisst_tempmin[[id]]
  }
  if(x$realm[1] != "Marine"){
    id <- which(names(cruts_tempmean) == as.character(x$year[1]))  
    tempmean_ras_obsyr <- cruts_tempmean[[id]]
    tempmax_ras_obsyr <- cruts_tempmax[[id]]
    tempmin_ras_obsyr <- cruts_tempmin[[id]]
  }
  
  temp_cell_obsyr <- data.frame(tempmean = round(tempmean_ras_obsyr[x$cell][,1], 2), 
                                tempmax = tempmax_ras_obsyr[x$cell][,1], 
                                tempmin = tempmin_ras_obsyr[x$cell][,1])
  return(temp_cell_obsyr)
}

# in the preliminary analyses using codes below, some cells in marine temperature data have no data but community were surveyed in those cells
# Use the temperature in neighboring cells (checked manually later) to replace missing values 
cell_147_42 <- cellFromXY(hadisst_tempmean[[1]], cbind(147.5, -42.5))
cell_146_41 <- cellFromXY(hadisst_tempmean[[1]], cbind(146.5, -41.5))
cell_150_23 <- cellFromXY(hadisst_tempmean[[1]], cbind(150.5, -23.5))
cell_10_58 <- cellFromXY(hadisst_tempmean[[1]], cbind(10.5, 58.5))
cell_10_59 <- cellFromXY(hadisst_tempmean[[1]], cbind(10.5, 59.5))

for(i in 1:length(hadisst_tempmean)){
  hadisst_tempmean[[i]][cell_147_42] <- hadisst_tempmean[[i]][cell_147_42 + 1] 
  hadisst_tempmean[[i]][cell_146_41] <- (hadisst_tempmean[[i]][cell_146_41 - 2] + hadisst_tempmean[[i]][cell_146_41 + 2])/2
  hadisst_tempmean[[i]][cell_150_23] <- hadisst_tempmean[[i]][cell_150_23 + 1]
  hadisst_tempmean[[i]][cell_10_59] <- hadisst_tempmean[[i]][cell_10_58]
  
  hadisst_tempmax[[i]][cell_147_42] <- hadisst_tempmax[[i]][cell_147_42 + 1] 
  hadisst_tempmax[[i]][cell_146_41] <- (hadisst_tempmax[[i]][cell_146_41 - 2] + hadisst_tempmax[[i]][cell_146_41 + 2])/2
  hadisst_tempmax[[i]][cell_150_23] <- hadisst_tempmax[[i]][cell_150_23 + 1]
  hadisst_tempmax[[i]][cell_10_59] <- hadisst_tempmax[[i]][cell_10_58] 
  
  hadisst_tempmin[[i]][cell_147_42] <- hadisst_tempmin[[i]][cell_147_42 + 1] 
  hadisst_tempmin[[i]][cell_146_41] <- (hadisst_tempmin[[i]][cell_146_41 - 2] + hadisst_tempmin[[i]][cell_146_41 + 2])/2
  hadisst_tempmin[[i]][cell_150_23] <- hadisst_tempmin[[i]][cell_150_23 + 1]
  hadisst_tempmin[[i]][cell_10_59] <- hadisst_tempmin[[i]][cell_10_58] 
}

# get the temperature for each combinations of realm, year and cell
dat_cell_yr_temp <- dat_cell_yr %>%
  mutate(temp = map(data, get_temp)) %>%
  unnest(c(data, temp)) %>%
  ungroup() %>%
  dplyr::select(-realm_year)

# to extract climate at one year before the sampling
dat_cell_yr_temp_preyr <- dat_cell_yr %>%
  mutate(temp = map(data, get_temp, previous_year = TRUE)) %>%
  unnest(c(data, temp)) %>%
  ungroup() %>%
  dplyr::select(-realm_year)

# calculate the average temperature at the year and year -1 for main analysis;
# use the temperature at the sampling year for sensitivity analyses
dat_cell_yr_temp <- dat_cell_yr_temp %>%
  left_join(dat_cell_yr_temp_preyr %>% 
              dplyr::select(study, year, cell, tempmean_preyr = tempmean, tempmax_preyr = tempmax, tempmin_preyr = tempmin)) %>%
  mutate(tempmean_sampyr = tempmean,
         tempmax_sampyr = tempmax,
         tempmin_sampyr = tempmin,
         tempmean = (tempmean_preyr + tempmean)/2,
         tempmax = (tempmax_preyr + tempmax)/2,
         tempmin = (tempmin_preyr + tempmin)/2) %>%
  dplyr::select(-c(tempmean_preyr, tempmax_preyr, tempmin_preyr))

# how many locations in each study have no climate data
#  these sites are failed to extract climate data for some cells
dat_cell_yr_temp %>% filter(is.na(tempmean)) %>% distinct(study, realm, cell, latitude, longitude,  nloc)

dat_missing_temp <- dat_cell_yr_temp %>% 
  distinct(study, cell, .keep_all = TRUE) %>%
  group_by(study) %>%
  mutate(sum_nloc = sum(nloc)) %>%
  filter(is.na(tempmean)) %>% 
  mutate(sum_nloc_missing = sum(nloc)) %>% 
  ungroup() %>%
  distinct(study, realm, sum_nloc, sum_nloc_missing) %>%
  mutate(p_nloc_missing = sum_nloc_missing/sum_nloc) %>% 
  left_join(dat_meta %>% select(study, study_name))
# note: 19 studies have few sites that were failed to extract temperatures in the first stage analyses
# very few sites lack temperatures after replacing missing temperatures in raster

# six study have no temperature data for many locations in the first stage analyses
# bt_428, mr_117, mr_407, mr_457, mr_495, mr_553
dat_missing_temp %>% 
  filter(p_nloc_missing > 0.4)

# bt_428
loc_bt_428 <- dat %>% filter(study == "bt_428") %>%
  distinct(study, sample, latitude, longitude)
e <- ext(c(0, 45, 40, 65))
e <- ext(c(5, 15, 55, 62))
plot(hadisst_tempmean[[1]], ext = e)
points(loc_bt_428$longitude ,loc_bt_428$latitude, col = "red", cex = 0.5, pch = 1)

# mr_117
loc_mr_117 <- dat %>% filter(study == "mr_117") %>%
  distinct(study, sample, latitude, longitude)
e <- ext(c(140, 150, -45, -35))
plot(hadisst_tempmean[[1]], ext = e)
points(loc_mr_117$longitude ,loc_mr_117$latitude, col = "red", cex = 0.5, pch = 1)

# mr_407
loc_mr_407 <- dat %>% filter(study == "mr_407") %>%
  distinct(study, sample, latitude, longitude)
e <- ext(c(145, 155, -30, -20))
plot(hadisst_tempmean[[1]], ext = e)
points(loc_mr_407$longitude ,loc_mr_407$latitude, col = "red", cex = 0.5, pch = 1)

# mr_457
loc_mr_457 <- dat %>% filter(study == "mr_457") %>%
  distinct(study, sample, latitude, longitude)
e <- ext(c(140, 150, -45, -35))
plot(hadisst_tempmean[[1]], ext = e)
points(loc_mr_457$longitude ,loc_mr_457$latitude, col = "red", cex = 0.5, pch = 1)

# mr_495
loc_mr_495 <- dat %>% filter(study == "mr_495") %>%
  distinct(study, sample, latitude, longitude)
e <- ext(c(140, 150, -45, -35))
plot(hadisst_tempmean[[1]], ext = e)
points(loc_mr_495$longitude ,loc_mr_495$latitude, col = "red", cex = 0.5, pch = 1)

# mr_553
loc_mr_553 <- dat %>% filter(study == "mr_553") %>%
  distinct(study, sample, latitude, longitude)
e <- ext(c(140, 150, -45, -35))
plot(hadisst_tempmean[[1]], ext = e)
points(loc_mr_553$longitude ,loc_mr_553$latitude, col = "red", cex = 0.5, pch = 1)

# distribution of all locations from dataset edgar_2022  
loc_mr_edgar_2022 <- dat %>% filter(study %in% c("mr_109", "mr_418", "mr_543", "mr_120", "mr_498", 
                                          "mr_556", "mr_407", "mr_457", "mr_117", "mr_495", "mr_553")) %>%
  distinct(study, sample, latitude, longitude)
e <- ext(c(120, 170, -50, -10))
e <- ext(c(140, 150, -45, -35))
plot(hadisst_tempmean[[1]], ext = e)
points(loc_mr_edgar_2022$longitude ,loc_mr_edgar_2022$latitude, col = "red", cex = 0.5, pch = 1)


# check the data of bt_191, which has a temperature decrease of 2.5 degree
dat_cell_yr_temp %>% filter(study == "bt_191")
id <- which(names(hadisst_tempmean) == 1953)
mean(values(hadisst_tempmean[[id]]), na.rm=T)
hadisst_tempmean[[id]][17394]
hadisst_tempmean[[id]][17030]
hadisst_tempmean[[id]][17390]

id <- which(names(hadisst_tempmean) == 1965)
mean(values(hadisst_tempmean[[id]]), na.rm=T)
hadisst_tempmean[[id]][17394]
hadisst_tempmean[[id]][17030]
hadisst_tempmean[[id]][17390]
loc_bt_191 <- dat %>% filter(study == "bt_191") %>%
  distinct(study, latitude, longitude)

e <- ext(c(-72, -65, 40, 47))
e <- ext(c(-90, -30, 20, 70))
e <- ext(c(-90, -60, 35, 55))

par(mfrow = c(2, 1))
plot(hadisst_tempmean[[84]], ext = e) # 84
points(loc_bt_191$longitude ,loc_bt_191$latitude, col = "red", cex = 1, pch = 1)

plot(hadisst_tempmean[[96]], ext = e) # 96
points(loc_bt_191$longitude ,loc_bt_191$latitude, col = "red", cex = 1, pch = 1)

years <- 1870:2022
names(hadisst_tempmean) <- years
brks <- seq(7, 18, by=1) 
x <- subset(hadisst_tempmean, years <= 1965 & years >= 1950) %>% rast()
plot(x, ext = e, breaks=brks, col= terrain.colors(11)) 


# remove NA values
dat_cell_yr_temp <- dat_cell_yr_temp %>% filter(!is.na(tempmean))

# calculate average temperature across sites in each region at each year; use it for displaying climate change
dat_region_year_temp <- dat_cell_yr_temp %>%
  group_by(study, realm, year) %>%
  summarise(tempmean = weighted.mean(tempmean, nloc),
            tempmax = weighted.mean(tempmax, nloc),
            tempmin = weighted.mean(tempmin, nloc),
            tempmean_sampyr = weighted.mean(tempmean_sampyr, nloc),
            tempmax_sampyr = weighted.mean(tempmax_sampyr, nloc),
            tempmin_sampyr = weighted.mean(tempmin_sampyr, nloc)) %>%
  ungroup()


### link temperature in surveyed locations with observed species-level thermal limits to calculate thermal positions

# remove species with only one 10-km occurrence and same estimated upper and lower thermal tolerance. No records removed.
splimit <- splimit %>% 
  filter(nocc_climate > 1) %>%
  filter(tempmean_range > 0)

# combinations of study, year and species
study_yr_sp <- dat %>% 
  distinct(study, year, species, specieskey) %>%
  group_by(study) %>%
  complete(year, nesting(species, specieskey)) %>%
  ungroup()

# remove data with no species-level thermal data and climate data. no data removed.
study_yr_sp <- study_yr_sp %>%
  # remove species with no species-level thermal data
  filter(specieskey %in% splimit$specieskey) %>%
  # remove combinations of study and year having no temperature data
  unite(study_yr, study, year, remove = FALSE) %>%
  filter(study_yr %in% 
           (dat_cell_yr_temp %>%
              unite(study_yr, study, year) %>%
              pull(study_yr))) %>%
  dplyr::select(-study_yr)


## calculate thermal position for each species in each region at each sampled year
dat_spe_year_tempos <- dat_cell_yr_temp %>%
  left_join(study_yr_sp, relationship = "many-to-many") %>%
  left_join(splimit) %>%
  # calculate thermal position for each sampled location
  mutate(tmeanpos = (tempmean - tempmean_low)/(tempmean_high - tempmean_low),
         tmaxpos = (tempmax - tempmax_low)/(tempmax_high - tempmax_low),
         tminpos = (tempmin - tempmin_low)/(tempmin_high - tempmin_low),
         tmeanpos98 = (tempmean - tempmean_low2)/(tempmean_high98 - tempmean_low2),
         tmeanpos_sampyr = (tempmean_sampyr - tempmean_low)/(tempmean_high - tempmean_low)) %>%
  # number of locations that thermal positions are out of 0 - 1 
  mutate(nloc_out_tmeanpos = nloc*(tmeanpos >1 | tmeanpos < 0),
         nloc_out_tmaxpos = nloc*(tmaxpos >1 | tmaxpos < 0),
         nloc_out_tminpos = nloc*(tminpos >1 | tminpos < 0),
         nloc_out_tmeanpos98 = nloc*(tmeanpos98 >1 | tmeanpos98 < 0),
         nloc_out_tmeanpos_sampyr = nloc*(tmeanpos_sampyr >1 | tmeanpos_sampyr < 0)) %>%
  # replace values <0 with 0 and values >1 with 1.
  mutate(tmeanpos = ifelse(tmeanpos<0, 0, ifelse(tmeanpos>1, 1, tmeanpos)),
         tmaxpos = ifelse(tmaxpos<0, 0, ifelse(tmaxpos>1, 1, tmaxpos)),
         tminpos = ifelse(tminpos<0, 0, ifelse(tminpos>1, 1, tminpos)),
         tmeanpos98 = ifelse(tmeanpos98<0, 0, ifelse(tmeanpos98>1, 1, tmeanpos98)),
         tmeanpos_sampyr = ifelse(tmeanpos_sampyr<0, 0, ifelse(tmeanpos_sampyr>1, 1, tmeanpos_sampyr))) %>%
  group_by(study, year, species, specieskey) %>%
  summarise(tmeanpos = weighted.mean(tmeanpos, nloc),
            tmaxpos = weighted.mean(tmaxpos, nloc),
            tminpos = weighted.mean(tminpos, nloc),
            tmeanpos98 = weighted.mean(tmeanpos98, nloc),
            tmeanpos_sampyr = weighted.mean(tmeanpos_sampyr, nloc),
            nloc = sum(nloc),
            nloc_out_tmeanpos = sum(nloc_out_tmeanpos),
            nloc_out_tmaxpos = sum(nloc_out_tmaxpos),
            nloc_out_tminpos = sum(nloc_out_tminpos),
            nloc_out_tmeanpos98 = sum(nloc_out_tmeanpos98),
            nloc_out_tmeanpos_sampyr = sum(nloc_out_tmeanpos_sampyr)) %>%
  ungroup()

# how many locations with thermal position are out of 0 - 1
sum(dat_spe_year_tempos$nloc)
sum(dat_spe_year_tempos$nloc_out_tmeanpos)
sum(dat_spe_year_tempos$nloc_out_tmeanpos)/sum(dat_spe_year_tempos$nloc) # 3.2%
sum(dat_spe_year_tempos$nloc_out_tmaxpos)/sum(dat_spe_year_tempos$nloc) # 2.4%
sum(dat_spe_year_tempos$nloc_out_tminpos)/sum(dat_spe_year_tempos$nloc) # 3.3%
sum(dat_spe_year_tempos$nloc_out_tmeanpos98)/sum(dat_spe_year_tempos$nloc)  # 5.9%
sum(dat_spe_year_tempos$nloc_out_tmeanpos_sampyr)/sum(dat_spe_year_tempos$nloc)  # 3.5%


# save data
save(dat_region_year_temp, dat_spe_year_tempos, file = "data/Assemblages_species_thermal_position.RDATA")
