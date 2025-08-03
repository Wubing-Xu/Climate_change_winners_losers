## submit GBIF occurrence download requests in PC
# debug r-scripts to download occurrences, and then clean occurrence, which will be run in HPC

rm(list = ls())
setwd("C:/Dropbox/iDiv/Populations_warming/analyses_202406")

library(tidyverse)
library(rgbif)
library(CoordinateCleaner)
library(rnaturalearth)
library(terra)

# GBIF account 
user <- "wubing" # gbif.org username 
pwd <- "xxxxxxxxxx" # gbif.org password
email <- "wbingxu@gmail.com" # email 

# load checklist
load("data/combined_checklists/splist_gbif.RDATA")

# the specieskey from gbif, which are used to download occurrences 
gbif_specieskey <- spgbif %>% pull(specieskey) %>% unique()

## submit download request for GBIF occurrence data
# divide all species into 40 groups, which are used for data requests separately
nsp <- length(gbif_specieskey)
nsp.bin <- ceiling(nsp/40) # number of species in each group
bins <- rep(1:40, times =nsp.bin) 
bins <- bins[1:nsp]

# submit data requests
keys <- vector()

for (i in 1:40){
  x <- occ_download(
    pred_in("taxonKey", gbif_specieskey[bins==i]),
    pred("hasCoordinate", TRUE),
    pred("hasGeospatialIssue", FALSE),
    format = "SIMPLE_CSV",
    user=user, pwd=pwd, email=email
  )
  keys[i] <- x[1]
  print(i)
  Sys.sleep(800)
}

dir.create("data/GBIF", showWarnings = FALSE)
#save(spgbif, gbif_specieskey, buffland, buffocean, keys, file = "data/GBIF/data_to_get_occurrences_20240603.RDATA")
save(spgbif, gbif_specieskey, keys, file = "data/GBIF/data_to_get_occurrences_20240603.RDATA")


####################
## the below lines were just used for debugging in PC
# to formal computation, run the 3b_get_gbif_occurrences_HPC.R in HPC

######################
load("data/GBIF/data_to_get_occurrences_20240603.RDATA")
gbif_specieskey_old <- gbif_specieskey

# new specieskey from gbif
load("data/combined_checklists/splist_gbif.RDATA")

spgbif <- splist.gbif %>% 
  filter(!is.na(specieskey)) %>% 
  dplyr::select(species, specieskey, realm) %>% 
  distinct()

gbif_specieskey <- spgbif %>% pull(specieskey) %>% unique()

# the species key with no estimates of range size in the previous extraction
gbif_specieskey_more <- gbif_specieskey[!gbif_specieskey %in% gbif_specieskey_old]


####################
# get GBIF occurrences for those species that have not been included in previous extraction

# the species that have been included in previous extraction
load("data/gbif/data_to_get_occurrences_20230417.RDATA")
load("data/gbif/Species_distribution_summary.RDATA")

# load species list
load("data/combined_checklists/spgbif_habitat.RDATA")

# all specieskey from gbif
gbif_specieskey <- spgbif %>% pull(specieskey) %>% unique()

# the species key with no estimates of range size in the previous extraction
gbif_specieskey_more <- gbif_specieskey[!gbif_specieskey %in% pull(spsuma, specieskey)]

## submit download request for GBIF occurrence data
# divide all species into several groups, which are used for data requests separately
nsp <- length(gbif_specieskey_more)
nsp.bin <- ceiling(nsp/10) # number of species in each group
bins <- rep(26:35, times =nsp.bin) 
bins <- bins[1:nsp]

# submit data requests
# keys <- vector()
for (i in 26:35){
  x <- occ_download(
    pred_in("taxonKey", gbif_specieskey_more[bins==i]),
    pred("hasCoordinate", TRUE),
    pred("hasGeospatialIssue", FALSE),
    format = "SIMPLE_CSV",
    user=user, pwd=pwd, email=email
  )
  keys[i] <- x[1]
  print(i)
  Sys.sleep(400)
}

save(spgbif, gbif_specieskey, buffland, buffocean, keys, file="data/GBIF/data_to_get_occurrences_20230417.RDATA")


table(pull(spsuma, specieskey) %in% gbif_specieskey)
table(gbif_specieskey %in% pull(spsuma, specieskey))
length(gbif_specieskey) - nrow(spsuma)



#################################
## download and clean GBIF occurrences

# load species list and land/ocean boundary
load("data/GBIF/data_to_get_occurrences_20240603.RDATA")

# GBIF download keys
keys <- occ_download_list(user=user, pwd=pwd, limit=310, start=0)
# key <- keys$results$key[1]
key <- keys$results$key[keys$results$totalRecords <200000 & keys$results$totalRecords >10000][1]

# get occurrences from GBIF
dir.create("data/GBIF/raw_distribution_records", showWarnings = FALSE)
occ <- try(occ_download_get(key, path = "data/GBIF/raw_distribution_records"), silent = TRUE)
occ <- occ_download_import(key=key, path = "data/GBIF/raw_distribution_records")
occ <- occ %>% 
  setNames(tolower(names(.))) %>%
  filter(specieskey %in% gbif_specieskey)

# summary of species information
spsuma_temp <- occ %>% count(specieskey) %>% rename(n_raw = n)
spsuma <- spgbif %>% inner_join(spsuma_temp)

# project coordinates to behrmann equal area projection, and remove records within 1 km grid-cells
behr <- "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=km +no_defs"
xy_behr <- terra::project(as.matrix(occ[, c("decimallongitude","decimallatitude")]), from = "+proj=longlat +datum=WGS84", to = behr)
xy_behr <- as.data.frame(xy_behr) %>% tibble::tibble() %>% setNames(c("x", "y"))
print(7)

occ <- bind_cols(occ, xy_behr) %>% 
  mutate(dplyr::across(c(x, y), round)) %>%
  distinct(specieskey, species, year, x, y, .keep_all = TRUE)

spsuma <- spsuma %>% left_join(occ %>% count(specieskey) %>% rename(n_nodupl = n)) 

# clean gbif occurrences
occ <- occ %>% 
  filter(occurrencestatus  == "PRESENT") %>%
  filter(!basisofrecord %in% c("FOSSIL_SPECIMEN","LIVING_SPECIMEN")) %>%
  filter(!basisofrecord == "MATERIAL_SAMPLE") %>%  #remove *most* metagenomics records 
  filter(!publishingorgkey == "ab733144-7043-4e88-bd4f-fca7bf858880") %>% 
  filter(!taxonrank == "UNRANKED") %>% 
  filter(coordinateuncertaintyinmeters <= 100000 | is.na(coordinateuncertaintyinmeters)) %>%
  filter(!coordinateuncertaintyinmeters %in% c(301,3036,999,9999)) %>%
  CoordinateCleaner::cc_gbif(lon = "decimallongitude", lat = "decimallatitude", buffer = 1000) %>% 
  CoordinateCleaner::clean_coordinates(tests = c("capitals","centroids","equal","institutions","zeros"), 
                                       lon = "decimallongitude", lat = "decimallatitude",
                                       capitals_rad = 1000, centroids_rad = 1000, centroids_detail = "both", value = "clean") %>%
  dplyr::select(specieskey, species, decimallongitude, decimallatitude, year, x, y)

# remove occurrences before 1950
table(is.na(occ$year))
print(sum(is.na(occ$year))/nrow(occ))
occ <- occ %>% filter(year >= 1950)

spsuma <- spsuma %>% left_join(occ %>% count(specieskey) %>% rename(n_noissues = n))  


# flag records in land or sea using buffered land and sea boundary
data("buffland")
data("buffsea")
buffland <- terra::vect(buffland)
buffsea <- terra::vect(buffsea)
crs(buffland) <- crs(buffsea)
occ$inland <- cc_sea(occ, lon = "decimallongitude", lat = "decimallatitude", ref = buffland, value = "flagged")
occ$inocean <- !cc_sea(occ, lon = "decimallongitude", lat = "decimallatitude", ref = buffsea, value = "flagged")

# number of records in buffered land and ocean
spsuma_temp <- occ %>% count(specieskey, inland) %>% filter(inland) %>% dplyr::select(specieskey, n_land = n)
spsuma <- spsuma %>% left_join(spsuma_temp)
spsuma_temp <- occ %>% count(specieskey, inocean) %>% filter(inocean) %>% dplyr::select(specieskey, n_ocean = n)
spsuma <- spsuma %>% left_join(spsuma_temp)

#If 95%> records in one habitat (land or ocean) and more than in the other habitat, we assume it as the main habitat
spsuma <- spsuma %>% 
  mutate(keep_new = ifelse(n_land/n_nodupl > 0.95  & n_land > n_ocean,"land", NA)) %>%
  mutate(keep_new = ifelse(n_ocean/n_nodupl > 0.95 & n_ocean > n_land, "ocean", keep_new)) %>% 
  mutate(keep_new = ifelse(is.na(keep_new), keep, keep_new)) %>%
  relocate(keep_new, .after = keep)

# remove records in sea for land species and remove records in land for sea species
occ <- occ %>% 
  filter(!(!inland & specieskey %in% pull(spsuma[spsuma$keep_new == "land",], specieskey))) %>%
  filter(!(!inocean & specieskey %in% pull(spsuma[spsuma$keep_new == "ocean",], specieskey))) %>% 
  dplyr::select(-c(inland, inocean))

spsuma <- spsuma %>% left_join(occ %>% count(specieskey) %>% rename(n_clean = n))


## count number of occurrences ignoring years, and exclude projected x and y to save storage space
spsuma <- spsuma %>% left_join(
  occ %>% distinct(specieskey, species, x, y) %>% count(specieskey) %>% rename(n_noyear = n))

occ <- occ %>% 
  dplyr::select(-c(x, y))


# save filtered occurrences from GBIF
dir.create("data/GBIF/clean_distribution_records", showWarnings = FALSE)
save(occ, spsuma, file = paste("data/GBIF/clean_distribution_records", paste0("occ_", key, ".RDATA"),sep="/"))


################
## rasterize occurrences to 10 km grid-cells and calculate number of 10-km grid-cells (aoo10)

# raster template in 10-km resolution
ras <- terra::rast(xmin=-180, xmax=180, ymin=-90, ymax=90, crs="+proj=longlat +datum=WGS84", resolution = 0.01)
ras10 <- terra::project(ras, behr)
res(ras10) <- 10

# rasterize occurrences and keep one record in each grid-cell
xy_behr <- terra::project(as.matrix(occ[, c("decimallongitude","decimallatitude")]), from = "+proj=longlat +datum=WGS84", to = behr)

occ_rasid <- occ %>% 
  mutate(ras10id = xy_behr %>% terra::cellFromXY(ras10,.)) %>% 
  dplyr::select(c(species, specieskey, ras10id)) %>% 
  distinct() 

occ_rasid <- occ %>% 
  mutate(ras10id = xy_behr %>% terra::cellFromXY(ras10,.)) %>% 
  group_by(species, specieskey, ras10id) %>% 
  summarise(nocc = n()) %>% 
  ungroup()

## calculate AOO
aoo <- occ_rasid %>% count(specieskey) %>% rename(aoo10 = n)
spsuma <- spsuma %>% left_join(aoo) 

sum(spsuma$n_clean)
sum(occ_rasid$nocc)

spsuma <- spsuma %>% mutate(nratio = n_clean/aoo10)
summary(spsuma$nratio)
colSums(spsuma[5:9], na.rm = TRUE)


# save results
dir.create("data/GBIF/rasterized_records", showWarnings = FALSE)
save(ras10, occ_rasid, spsuma,
     file = paste("data/GBIF/rasterized_records", paste0("occ_10km_", key, ".RDATA"), sep="/"))
