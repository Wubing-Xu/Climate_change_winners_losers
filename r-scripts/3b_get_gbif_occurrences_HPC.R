# download occurrences, and then clean occurrence
# I submit array jobs in HPC and save clean occurrences for each individual job

rm(list = ls())

# Set user dependent working directories
print(getwd())
path2wd <- "/gpfs1/work/wubing/Populations_warming/analyses_202406"
setwd(path2wd)

# load packages
packages <- c("tidyverse", "rgbif", "terra", "rnaturalearth", "CoordinateCleaner")

for(x in packages){
  if(!require(x, character.only=TRUE)){
    if(!require(x,character.only=TRUE, lib.loc = "/gpfs1/schlecker/home/wubing/R/x86_64-pc-linux-gnu-library/4.2")){
      install.packages(x, repos = "http://cran.us.r-project.org", lib.loc = "/gpfs1/schlecker/home/wubing/R/x86_64-pc-linux-gnu-library/4.2", dependencies = TRUE)
      require(x, character.only = TRUE, lib.loc = "/gpfs1/schlecker/home/wubing/R/x86_64-pc-linux-gnu-library/4.2")					
    }
  }
}

# GBIF account 
print(1)
user <- "wubing" # gbif.org username 
pwd <- "xxxxxxxxxx" # gbif.org password
email <- "wbingxu@gmail.com" # email 

# load species list 
load("data/GBIF/data_to_get_occurrences_20240603.RDATA")

# download keys
## get task id from array jobs submitted in HPC cluster 
task <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
key <- keys[task]
print(task)
print(key)
# key = "0127999-210914110416597"

# A function to unzip large files (>4G): from https://stackoverflow.com/questions/42740206/r-possible-truncation-of-4gb-file
decompress_file <- function(directory, file, .file_cache = FALSE) {
  
  if (.file_cache == TRUE) {
    print("decompression skipped")
  } else {
    
    # Set working directory for decompression
    wd <- getwd()
    setwd(directory)
    
    # Run decompression
    decompression <-
      system2("unzip",
              args = c("-o", # include override flag
                       file),
              stdout = TRUE)
    
    # Reset working directory
    setwd(wd); rm(wd)
    
    # Test for success criteria
    # change the search depending on 
    # your implementation
    if (grepl("Warning message", tail(decompression, 1))) {
      print(decompression)
    }
  }
} 

# A self_defined function of occ_download_import for large files
occ_download_import_self <- function (key = NULL, path = ".", fill = FALSE, 
                                      encoding = "UTF-8", ...) {
  decompress_file(directory = path, file = paste0(key, ".zip"))
  
  targetpath <- sprintf("%s/%s.csv", path, key)
  df <- data.table::fread(targetpath, data.table = FALSE, fill = fill, 
                          encoding = encoding)
  df$countryCode[is.na(df$countryCode)] <- "NA"
  df <- structure(tibble::as_tibble(df), type = "single")
  
  file.remove(targetpath)
  return(df)
}
print(2)

# get occurrences from GBIF
dir.create("data/GBIF/raw_distribution_records", showWarnings = FALSE)
occ <- try(occ_download_get(key, path = "data/GBIF/raw_distribution_records"), silent = TRUE)
print(3)
occ <- occ_download_import_self(key=key, path = "data/GBIF/raw_distribution_records")
print(4)
occ <- occ %>% 
  setNames(tolower(names(.))) %>%
  filter(specieskey %in% gbif_specieskey)
print(5)

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
print(8)

# remove occurrences before 1950
table(is.na(occ$year))
print(sum(is.na(occ$year))/nrow(occ))
occ <- occ %>% filter(year >= 1950)

spsuma <- spsuma %>% left_join(occ %>% count(specieskey) %>% rename(n_clean = n))  

## count number of occurrences ignoring years, and exclude projected x and y to save storage space
spsuma <- spsuma %>% left_join(
  occ %>% distinct(specieskey, species, x, y) %>% count(specieskey) %>% rename(n_noyear = n))

occ <- occ %>% 
  dplyr::select(-c(x, y))

# save filtered occurrences from GBIF
# dir.create("data/GBIF/clean_distribution_records", showWarnings = FALSE)
# save(occ, spsuma, file = paste("data/GBIF/clean_distribution_records", paste0("occ_", key, ".RDATA"),sep="/"))
print(9)


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
  group_by(species, specieskey, ras10id) %>% 
  summarise(nocc = n()) %>% 
  ungroup()


## calculate AOO
aoo <- occ_rasid %>% count(specieskey) %>% rename(aoo10 = n)
spsuma <- spsuma %>% left_join(aoo) 
print(10)

# save results
dir.create("data/GBIF/rasterized_records", showWarnings = FALSE)
save(ras10, occ_rasid, spsuma,
     file = paste("data/GBIF/rasterized_records", paste0("occ_10km_", key, ".RDATA"), sep="/"))
print(11)
