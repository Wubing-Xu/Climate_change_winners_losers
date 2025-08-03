####################################
## Combine GBIF occurrences (raw and rasterized), alpha hulls and estimates of range size from individual jobs run HPC

rm(list = ls())

# Set user dependent working directories
print(getwd())
path2wd <- "/gpfs1/work/wubing/Populations_warming/analyses_202406"
setwd(path2wd)

# load packages
packages <- c("tidyverse","dplyr")

for(x in packages){
  if(!require(x, character.only=TRUE)){
    if(!require(x,character.only=TRUE, lib.loc = "/gpfs1/schlecker/home/wubing/R/x86_64-pc-linux-gnu-library/4.2")){
      install.packages(x, repos = "http://cran.us.r-project.org", lib.loc = "/gpfs1/schlecker/home/wubing/R/x86_64-pc-linux-gnu-library/4.2", dependencies = TRUE)
      require(x, character.only = TRUE, lib.loc = "/gpfs1/schlecker/home/wubing/R/x86_64-pc-linux-gnu-library/4.2")					
    }
  }
}

# load species list 
load("data/GBIF/data_to_get_occurrences_20240603.RDATA")

# combine species metadata and distributions at the resolution of 10-km
occ_raster <- lapply(keys, function(x) mget(load(paste0("data/GBIF/rasterized_records/occ_10km_", x, ".RDATA"))))

ras10 <- occ_raster[[1]][[1]]
occ_rasid <- do.call(bind_rows, lapply(occ_raster, "[[", 2))
spsuma <- do.call(bind_rows, lapply(occ_raster, "[[", 3))

# save the combined results
save(spsuma, file = "data/GBIF/Species_distribution_summary.RDATA")
save(ras10, occ_rasid, file = "data/GBIF/Species_distributions_gridcells.RDATA")
