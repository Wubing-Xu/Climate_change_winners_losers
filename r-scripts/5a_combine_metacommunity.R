## Combine the filtered datasets from four database (BioTIME, RivFishTime, InsectChange, Metacommunity Resurvey) and unify metadata  
## and add information such as how many species in each database have thermal estimates

rm(list = ls())

# Set user dependent working directories
user <- Sys.info()["nodename"]
path2wd <- switch(user,
                  "IDIVNB341" = "C:/Dropbox/iDiv/Populations_warming/analyses_202406",
                  "IDIVTS01" = "H:/wubing/iDiv/Populations_warming/analyses_202406")
setwd(path2wd)

# load packages
needed_libs <- c("tidyverse","vegan", "reshape2", "sf", "rgbif")

usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1])) {   
    install.packages(p, dep = TRUE)
  }
  require(p, character.only = TRUE)
}

sapply(needed_libs, usePackage)
rm(usePackage)


load("data/Metacommunity_Resurvey/metacommunityResurvey_filtered.RDATA")
load("data/BioTIME/BioTIME_4locations_10years_filtered_final.RDATA")
load("data/RivFishTIME/RivFishTIME_4sites_10years_filtered.RDATA")
load("data/InsectChange/InsectChange_filtered.RDATA")
load("data/combined_checklists/splist_gbif.RDATA")
# load("data/GBIF/Species_distribution_summary.RDATA")
load("data/Species_thermal_limits.RDATA")


##############
## BioTIME data

# update species with GBIF names and include number of 10-km records
bt <- bt_filtered %>% 
  # filter(!STUDY_ID %in% c(70, 458:465)) %>% # remove problematic studies 
  left_join(splist.gbif %>% dplyr::select(input.species, species, specieskey), by = c("GENUS_SPECIES" = "input.species")) %>% 
  mutate(species = ifelse(is.na(species), GENUS_SPECIES, species)) %>%
  dplyr::select(-GENUS_SPECIES) %>%
  filter((sum.allrawdata.ABUNDANCE > 0 & !is.na(sum.allrawdata.ABUNDANCE)) | (sum.allrawdata.BIOMASS > 0 & !is.na(sum.allrawdata.BIOMASS))) %>% 
  left_join(splimit %>% dplyr::select(specieskey, keep_new, nocc_climate), by ="specieskey")

# select the needed columns and calculate number of samples and duration
# choose studies at least 4 samples in each year, studies with >=2 years, and duration > 10 years,
bt_input <- bt %>% dplyr::select(STUDY_ID, YEAR, location, species, specieskey, keep_new, nocc_climate, LATITUDE, LONGITUDE,
                                 sum.allrawdata.ABUNDANCE, AB_BIO, ABUNDANCE_TYPE) %>% 
  rename(studyID = STUDY_ID, year = YEAR, sample = location, latitude = LATITUDE, longitude = LONGITUDE, abundance = sum.allrawdata.ABUNDANCE,
         metric = AB_BIO, abundance_type = ABUNDANCE_TYPE) %>% # abundance type-- use: Count, Density; have NA values
  group_by(studyID, year, sample, specieskey) %>%
  mutate(abundance = sum(abundance)) %>%
  distinct() %>%
  group_by(studyID) %>%
  mutate(n_samp = n_distinct(sample),
         n_years = n_distinct(year),
         duration = max(year) - min(year) + 1) %>%
  ungroup() %>%
  filter(n_samp > 3 & n_years >= 2 & duration >= 10)

## check
bt_input %>% filter(studyID == 152) %>% distinct(studyID, longitude, latitude) # all longitude values are 180. remove this study 
bt_input <- bt_input %>% filter(studyID != 152)

# check how many studies and their attributes
bt_studies <- bt_input %>% distinct(studyID, n_samp, n_years, duration) # 86 studies
table(bt_studies$n_samp) # 63 studies >= 10
table(bt_studies$n_years)  # 59 studies <= 3; 13 studies >= 10
table(bt_studies$duration)


## meta data for selected studies from BioTIME
bt_input_meta <- bt %>% 
  rename(studyID = STUDY_ID, year = YEAR, sample = location, study_name = TITLE) %>% 
  unite(col = event, studyID, year, sample, remove=FALSE) %>%
  filter(event %in% (bt_input %>%  unite(col = event, studyID, year, sample) %>% pull(event))) %>%
  distinct(studyID, sample, LATITUDE, LONGITUDE, CENT_LONG, CENT_LAT, REALM, TAXA, CLIMATE, GRAIN_SQ_KM, AREA_SQ_KM, study_name) %>%
  distinct(studyID, sample, .keep_all =TRUE) %>%
  group_by(studyID) %>%
  mutate(min_lat = min(LATITUDE, na.rm = TRUE),
         max_lat = max(LATITUDE, na.rm = TRUE),
         min_long = min(LONGITUDE, na.rm = TRUE),
         max_long = max(LONGITUDE, na.rm = TRUE),
         cent_lat = mean(LATITUDE, na.rm = TRUE),
         cent_long = mean(LONGITUDE, na.rm = TRUE)) %>%
  ungroup()

# update central longitude for several studies spanning 180 degree
bt_input_meta <- bt_input_meta %>% 
  mutate(long.new = ifelse((max_long - min_long) > 180, LONGITUDE, NA)) %>% 
  mutate(long.new = ifelse(long.new < 0, long.new + 180, long.new - 180)) %>% 
  group_by(studyID) %>%
  mutate(cent_long_new = mean(long.new, na.rm = TRUE)) %>% 
  ungroup() %>%
  mutate(cent_long_new = ifelse(cent_long_new < 0, cent_long_new + 180, cent_long_new - 180)) %>%
  # distinct(studyID, CENT_LONG, cent_long, cent_long_new)
  mutate(cent_long = ifelse(!is.na(cent_long_new), cent_long_new, cent_long)) %>%
  dplyr::select(-c(long.new, cent_long_new))


## Calculate the area of samples within studies
# the spatial points
bt_input_extent <- bt_input_meta %>% 
  distinct(studyID, LATITUDE , LONGITUDE, min_long, max_long) %>%
  # update longitude for several studies spanning 180 degree
  mutate(long.new = ifelse((max_long - min_long) > 180, LONGITUDE, NA)) %>% 
  mutate(long.new = ifelse(long.new < 0, long.new + 180, long.new - 180)) %>% 
  mutate(LONGITUDE = ifelse(!is.na(long.new), long.new, LONGITUDE)) %>% 
  # determine spatial objects
  distinct(studyID, LATITUDE , LONGITUDE) %>%
  st_as_sf(coords = c('LONGITUDE', 'LATITUDE'))  %>% 
  st_set_crs(4326) %>%
  group_by(studyID) %>% 
  summarise(geometry = st_combine(geometry))

# areas are defined with convex hull
bt_input_area <- bt_input_extent %>% 
  st_convex_hull() %>% 
  st_area() %>%
  as.numeric()
bt_input_area <- bind_cols(studyID = bt_input_extent$studyID, area = round(bt_input_area/10^6,1)) #the unit of area is km2
summary(bt_input_area$area)

# Compare the area and coordinates between the determined based on sites and the raw values
bt_input_meta_check <-  bt_input_meta %>% 
  distinct(studyID, min_lat, max_lat, min_long, max_long, CENT_LONG, CENT_LAT, cent_lat, cent_long, GRAIN_SQ_KM, AREA_SQ_KM) %>% 
  left_join(bt_input_area)

ggplot(bt_input_meta_check) + geom_point(aes(area, AREA_SQ_KM)) + 
  geom_abline(intercept = 0, slope = 1) +
  scale_x_log10() + scale_y_log10()

ggplot(bt_input_meta_check) + geom_point(aes(cent_lat, CENT_LAT)) + 
  geom_abline(intercept = 0, slope = 1) 

ggplot(bt_input_meta_check) + geom_point(aes(cent_long, CENT_LONG)) + 
  geom_abline(intercept = 0, slope = 1) 

ggplot(bt_input_meta_check) + 
  geom_text(aes(cent_long, cent_lat,  label = studyID, colour = "new")) +
  geom_text(aes(CENT_LONG, cent_lat, label = studyID, colour = "raw")) +
  xlim(-180, 180) + ylim(-75, 75)


# choose the determined area based on sites 
bt_input_area <- bt_input_meta_check %>% 
  mutate(mean_grain_m2 = GRAIN_SQ_KM*10^6,
         sd_grain_m2 = NA,
         extent_km2 = area) %>% 
  mutate(mean_grain_m2 = ifelse(mean_grain_m2 == 0, NA, mean_grain_m2)) %>%
  dplyr::select(studyID, mean_grain_m2, sd_grain_m2, extent_km2)

# the start and end years of studies
bt_input_year <- bt_input %>% 
  group_by(studyID) %>% 
  summarise(start_year = min(year),
            end_year = max(year),
            n_years = n_distinct(year),
            duration = end_year - start_year + 1)

# combined meta data with area and start and end years
bt_input_meta <- bt_input_meta %>% 
  dplyr::select(studyID, taxon = TAXA, realm = REALM, climate = CLIMATE, study_name, min_lat, max_lat, min_long, max_long, cent_lat, cent_long) %>%
  distinct() %>%
  inner_join(bt_input_area) %>%
  inner_join(bt_input_year)

# 7 studies with taxa recorded as "All"
bt_studies_taxon_all <- bt_input_meta %>% filter(taxon == "All") %>% distinct(studyID) %>% pull()
bt %>% filter(STUDY_ID %in% bt_studies_taxon_all) %>% distinct(STUDY_ID, REALM, TAXA, ORGANISMS)

# 10 studies with taxa recorded as "Benthos"
bt_studies_taxon_benthos <- bt_input_meta %>% filter(taxon == "Benthos") %>% distinct(studyID) %>% pull()
bt %>% filter(STUDY_ID %in% bt_studies_taxon_benthos) %>% distinct(STUDY_ID, REALM, TAXA, ORGANISMS)



##########
# metacommunity resurvey data

#  update species with GBIF names and indicate whether estimates of range size were calculated
mr <- mr_filtered %>% 
  rename(input_species = species) %>%
  left_join(splist.gbif %>% dplyr::select(input.species, species, specieskey), by = c("input_species" = "input.species")) %>% 
  mutate(species = ifelse(is.na(species), input_species, species)) %>%
  dplyr::select(-input_species) %>%
  left_join(splimit %>% dplyr::select(specieskey, keep_new, nocc_climate), by ="specieskey")

# select the needed columns and calculate number of samples and duration
# choose studies with at least 4 samples in each year, studies with >=2 years, and duration > 10 years,
mr_input <- mr %>% 
  dplyr::select(studyID = ID, year, sample, species, specieskey, keep_new, nocc_climate, latitude, longitude,
                abundance = value, abundance_type = metric) %>% 
  mutate(metric = ifelse(abundance_type == "cover", "B", ifelse(abundance_type == "pa", "PA", "A"))) %>%
  group_by(studyID, year, sample, specieskey) %>%
  mutate(abundance = sum(abundance)) %>%
  distinct() %>%
  group_by(studyID) %>%
  mutate(n_samp = n_distinct(sample),
         n_years = n_distinct(year),
         duration = max(year) - min(year) + 1) %>%
  ungroup() %>%
  filter(n_samp > 3 & n_years >= 2 & duration >= 10)


# check how many studies and their attributes
mr_studies <- mr_input %>% distinct(studyID, n_samp, n_years, duration) # 331 studies
table(mr_studies$n_samp) # 190 studies >= 10
table(mr_studies$n_years)  # 169 studies with n_years = 2, 45 studies >= 10
table(mr_studies$duration)


#### meta data for selected studies
mr_input_meta <- mr %>% 
  unite(col = event, ID, year, sample, remove=FALSE) %>%
  filter(event %in% (mr_input %>%  unite(col = event, studyID, year, sample) %>% pull(event))) %>%
  unite(col = study_name, dataset_id, regional, sep="_",remove=FALSE) %>% 
  rename(studyID = ID) %>% 
  distinct(studyID, study_name, dataset_id, sample, latitude, longitude, realm, taxon, 
           alpha_grain_m2, gamma_extent_km2 = gamma_bounding_box_km2, alpha_grain_type, gamma_extent_type = gamma_bounding_box_type) %>% 
  distinct(studyID, sample, .keep_all =TRUE) %>%
  group_by(studyID) %>%
  mutate(min_lat = min(latitude, na.rm = TRUE),
         max_lat = max(latitude, na.rm = TRUE),
         min_long = min(longitude, na.rm = TRUE),
         max_long = max(longitude, na.rm = TRUE),
         cent_lat = mean(latitude, na.rm = TRUE),
         cent_long = mean(longitude, na.rm = TRUE)) %>%
  ungroup()

# 147 studies have one coordinate (local or central coordinates), and 33 of them have extent > 100 km2, and 34 of them have no extent 
mr_input_meta %>% 
  filter(min_lat == max_lat & min_long == max_long) %>%
  distinct(studyID, study_name, min_lat, max_lat, min_long, max_long)

mr_input_meta %>% 
  filter(min_lat == max_lat & min_long == max_long & gamma_extent_km2 > 100) %>%
  distinct(studyID, study_name, min_lat, max_lat, min_long, max_long)

mr_input_meta %>% 
  filter(min_lat == max_lat & min_long == max_long & is.na(gamma_extent_km2)) %>%
  distinct(studyID, study_name, min_lat, max_lat, min_long, max_long)

# no studies with longitudinal spans >180
mr_input_meta %>% filter(max_long -min_long > 180) %>% distinct(studyID)

## Calculate the area of samples within studies
# the spatial points
mr_input_extent <- mr_input_meta %>% 
  distinct(studyID, latitude , longitude, min_long, max_long) %>%
  # update longitude for several studies spanning 180 degree
  mutate(long.new = ifelse((max_long - min_long) > 180, longitude, NA)) %>% 
  mutate(long.new = ifelse(long.new < 0, long.new + 180, long.new - 180)) %>% 
  mutate(longitude = ifelse(!is.na(long.new), long.new, longitude)) %>% 
  # determine spatial objects
  filter(!is.na(latitude) & !is.na(longitude)) %>%
  distinct(studyID, latitude, longitude) %>%
  st_as_sf(coords = c('longitude', 'latitude'))  %>% 
  st_set_crs(4326) %>%
  group_by(studyID) %>% 
  summarise(geometry = st_combine(geometry))

# areas are defined with convex hull
mr_input_area <- mr_input_extent %>% 
  st_convex_hull() %>% 
  st_area() %>%
  as.numeric()
mr_input_area <-bind_cols(studyID = mr_input_extent$studyID, area = round(mr_input_area/10^6,1)) #the unit of area is km2

# add gamma extent from raw meta data for studies missing values
mr_input_area_inmeta <- mr_input_meta %>% 
  group_by(studyID, study_name, dataset_id) %>% 
  summarise(mean_grain_m2 = mean(alpha_grain_m2, na.rm=TRUE),
            sd_grain_m2 = sd(alpha_grain_m2, na.rm=TRUE),
            extent_km2 = mean(gamma_extent_km2, na.rm=TRUE),
            nsamp = n_distinct(sample))

mr_input_area <- full_join(mr_input_area, mr_input_area_inmeta)

# Compare calculated area and the extent provided by database
ggplot(mr_input_area) + geom_point(aes(area, extent_km2)) + 
  geom_abline(intercept = 0, slope = 1) +
  scale_x_log10() + scale_y_log10()

# Because many studies didn't provide coordinates and thus can't calculated area based on locations, use the extent provided by database
mr_input_area <- mr_input_area %>% 
  mutate(extent_km2 = ifelse(area == 0, extent_km2, area)) %>%
  dplyr::select(-area)

# check whether some studies have no extent
mr_input_area %>% filter(is.na(extent_km2))
mr_input_area %>% filter(is.na(extent_km2)) %>% distinct(dataset_id)
## note: Jandt_2022: it have local coordinates of sites or plots. Multiple plant plots were resurveyed at the same site (coordinate) of some ReSurveyGermany projects. 
#                 spatial extent is small, 1km?
## quimbayo_2022: it have local coordinates of each site. Multiple marine transects (40m2, local) were surveyed at the same site (coordinate).
#                spatial extent is small, 2-4 km?
## rennie_2017_butterflies: it have coordinates of each transect. Each 1-2km long transect was split in up to 15 fixed sections based on habitats. 
#                Local scale is a section of a transect. spatial extent is small, 2 km?
## alves_2022: it has coordinates of each site. Four to ten, 30m*2m belt transects were surveyed for reef community compositions.  
#                Local scale is a belt transect. spatial extent is small, 2-4 km?

# add spatial extent manually based on number of samples and source data information
mr_input_area <- mr_input_area %>% 
  mutate(extent_km2 = ifelse(dataset_id == "quimbayo_2022" & is.na(extent_km2), 0.2*nsamp , extent_km2)) %>% 
  mutate(extent_km2 = ifelse(dataset_id == "rennie_2017_butterflies" & is.na(extent_km2), 2 , extent_km2)) %>% 
  mutate(extent_km2 = ifelse(dataset_id == "alves_2022" & is.na(extent_km2), 0.2*nsamp , extent_km2)) %>% 
  mutate(extent_km2 = ifelse(grepl("jandt_2022", dataset_id) & is.na(extent_km2), round(0.2*sqrt(nsamp), 1) , extent_km2)) 

# the start and end years of studies
mr_input_year <- mr_input %>% 
  group_by(studyID) %>% 
  summarise(start_year = min(year),
            end_year = max(year),
            n_years = n_distinct(year),
            duration = end_year - start_year + 1)

# combined meta data with area and start and end years
mr_input_meta <- mr_input_meta %>% 
  dplyr::select(studyID, study_name, dataset_id, taxon, realm, min_lat, max_lat, min_long, max_long, cent_lat, cent_long) %>%
  distinct() %>%
  inner_join(mr_input_area) %>%
  inner_join(mr_input_year)

# check whether some studies have no central coordinates
mr_input_meta %>% 
  filter(is.na(cent_lat) | is.na(cent_long)) %>% 
  distinct(studyID, study_name, dataset_id)


#########################
## RivFishTIME data

# update species, and determine whether range size data exists
ft <- ft_filtered %>% 
  left_join(splist.gbif %>% dplyr::select(input.species, species, specieskey), by = c("Species" = "input.species")) %>% 
  mutate(species = ifelse(is.na(species), Species, species)) %>%
  dplyr::select(-Species) %>%
  left_join(splimit %>% dplyr::select(specieskey, keep_new, nocc_climate), by ="specieskey")
 
# select the needed columns and calculate number of samples and duration
# choose studies with at least 4 samples in each year, studies with >=2 years, and duration > 10 years,
ft_input <- ft %>% dplyr::select(SourceID, Year, TimeSeriesID, species, specieskey, keep_new, nocc_climate, Latitude, Longitude, Abundance, UnitAbundance) %>%  #, Latitude, Longitude
  rename(studyID = SourceID, year = Year, sample = TimeSeriesID, latitude = Latitude, longitude = Longitude, 
         abundance = Abundance, abundance_type = UnitAbundance) %>% # abundance type-- use: Count, CPUE, Ind.100m2; no NA values
  mutate(metric = "A") %>%
  group_by(studyID, year, sample, specieskey) %>%
  mutate(abundance = sum(abundance)) %>%
  distinct() %>%
  group_by(studyID) %>%
  mutate(n_samp = n_distinct(sample),
         n_years = n_distinct(year),
         duration = max(year) - min(year) + 1) %>%
  ungroup() %>%
  filter(n_samp > 3 & n_years >= 2 & duration >= 10)

# check how many studies and their attributes
ft_studies <- ft_input %>% distinct(studyID, n_samp, n_years, duration) #35 studies
table(ft_studies$n_samp) # # 24 studies >= 10
table(ft_studies$n_years)  # 21 studies <=3; 8 studies >= 10
table(ft_studies$duration)


#### meta data for selected studies from ft
ft_input_meta <- ft %>% 
  rename(studyID = SourceID, year = Year, sample = TimeSeriesID) %>%
  unite(col = event, studyID, year, sample, remove=FALSE) %>%
  filter(event %in% (ft_input %>%  unite(col = event, studyID, year, sample) %>% pull(event))) %>%
  distinct(studyID, sample, Latitude , Longitude) %>%
  filter(!duplicated(.[,c("studyID", "sample")])) %>%
  group_by(studyID) %>%
  mutate(min_lat = min(Latitude, na.rm = TRUE),
         max_lat = max(Latitude, na.rm = TRUE),
         min_long = min(Longitude, na.rm = TRUE),
         max_long = max(Longitude, na.rm = TRUE),
         cent_lat = mean(Latitude, na.rm = TRUE),
         cent_long = mean(Longitude, na.rm =TRUE)) %>%
  ungroup() %>%
  mutate(taxon = "fish",
         realm = "freshwater")

# no study have longitudinal spans >180
ft_input_meta %>% filter(max_long -min_long > 180) %>% distinct(studyID)


## Calculate the area of samples within studies
# the spatial points
ft_input_extent <- ft_input_meta %>% 
  distinct(studyID, Latitude, Longitude) %>%
  st_as_sf(coords = c('Longitude', 'Latitude'))  %>% 
  st_set_crs(4326) %>%
  group_by(studyID) %>% 
  summarise(geometry = st_combine(geometry))

# areas are defined with convex hull
ft_input_area <- ft_input_extent %>% 
  st_convex_hull() %>% 
  st_area() %>%
  as.numeric()
ft_input_area <- bind_cols(studyID = ft_input_extent$studyID, area = round(ft_input_area/10^6,1)) #the unit of area is km2

# add columns in the meta data in other two datasets
ft_input_area <- ft_input_area %>% 
  mutate(mean_grain_m2 =NA,
         sd_grain_m2 =NA,
         extent_km2 = area) %>%
  dplyr::select(-area)
  
# the start and end years of studies
ft_input_year <- ft_input %>% 
  group_by(studyID) %>% 
  summarise(start_year = min(year),
            end_year = max(year),
            n_years = n_distinct(year),
            duration = end_year - start_year + 1)

# combined area with other meta data
ft_input_meta <- ft_input_meta %>% 
  dplyr::select(studyID, taxon, realm, min_lat, max_lat, min_long, max_long, cent_lat, cent_long) %>%
  distinct() %>%
  inner_join(ft_input_area) %>%
  inner_join(ft_input_year)


#########################
## InsectChange dataset

# update species, and determine whether range size data exists
ic <- ic_filtered %>% 
  rename(raw_species = species) %>%
  left_join(splist.gbif %>% dplyr::select(input.species, species, specieskey), by = c("raw_species" = "input.species")) %>% 
  mutate(species = ifelse(is.na(species), raw_species, species)) %>%
  dplyr::select(-raw_species) %>%
  left_join(splimit %>% dplyr::select(specieskey, keep_new, nocc_climate), by ="specieskey") %>%
  filter(Number > 0 & !is.na(Number)) %>% # remove the true absence species
  filter(Datasource_ID != "1533") # remove data set 1533 as abundance of most of records are zero (detected before)

  
# select the needed columns and calculate number of samples and duration
# choose studies with at least 4 samples in each year, studies with >=2 years, and duration > 10 years,
ic_input <- ic %>% dplyr::select(Datasource_ID, Year, Plot_ID, species, specieskey, keep_new, nocc_climate, Latitude, Longitude, 
                                 Number, Abundance.Biomass, Unit.y) %>%  #, Latitude, Longitude
  rename(studyID = Datasource_ID, year = Year, sample = Plot_ID, latitude = Latitude, longitude = Longitude, 
         abundance = Number, metric = Abundance.Biomass, abundance_type = Unit.y) %>% # abundance type-- use: abundance, density; no NA values
  group_by(studyID, year, sample, specieskey) %>%
  mutate(abundance = sum(abundance)) %>%
  distinct() %>%
  group_by(studyID) %>%
  mutate(n_samp = n_distinct(sample),
         n_years = n_distinct(year),
         duration = max(year) - min(year) + 1) %>%
  ungroup() %>%
  filter(n_samp > 3 & n_years >= 2 & duration >= 10)

# check how many studies and their attributes
ic_studies <- ic_input %>% distinct(studyID, n_samp, n_years, duration) # 22 studies
table(ic_studies$n_samp) # 16 studies >= 9
table(ic_studies$n_years) # 9 studies <=3; 6 studies > 10
table(ic_studies$duration)


#### meta data for selected studies from it
ic_input_meta <- ic %>% 
  rename(studyID = Datasource_ID, study_name = Datasource_name, year = Year, sample = Plot_ID) %>%
  unite(col = event, studyID, year, sample, remove=FALSE) %>%
  filter(event %in% (ic_input %>%  unite(col = event, studyID, year, sample) %>% pull(event))) %>%
  distinct(studyID, study_name, sample, Latitude , Longitude, Realm, Climate_zone) %>%
  filter(!duplicated(.[,c("studyID", "sample")])) %>%
  group_by(studyID) %>%
  mutate(min_lat = min(Latitude, na.rm = TRUE),
         max_lat = max(Latitude, na.rm = TRUE),
         min_long = min(Longitude, na.rm = TRUE),
         max_long = max(Longitude, na.rm = TRUE),
         cent_lat = mean(Latitude, na.rm = TRUE),
         cent_long = mean(Longitude, na.rm =TRUE)) %>%
  ungroup() %>%
  mutate(taxon = "invertebrates")

# no study have longitudinal spans >180
ic_input_meta %>% filter(max_long -min_long > 180) %>% distinct(studyID)


## Calculate the area of samples within studies
# the spatial points
ic_input_extent <- ic_input_meta %>% 
  distinct(studyID, Latitude, Longitude) %>%
  st_as_sf(coords = c('Longitude', 'Latitude'))  %>% 
  st_set_crs(4326) %>%
  group_by(studyID) %>% 
  summarise(geometry = st_combine(geometry))

# areas are defined with convex hull
ic_input_area <- ic_input_extent %>% 
  st_convex_hull() %>% 
  st_area() %>%
  as.numeric()
ic_input_area <- bind_cols(studyID = ic_input_extent$studyID, area = round(ic_input_area/10^6,1)) #the unit of area is km2

# add columns in the meta data in other datasets
ic_input_area <- ic_input_area %>% 
  mutate(mean_grain_m2 =NA,
         sd_grain_m2 =NA,
         extent_km2 = area) %>%
  dplyr::select(-area)

# two studies have no extent. Find the extents manually
ic_input_area %>% filter(extent_km2 == 0)

# the area of 1367: a coarse estimate based on map with locations;
# the area of 1408: a coarse estimate based on summed area of catchments
ic_input_area_manul <- data.frame(studyID = c("1367", "1408"),
                                  extent_km2 = c(10, 350))
id <- match(ic_input_area_manul[, "studyID"], pull(ic_input_area, studyID))
ic_input_area[id[!is.na(id)],c("extent_km2")] <- ic_input_area_manul[!is.na(id), c("extent_km2")]


# the start and end years of studies
ic_input_year <- ic_input %>% 
  group_by(studyID) %>% 
  summarise(start_year = min(year),
            end_year = max(year),
            n_years = n_distinct(year),
            duration = end_year - start_year + 1)

# combined area with other meta data
ic_input_meta <- ic_input_meta %>% 
  dplyr::select(studyID, study_name, taxon, realm = Realm, climate = Climate_zone, min_lat, max_lat, min_long, max_long, cent_lat, cent_long) %>%
  distinct() %>%
  inner_join(ic_input_area) %>%
  inner_join(ic_input_year)



#######
## combine datasets from four database: BioTIME, RivFishTIME, InsectChange, Metacommunity Resurvey

# community data
dat <- bind_rows(
  bt_input %>% mutate(database = "bt", studyID = as.character(studyID)),
  ft_input %>% mutate(database = "ft", studyID = as.character(studyID)),
  ic_input %>% mutate(database = "ic", studyID = as.character(studyID)),
  mr_input %>% mutate(database = "mr", studyID = as.character(studyID)),
  ) %>%
  unite(col = "study", database, studyID, remove = FALSE) %>% 
  relocate(database, .before = studyID)

## remove migratory birds 
migrabirds <- read_csv("data/combined_checklists/Birdlife_migratory_birds_2024.06.csv")
migrabirds_gbif <- lapply(1:nrow(migrabirds),function(x) name_backbone(name=migrabirds[x, "scientific_name"], class = "Aves"))
migrabirds_gbif <- bind_rows(migrabirds_gbif)
# migrabirds_gbif <- read_csv("data/combined_checklists/Birdlife_migratory_birds_AddGbifSpeciesNames.csv")
migrabirds_gbif <- migrabirds_gbif %>% dplyr::select(species, speciesKey) %>% distinct()

dat <- dat %>%
  filter(!specieskey %in% migrabirds_gbif$speciesKey)

## add information whether a species have thermal limit data
species_limit <-   splimit %>% 
  filter(nocc_climate >= 20 & !is.na(tempmean_low)) %>%
  select(species, specieskey) 

dat <- dat %>%
  mutate(thermal_limit = ifelse(specieskey %in% species_limit$specieskey, TRUE, FALSE))

# meta data
dat_meta <- bind_rows(
  bt_input_meta %>% 
    mutate(database = "bt",
           studyID = as.character(studyID)),
  ft_input_meta %>% 
    mutate(database = "ft",
           studyID = as.character(studyID)),
  ic_input_meta %>% 
    mutate(database = "ic",
           studyID = as.character(studyID)),
  mr_input_meta %>% 
    mutate(database = "mr",
           studyID = as.character(studyID))) %>%
  relocate(database) %>%
  relocate(climate, .after = realm)  %>%
  mutate(realm = tolower(realm)) %>%
  unite(col = "study", database, studyID, remove = FALSE)

# remove data sets without local coordinates and their spatial extent > 100 km2
dat_meta <- dat_meta %>% 
  filter(!(min_lat == max_lat & min_long == max_long & extent_km2 > 100)) 

dat <- dat %>% filter(study %in% dat_meta$study)

# determine the climate based on central coordinates
dat_meta <- dat_meta %>% 
  mutate(climate = ifelse(abs(cent_lat) <= 23.5, "Tropical", NA),
         climate = ifelse(abs(cent_lat) > 23.5 & abs(cent_lat)<= 60, "Temperate", climate),
         climate = ifelse(abs(cent_lat) > 60, "Polar", climate)) %>%
  relocate(climate, .after = realm)

## indicate whether species level realm information (land or ocean) is consistent with the realm where the community located 
dat <- dat %>%
  left_join(dat_meta %>% dplyr::select(study, realm)) %>% 
  relocate(realm, .after = keep_new) %>%
  mutate(realm_consistent = ifelse((keep_new == "land" & realm %in% c("terrestrial", "freshwater")) | 
                  (keep_new == "ocean" & realm %in% c("marine")), TRUE, FALSE))

# the species with inconsistent information
dat %>% 
  left_join(dat_meta %>% dplyr::select(study, taxon)) %>% 
  filter(!realm_consistent) %>% 
  distinct(species, specieskey, keep_new, realm, realm_consistent, taxon)


# how many species have inconsistent information for a given study and a total of species within the study  
dat %>% 
  distinct(study, specieskey, keep_new, realm, realm_consistent) %>% 
  filter(!is.na(specieskey)) %>%
  group_by(study) %>%
  mutate(n_spp_all = n_distinct(specieskey)) %>%
  filter(!realm_consistent) %>%
  group_by(study, keep_new, realm, n_spp_all) %>%
  summarise(n_spp_realm = n_distinct(specieskey)) %>%
  print(n = 30)
# note: very few species have inconsistent information
# bt_467: marine fish; because some species also distributed in the freshwater, their species-level realm is not certain
# mr_1514: freshwater phytoplankton; many Chromista species also distributed in the ocean.

# determine the proportion of species that have thermal estimates and consistent species level and site-level realm information
data_meta_sprich <- dat %>% 
  group_by(study) %>%
  summarise(sprich = n_distinct(species),
            sprich_gbif = n_distinct(specieskey[!is.na(specieskey)]),
            sprich_thermal = n_distinct(specieskey[thermal_limit & realm_consistent]),
            Psprich_keep = round(sprich_thermal/sprich, 4))

# 52 studies have < 10 species
data_meta_sprich %>% filter(sprich_thermal < 10)
# 9 studies have at least 10 species but less than 50% species were standardized and having thermal limit information
data_meta_sprich %>% filter(sprich_thermal >= 10 & Psprich_keep < 0.5)

# keep data sets with at least 10 species and 50% of all sampled species that have thermal limit information
dat_meta <- dat_meta %>% 
  inner_join(data_meta_sprich) %>% 
  filter(sprich_thermal >= 10 & Psprich_keep >= 0.5)

dat <- dat %>% filter(study %in% dat_meta$study)

#  total number of species: 14,925 species 
dat %>% distinct(species) %>% nrow() 
# number of species that are matched to GBIF backbone: 12,715 species
dat %>% filter(species %in% splist.gbif$species) %>% distinct(species) %>% nrow()
# number of species have thermal estimates and consistent species level and site-level realm information; 11,948 species
dat %>% filter(thermal_limit & realm_consistent) %>% distinct(species) %>% nrow() 


# Update community: select only species useful thermal information
dat <- dat %>% 
  filter(thermal_limit & realm_consistent) %>%
  group_by(study) %>%
  mutate(n_samp = n_distinct(sample),
         n_years = n_distinct(year),
         duration = max(year) - min(year) + 1,
         sprich_keep = n_distinct(species)) %>% 
  # choose studies with number of species >= 10
  filter(n_samp > 3 & n_years >= 2 & duration >= 10 & sprich_keep >= 10) %>%
  group_by(study) %>%
  mutate(period = ifelse(year == min(year), "first", 
                         ifelse(year == max(year), "last", "intermediate"))) %>%
  relocate(period, .after = year) %>%
  ungroup() %>% 
  # add study type
  inner_join(dat_meta %>% dplyr::select(study, taxon, study_name)) %>%
  dplyr::select(-c(thermal_limit, realm_consistent))

# update meta data using updated community data
dat_meta_update <- dat %>% 
  group_by(study) %>%
  summarise(start_year = min(year),
            end_year = max(year),
            n_years = n_distinct(year),
            duration = end_year - start_year + 1,
            n_site = n_distinct(sample))

dat_meta <- dat_meta %>% 
  dplyr::select(study:extent_km2, dataset_id, sprich: Psprich_keep) %>% 
  inner_join(dat_meta_update)


## indicate the subset which can be used for abundance change analyses
dat <- dat %>%
  mutate(subset_abundance = ifelse(database == "bt" & metric != "B" & abundance_type %in% c("Count", "Density"), TRUE, FALSE)) %>%
  mutate(subset_abundance = ifelse(database == "mr" & metric == "A" & abundance_type %in% c("abundance", "density"), TRUE, subset_abundance)) %>%
  mutate(subset_abundance = ifelse(database == "ft" & metric == "A" & abundance_type %in% c("Count", "CPUE", "Ind.100m2"), TRUE, subset_abundance)) %>%
  mutate(subset_abundance = ifelse(database == "ic" & abundance_type %in% c("abundance", "density"), TRUE, subset_abundance))

dat_meta <- dat_meta %>% 
  left_join(dat %>% distinct(study, subset_abundance)) %>%
  relocate(subset_abundance, .after = studyID)

table(dat_meta$subset_abundance) # 244 studies have abundance information


## generalize taxon groups
dat_meta <- dat_meta %>%
  mutate(taxon = tolower(taxon))
table(dat_meta$taxon)

# label taxonomic groups manually
taxon_manual <- data.frame("raw" = c("amphibians", "all", 'benthos', 'birds', 'fish', 
                                     'freshwater invertebrates', 'herpetofauna', 'invertebrates', 'mammals', 'marine invertebrates', 
                                     'marine plants', 'multiple taxa', 'plants', 'terrestrial invertebrates', 'terrestrial plants'),
                           "new" = c("Amphibians and reptiles", "Multiple taxa", "Benthos", "Birds", "Fish", 
                                     "Invertebrates", "Amphibians and reptiles", "Invertebrates", "Mammals", "Invertebrates", 
                                     "Plants", "Multiple taxa", "Plants", "Invertebrates", "Plants")) 
unique(dat_meta$taxon)[!unique(dat_meta$taxon) %in% taxon_manual[,1]]

# add updated taxonomic groups
id <- match(pull(dat_meta, taxon), taxon_manual[,1], nomatch = NA)
table(is.na(id)) #all matched

dat_meta <- dat_meta %>%
  mutate(taxon_new = taxon_manual[id, 2]) %>%
  relocate(taxon_new, .after = taxon)

# transfer the first letter of realms to be capital 
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
dat_meta <- dat_meta %>% mutate(realm = firstup(realm))


# check number of studies across database, study_types,realms, taxonomic groups, climate
table(dat_meta[,c( "database")]) # bt-64, ft-31, ic-21, mr-263 (a total: 379)
table(dat_meta[,c( "realm")]) # terrestrial-183, freshwater-53, marine-143
table(dat_meta[,c("realm", "database")]) 
table(dat_meta[,c("climate", "database")])
table(dat_meta[,c("taxon_new", "realm")])
table(dat_meta[,c("taxon_new", "realm", "database")])
table(dat_meta$n_years == 2) # 193 studies vs. 186; 230 studies <=3
table(dat_meta$n_years < 10) # 321 studies vs. 58 (>= 10 years)


save(dat, dat_meta, file = "data/Combined_assemblages.RDATA")
