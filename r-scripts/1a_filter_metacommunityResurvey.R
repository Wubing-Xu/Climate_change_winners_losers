## filter metacommunity-resurvey dataset
## to choose subsets of datasets with number of locations in each year >=4, duration >=10, 
## and matching locations across years to keep sites in same locations

rm(list = ls())
setwd("C:/Dropbox/iDiv/Populations_warming/analyses_202406")

library(tidyverse)
library(vegan)
library(reshape2)
library(magrittr)

mr <- read_csv("data/Metacommunity_Resurvey/metacommunity-survey_communities-standardised.csv")
mr_meta <- read_csv("data/Metacommunity_Resurvey/metacommunity-survey_metadata-standardised.csv")
mr_species <- read_csv("data/Metacommunity_Resurvey/manual_community_species_filled_20240530.csv")

mr <- mr %>% dplyr::select(ID, dataset_id, regional, local,  year, species, species_original, value, metric, unit)


## use the manually checked species names to replace the raw one
mr <- mr %>% left_join(mr_species %>% distinct(dataset_id, species_original, species.new))
# all species are matched
mr %>% filter(is.na(species.new)) %>% distinct(dataset_id, species, species_original, species.new)

mr <- mr %>%
  mutate(species = ifelse(is.na(species.new), species, species.new)) %>%
  dplyr::select(- species.new)


# magalhaes_2020: the baseline survey of 5 sites were performed at 4 years. Use the average year to replace it.
mr %>% filter(dataset_id == "magalhaes_2020") %>% distinct(regional, local, year)
id <- mr$dataset_id == "magalhaes_2020" & mr$year %in% c(2003:2006)
mr$year[id] <- 2005

id <- mr_meta$dataset_id == "magalhaes_2020" & mr_meta$year %in% c(2003:2006)
mr_meta$year[id] <- 2005


## check the sites with different coordinates for same locations sampled in different years (coordintes have some uncertainty)
mr_loc_year <- mr_meta %>% distinct(ID, dataset_id, regional, local, latitude, longitude, year)

mc_loc <- mr_loc_year %>% 
  group_by(ID, dataset_id, regional, local) %>% 
  mutate(latitude_median = median(latitude),
         longitude_median = median(longitude),
         latitude_sd = sd(latitude),
         longitude_sd = sd(longitude)) %>% 
  ungroup()

mc_loc %>% filter(latitude_sd > 0 | longitude_sd > 0) %>% distinct(dataset_id) # 9 datasets with varying coordinates
mc_loc %>% filter(latitude_sd > 0 | longitude_sd > 0) %>% distinct(ID, dataset_id, regional) # 16 regions with varying coordinates

#  for the same site id with varying coordinates across years, choose the sites with uncertainty smaller than 0.01 degrees in latitude and longitude, and use the mean values to replace them  
mr_loc_varying_mean <- mc_loc  %>% 
  filter(latitude_sd > 0 | longitude_sd > 0) %>% 
  filter(abs(latitude - latitude_median) < 0.01 & abs(longitude - longitude_median) < 0.01) %>% 
  group_by(ID, dataset_id, regional, local) %>% 
  mutate(latitude_mean = mean(latitude),
         longitude_mean = mean(longitude),) %>% 
  ungroup() %>% 
  dplyr::select(-c(latitude, longitude, latitude_median:longitude_sd)) %>% 
  rename(latitude = latitude_mean, longitude = longitude_mean)

# combine the sites with no uncertainty in coordinates with sites with small uncertainty
mc_loc_clean <- mc_loc %>% 
  filter(latitude_sd == 0 & longitude_sd == 0) %>% 
  dplyr::select(-c(latitude, longitude, latitude_sd,longitude_sd)) %>% 
  rename(latitude = latitude_median, longitude = longitude_median) %>% 
  bind_rows(mr_loc_varying_mean) %>% 
  arrange(ID)

# keep clean locality and years, and update the coordinates
mr_meta <- mc_loc_clean %>% 
  left_join(mr_meta %>% dplyr::select(-c(latitude, longitude)))

mr <- mr %>% inner_join(mc_loc_clean %>% dplyr::select(-c(latitude, longitude)) )


# check whether some datasets have multiple taxon types and reams in the same region 
mr_meta %>% distinct(dataset_id, regional, realm, taxon) %>% filter(duplicated(.[,c("dataset_id", "regional")]))
mr_meta %>% distinct(dataset_id, realm, taxon) %>% filter(duplicated(.[,c("dataset_id")])) 

# check where some samples have duplicated information in the metadata
mr_meta %>% distinct() %>% filter(duplicated(.[,c("dataset_id","year", "regional" ,"local")])) 
mr_meta %>% distinct() %>% filter(duplicated(.[,c("dataset_id","year", "regional" ,"local")]))  %>% distinct(dataset_id)


# combine community and meta data
mr <- mr %>% left_join(mr_meta %>% 
                         dplyr::select(-c(data_pooled_by_authors, data_pooled_by_authors_comment, alpha_grain_comment,
                                          gamma_bounding_box_comment,  gamma_sum_grains_comment, comment, comment_standardisation)) %>% 
                         distinct(), 
                       by =c("ID", "dataset_id", "year", "regional", "local"))


# remove studies that have been included in other database
mr_meta <- mr_meta %>%
  # remove three studies that have included in InsectChange
  # valtonen_2018 = Hungary moths;  schuch_2011 = Germany Marchand Schuch, magnuson_2020 = LTER NTL Macroinvertebrates
  # the study "willig_2010" has been included in BioTIME (StudyID = 54), but BioTIME doesn't provide coordinates or plotID for this study. keep it here 
  filter(! dataset_id %in% c("valtonen_2018", "schuch_2011", "magnuson_2020")) 

mr <- mr %>% 
  filter(dataset_id %in% mr_meta$dataset_id) %>%
  rename(sample = local)

# remove studies with few sites and short duration
mr_4loc <- mr %>%
  group_by(ID, year) %>%
  mutate(n_samp = n_distinct(sample)) %>%
  filter(n_samp > 3) %>%
  group_by(ID) %>%
  mutate(duration = max(year) - min(year) +1) %>%
  ungroup() %>%
  filter(duration >= 10) %>%
  dplyr::select(!c(n_samp, duration))


###########################
## check abundance information 
# check whether samples within the same study have different abundance metrics and units 
mr_4loc %>%
   distinct(ID, sample, metric, unit) %>% 
  group_by(ID, metric, unit) %>%
  summarise(n_metric = n_distinct(metric),
            n_unit = n_distinct(unit)) %>%
  filter(n_metric > 1 | n_unit > 1)
  
table(mr_4loc$metric) #abundance, cover, density, pa, relative abundance
table(mr_4loc$unit) #count,  ind per 250m2 individuals, per transect, pa, percent cover

mr_4loc %>% filter(metric == "cover") %$% hist(value)
mr_4loc %>% filter(metric == "cover") %>% distinct(ID, dataset_id, regional)

mr_4loc %>% filter(metric == "relative abundance") %$% hist(value)
mr_4loc %>% filter(metric == "relative abundance") %>% distinct(ID, dataset_id, regional)

# check datasets with small min and max abundance values 
mr_study_abundance_range <- mr_4loc %>% 
  filter(metric %in% c("abundance", "density")) %>% 
  group_by(ID, dataset_id, regional, metric) %>%
  summarise(max_abundance = max(range(value)),
            min_abundance = min(range(value)),
            range_abundance = diff(range(value))) 

mr_study_abundance_range %>% filter(max_abundance <= 10 | min_abundance < 1)

mr_4loc %>% filter(ID == 58) %>% ggplot(aes(x = value)) + geom_histogram() + scale_x_log10()
mr_4loc %>% filter(ID == 100) %>% ggplot(aes(x = value)) + geom_histogram()
mr_4loc %>% filter(ID == 121) %>% ggplot(aes(x = value)) + geom_histogram()
mr_4loc %>% filter(ID == 767) %>% ggplot(aes(x = value)) + geom_histogram()
mr_4loc %>% filter(ID == 768) %>% ggplot(aes(x = value)) + geom_histogram()
mr_4loc %>% filter(ID == 769) %>% ggplot(aes(x = value)) + geom_histogram() 
mr_4loc %>% filter(ID == 1401) %>% ggplot(aes(x = value)) + geom_histogram()
mr_4loc %>% filter(ID == 1438) %>% ggplot(aes(x = value)) + geom_histogram() + scale_x_log10()
mr_4loc %>% filter(ID == 1514) %>% ggplot(aes(x = value)) + geom_histogram() + scale_x_log10()
mr_4loc %>% filter(ID == 1515) %>% ggplot(aes(x = value)) + geom_histogram() + scale_x_log10()
mr_4loc %>% filter(ID == 1516) %>% ggplot(aes(x = value)) + geom_histogram() + scale_x_log10()
mr_4loc %>% filter(ID == 1529) %>% ggplot(aes(x = value)) + geom_histogram() + scale_x_log10()
mr_4loc %>% filter(ID == 1531) %>% ggplot(aes(x = value)) + geom_histogram() + scale_x_log10()
mr_4loc %>% filter(ID == 1558) %>% ggplot(aes(x = value)) + geom_histogram() 
mr_4loc %>% filter(ID == 1592) %>% ggplot(aes(x = value)) + geom_histogram()  # exclude this dataset; trees and shrubs presence data
mr_4loc %>% filter(ID == 1593) %>% ggplot(aes(x = value)) + geom_histogram() # exclude this dataset; trees and shrubs presence data
mr_4loc %>% filter(ID == 1629) %>% ggplot(aes(x = value)) + geom_histogram() + scale_x_log10()
mr_4loc %>% filter(ID == 1635) %>% ggplot(aes(x = value)) + geom_histogram()
mr_4loc %>% filter(ID == 1637) %>% ggplot(aes(x = value)) + geom_histogram() 

# note: use metric values as "abundance" and "density" (delete "cover" and "pa") 
# for datasets 1592 and 1593, exclude them in abundance analyses
mr_4loc <- mr_4loc %>% 
  mutate(metric = ifelse(ID %in% c(1592, 1593), "abundance_delete", metric)) 


########################
## the number and locations of sites are somewhat different across years
# filter the dataset to maximum the number of same samples for at least two years within a region
mr_filtered <- NULL
for(i in 1:length(unique(mr_4loc$ID))){
  # perform loop for each study
  study <- mr_4loc %>% 
    filter(ID == unique(mr_4loc$ID)[i])
  
  # get a matrix with row are years and colunmns is samples
  year_sample <- as.matrix(xtabs(~ year + sample, data = unique(study[,c("year","sample")]), sparse=TRUE))
  
  # remove samples that are rarely resurveyed (less than half mean number of years resurveyed across samples)
  sample_nyrs <- colSums(year_sample)
  sample_resurved <- sample_nyrs > 0.5 * mean(sample_nyrs)
  year_sample <- year_sample[, sample_resurved]
  
  # number of co-occurred samples between years
  co_samp <- as.matrix(designdist(year_sample, method = "J", terms= "binary"))
  
  #  find which two years (year-pair) have the maximum number of co-occurred samples and duration
  max_co_samp <- melt(co_samp) %>% 
    as_tibble() %>%
    purrr::set_names("year1","year2","n_samp") %>%
    mutate(year1 = as.numeric(year1),
           year2 = as.numeric(year2),
           duration = year2 - year1 + 1) %>%
    filter(duration > 9 & n_samp > 3)
  
  if(nrow(max_co_samp) == 0) {next}
  
  max_co_samp <-  filter(max_co_samp, n_samp >= 0.9*max(n_samp))
  max_co_samp <- filter(max_co_samp, duration == max(duration))
  max_co_samp <-  filter(max_co_samp, n_samp == max(n_samp))
  
  # keep samples that are shared in the two determined years
  samps_year1 <- year_sample[rownames(year_sample) %in% unlist(max_co_samp[1,1]), ]
  samps_year2 <- year_sample[rownames(year_sample) %in% unlist(max_co_samp[1,2]), ]
  samps_shared <- samps_year1 > 0 & samps_year2 > 0
  year_sample <- year_sample[, samps_shared, drop=FALSE]
  
  # Other years except the two priority years will be compared to the two years. Keep only the years have the same sites with the priority years
  id <- rowSums(year_sample) == ncol(year_sample)
  year_sample <-  year_sample[id, ]
  
  # only keep the selected co-occurred samples and years 
  study_filtered <- study %>% 
    filter(sample %in% colnames(year_sample) & year %in% rownames(year_sample))
  
  mr_filtered <- bind_rows(mr_filtered, study_filtered)
}

ggplot(data = study_filtered, aes(longitude, latitude)) + 
  facet_wrap(~year) +
  geom_point(size = 1, alpha = 0.7) + 
  coord_fixed()

# keep only years with at least 4 samples,
# and keep studies with at least 2 time points and duration >10 years
mr_filtered <- mr_filtered %>%
  group_by(ID, year) %>%
  mutate(n_samp = n_distinct(sample)) %>%
  filter(n_samp > 3) %>%
  group_by(ID) %>%
  mutate(n_samples = n_distinct(sample),
         min_samp = min(n_samp),
         n_years = n_distinct(year),
         duration = max(year) - min(year) +1) %>%
  ungroup() %>%
  filter(n_years >= 2 & duration >= 10)

# check how many studies and their attributes
mr_studies <- mr_filtered %>% distinct(ID, n_samples, n_years, duration) # 331 studies
table(mr_studies$n_samples) # 190 studies >= 10
table(mr_studies$n_years)  # 205 studies with n_years <=3, 45 studies >= 10
table(mr_studies$duration)

# the taxon of dataset perry_2023 should be , rather than Invertebrates. Change it.
mr_filtered <- mr_filtered %>% 
  mutate(taxon = ifelse(dataset_id == "perry_2023", "Phytoplankton", taxon))

# save the filtered data
mr_filtered <- mr_filtered %>% dplyr::select(-(n_samp:duration)) 
save(mr_filtered, file = "data/Metacommunity_Resurvey/metacommunityResurvey_filtered.RDATA")


# locations of the filtered dataset
mr_filtered_locations <- mr_filtered %>% 
  distinct(ID, sample, year, latitude, longitude)

# locations of the full dataset
mr_4loc_locations <- mr_4loc %>% 
  distinct(ID, sample, year, latitude, longitude) %>% 
  left_join(mr_filtered_locations %>%
              mutate(keep = "yes")) %>%
  mutate(keep = ifelse(is.na(keep), "no", keep)) %>% 
  group_by(ID) %>%
  mutate(n_coord = n_distinct(latitude, longitude)) %>% 
  ungroup() %>%
  filter(!is.na(latitude) & n_coord >2)


# plot distributions of samples and indicate which records will be removed
pdf('data/Metacommunity_Resurvey/metacommunityResurvey_filtered.pdf', width = 12, height = 10)
id_study <- unique(mr_4loc_locations$ID)
for(i in 1:length(id_study)){
  study <- mr_4loc_locations %>% 
    filter(ID %in% id_study[i])
  
  p <- ggplot(data = study, aes(longitude, latitude)) + 
    facet_wrap(~year) +
    geom_point(aes(colour = keep), size = 1, alpha = 0.7) + 
    labs(title = id_study[i]) +
    theme(legend.position = "top", legend.text = element_text(size=12)) + 
    coord_fixed() +
    scale_color_manual(values=c("yes" = "deepskyblue", "no" = "coral"))
  
  print(p)
}
dev.off()
