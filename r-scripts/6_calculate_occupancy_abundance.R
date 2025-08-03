## calculate occupancy in the first and last year and
## the average ones for early and late periods that pool multiple years

rm(list = ls())

# Set user dependent working directories
user <- Sys.info()["nodename"]
path2wd <- switch(user,
                  "IDIVNB341" = "C:/Dropbox/iDiv/Populations_warming/analyses_202406",
                  "IDIVTS01" = "H:/wubing/iDiv/Populations_warming/analyses_202406")
setwd(path2wd)

library(tidyverse)

load("data/Combined_assemblages.RDATA")


# A function to calculate occupancy
get_occupancy <- function(data, fillspecies = TRUE){
  # data should have four columns: study_id, year, sample and species
  colnames(data) <- c("study","year","sample","species")
  
  nsamples <- data %>%
    group_by(study) %>%
    summarise(n_samp = n_distinct(sample)) %>%
    ungroup()
  
   data <- inner_join(data, nsamples) %>%
     group_by(study) %>% 
     filter(n_distinct(year) >= 2) %>% 
     ungroup()
   
   study_meta <- data %>%
     group_by(study) %>%
     summarise(sprich = n_distinct(species),
               #n_years = n_distinct(year),
               start_year = min(year),
               end_year = max(year),
               duration = max(year) - min(year) + 1)
   
     n_occ <- data %>%
       group_by(study, year, species) %>% 
       summarise(n_occ = n_distinct(sample)) %>% 
       ungroup()
     
     if(fillspecies){
       n_occ <- group_by(n_occ, study) %>%
         complete(year, species, fill = list(n_occ = 0)) %>%
         ungroup()
     }
   
   occupancy <- study_meta %>% 
     left_join(nsamples) %>% 
     left_join(n_occ) %>%
     mutate(occupancy = n_occ/n_samp) %>%
     arrange(study, year)
   
   return(occupancy)
}



#######
## Use the first and last years to calculate occupancy and abundance and their changes 

# calculate occupancy
dat_occupancy_2yr <- dat %>% 
  filter(period %in% c("first", "last")) %>% 
  dplyr::select(study, year, sample, specieskey)

occupancy_2yr <- get_occupancy(data = dat_occupancy_2yr, fillspecies = TRUE)

occupancy_2yr <- occupancy_2yr %>% 
  rename(specieskey = species, n_site = n_samp) %>%
  left_join(dat %>% distinct(study, database, studyID, study_name, year, period), by =c("study", "year")) %>%
  left_join(dat %>% distinct(species, specieskey), by =c("specieskey")) %>%
  mutate(n_samp = n_site) %>%
  relocate(database, studyID, study_name, .after = study) %>% 
  relocate(period, .after = year) %>% 
  relocate(species, .after = specieskey) %>% 
  relocate(n_samp, .after = n_site)


# calculate abundance
abundance_2yr <- dat %>% 
  filter(period %in% c("first", "last") & subset_abundance) %>%
  group_by(study, year, specieskey) %>%
  summarise(abundance = sum(abundance)) %>%
  group_by(study) %>%
  complete(year, specieskey, fill = list(abundance = 0)) %>%
  ungroup()

abundance_2yr <- abundance_2yr %>%
  left_join(occupancy_2yr %>% dplyr::select(-c(n_occ, occupancy)) %>% distinct()) %>%
  relocate(study, database, studyID, study_name, year, period, specieskey, species)
  

# calculate occupancy change
occupancy_2yr_occup <- occupancy_2yr %>%
  dplyr::select(-c(year, n_occ)) %>%
  group_by(study, specieskey) %>% 
  pivot_wider(names_from = period, names_prefix = "occup_", values_from = occupancy) %>%
  ungroup()

occupancy_2yr_nocc <- occupancy_2yr %>%
  dplyr::select(-c(year, occupancy)) %>%
  group_by(study, specieskey) %>% 
  pivot_wider(names_from = period, names_prefix = "nocc_", values_from = n_occ) %>%
  ungroup()

occupancy_change_2yr <- occupancy_2yr_occup %>% 
  left_join(occupancy_2yr_nocc) %>% 
  rename(occup_early = occup_first, occup_late = occup_last, nocc_early = nocc_first, nocc_late = nocc_last) %>%
  mutate(dynamic = ifelse(occup_early > 0 & occup_late > 0, "persistent", 
                          ifelse(occup_early > 0 & occup_late == 0, "extinction", "colonization")),
         occup_change = occup_late - occup_early) %>%
  relocate(dynamic, occup_early, occup_late, occup_change, nocc_early, nocc_late, .after = last_col()) %>%
  filter(occup_early > 0 | occup_late > 0) 

table(occupancy_change_2yr[, c("dynamic", "database")])


# calculate abundance change
abundance_change_2yr <- abundance_2yr %>%
  dplyr::select(-c(year)) %>%
  group_by(study, specieskey) %>% 
  pivot_wider(names_from = period, names_prefix = "N_", values_from = abundance) %>%
  ungroup() %>%
  rename(N_early = N_first, N_late = N_last) %>%
  mutate(N_change = N_late - N_early,
         dynamic = ifelse(N_early > 0 & N_late > 0, "persistent", 
                          ifelse(N_early > 0 & N_late == 0, "extinction", "colonization"))) %>%
  filter(N_early > 0 | N_late > 0) # 20 species in bt_176 have 0-0 records



#######
## Use the same number of years before and after the median years to calculate occupancy and abundance and their changes 

## to determine which years will be combined for early and late periods
years_period <- dat %>% 
  distinct(study, duration, year) %>%
  group_by(study) %>%
  mutate(median_year = (max(year) + min(year))/2,
         nyears_early_median =  sum(year < unique(median_year)),
         nyears_late_median =  sum(year > unique(median_year)),
         nyears_combined = min(nyears_early_median, nyears_late_median),
         period_early = year <= sort(year)[nyears_combined],
         period_late = year >= sort(year, decreasing = TRUE)[nyears_combined],
         period = ifelse(period_early, "early", ifelse(period_late, "late", "intermediate"))) %>%
  ungroup()

years_period <- years_period %>% 
  group_by(study) %>%
  mutate(year_period = case_when(period == "early" ~ min(year),
                                 period == "late" ~ max(year),
                                 period == "intermediate" ~ mean(year))) %>%
  group_by(study, period) %>%
  mutate(year_mean_period = mean(year)) %>% 
  group_by(study) %>%
  mutate(duration_mean = max(year_mean_period) - min(year_mean_period) + 1) %>%
  ungroup()

years_period <- years_period %>%
  dplyr::select(study, year, median_year, nyears_combined, period, year_period, year_mean_period, duration_mean)

# add information about year combined to meta data
dat_meta <- dat_meta %>%
  left_join(years_period %>% distinct(study, nyears_combined, duration_mean))


## calculate occupancy
# only use the time points to be pooled to calculate occupancy
dat_occupancy_period <- dat %>% 
  dplyr::select(- period) %>%
  left_join(years_period, by = c("study", "year")) %>% 
  filter(period %in% c("early", "late")) %>% 
  dplyr::select(study, year, sample, specieskey) %>%
  distinct()

occupancy <- get_occupancy(data = dat_occupancy_period, fillspecies = TRUE)

occupancy <- occupancy %>% 
  rename(specieskey = species, n_site = n_samp) %>%
  left_join(dat %>% distinct(study, database, studyID, study_name), by =c("study")) %>% 
  left_join(years_period, by = c("study", "year")) %>%
  left_join(dat %>% distinct(species, specieskey), by =c("specieskey")) %>%
  mutate(n_samp = n_site*nyears_combined) %>%
  dplyr::select(-c(year)) %>%
  relocate(database, studyID, study_name, .after = study) %>% 
  relocate(year_period, period, .after = duration) %>% 
  relocate(species, .after = specieskey) %>% 
  relocate(n_samp, .after = n_site)

# the mean occupancy for the early and late periods
mean_occupancy_period  <- occupancy  %>% 
  group_by(study, specieskey, period) %>%
  summarise(n_occ = sum(n_occ),
            occupancy = mean(occupancy)) %>%
  ungroup()

occupancy_period <- occupancy %>% 
  dplyr::select(-c(n_occ, occupancy)) %>% 
  distinct() %>% 
  left_join(mean_occupancy_period, by = c("study", "specieskey", "period"))


## calculate abundance
abundance_period <- dat %>% 
  filter(subset_abundance) %>% # subset for abundance analyses
  dplyr::select(- period) %>%
  left_join(years_period, by = c("study", "year")) %>% 
  filter(period %in% c("early", "late")) %>% 
  group_by(study, year, period, specieskey) %>%
  summarise(abundance = sum(abundance)) %>%
  group_by(study,period, specieskey) %>%
  summarise(abundance = mean(abundance)) %>%
  group_by(study) %>%
  complete(period, specieskey, fill = list(abundance = 0)) %>%
  ungroup()

abundance_period <- abundance_period %>%
  left_join(occupancy_period %>% dplyr::select(-c(n_occ, occupancy)) %>% distinct()) %>%
  relocate(study, database, studyID, study_name, year_period, period, specieskey, species)


## calculate occupancy change
occupancy_period_occup <- occupancy_period %>%
  dplyr::select(-c(year_period, year_mean_period, n_occ)) %>%
  group_by(study, specieskey) %>% 
  pivot_wider(names_from = period, names_prefix = "occup_", values_from = occupancy) %>%
  ungroup()

occupancy_period_nocc <- occupancy_period %>%
  dplyr::select(-c(year_period, year_mean_period, occupancy)) %>%
  group_by(study, specieskey) %>% 
  pivot_wider(names_from = period, names_prefix = "nocc_", values_from = n_occ) %>%
  ungroup()

occupancy_change_period <- occupancy_period_occup %>% 
  left_join(occupancy_period_nocc) %>% 
  mutate(dynamic = ifelse(occup_early > 0 & occup_late > 0, "persistent", 
                          ifelse(occup_early > 0 & occup_late == 0, "extinction", "colonization")),
         occup_change = occup_late - occup_early) %>%
  relocate(dynamic, occup_early, occup_late, occup_change, nocc_early, nocc_late, .after = last_col()) %>%
  filter(occup_early > 0 | occup_late > 0) 

table(occupancy_change_period[, c("dynamic", "database")])


## calculate abundance change
abundance_change_period <- abundance_period %>%
  dplyr::select(-c(year_period, year_mean_period)) %>%
  group_by(study, specieskey) %>% 
  pivot_wider(names_from = period, names_prefix = "N_", values_from = abundance) %>%
  ungroup() %>%
  mutate(N_change = N_late - N_early,
         dynamic = ifelse(N_early > 0 & N_late > 0, "persistent", 
                          ifelse(N_early > 0 & N_late == 0, "extinction", "colonization"))) %>%
  filter(N_early > 0 | N_late > 0) # 20 species in bt_176 have 0-0 records


# save occupancy and abundance change data
dir.create("intermediate_results", showWarnings = FALSE)
save(occupancy_2yr, occupancy_period, occupancy_change_2yr, occupancy_change_period, 
     abundance_2yr, abundance_period, abundance_change_2yr, abundance_change_period, 
     years_period, dat_meta, 
     file = "intermediate_results/Assemblage_occupancy_abundance_change.RDATA")
