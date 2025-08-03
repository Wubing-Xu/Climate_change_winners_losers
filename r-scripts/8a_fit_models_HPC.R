######### organize data and then fit models in HPC

rm(list = ls())

# Set user dependent working directories
path2wd <- "/work/wubing/Populations_warming/analyses_202406"
setwd(path2wd)

# load packages
packages <- c("tidyverse", "brms")

for(x in packages){
  if(!require(x, character.only=TRUE)){
    if(!require(x,character.only=TRUE, lib.loc = "/gpfs1/schlecker/home/wubing/R/x86_64-pc-linux-gnu-library/4.2")){
      install.packages(x, repos = "http://cran.us.r-project.org", lib.loc = "/gpfs1/schlecker/home/wubing/R/x86_64-pc-linux-gnu-library/4.2", dependencies = TRUE)
      require(x, character.only = TRUE, lib.loc = "/gpfs1/schlecker/home/wubing/R/x86_64-pc-linux-gnu-library/4.2")					
    }
  }
}

# task index
index <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))

# load occupancy change and meta data
load("intermediate_results/Assemblage_occupancy_abundance_change.RDATA")
load("intermediate_results/Assemblage_thermal_position_tempChange.RDATA")
load("data/Assemblages_regions.RDATA")
load("data/Assemblages_taxa.RDATA")


# combine meta data 
dat_meta <- dat_meta %>% 
  left_join(region %>% dplyr::select(study, region)) %>% 
  left_join(temp_early_change_period) %>%
  left_join(temp_first_change_2yr %>%  setNames(c("study", "realm", "tempmean_2yr", "tempmax_2yr", "tempmin_2yr", 
                                                  "tempmean_sampyr_2yr",  "change_tempmean_2yr", "change_tempmax_2yr", 
                                                  "change_tempmin_2yr", "change_tempmean_sampyr_2yr"))) %>%
  dplyr::select(!taxon_new) %>%
  left_join(dat_meta_taxa %>% dplyr::select(study, taxon_new, taxon_final)) %>%
  relocate(taxon_new, taxon_final, .after = taxon)


# data for occupancy change between periods
oc_period <- occupancy_change_period %>%
  # calculate occupancy change as the change in the logit of occupancy between two periods (e.g. log odds ratio of occupancy; standardized by duration)
  mutate(occup_early_logit = ifelse(n_samp<100, log((occup_early + 0.005)/(1- occup_early + 0.005)), 
                                    log((occup_early + 0.5/n_samp)/(1- occup_early + 0.5/n_samp))),
         occup_late_logit = ifelse(n_samp<100, log((occup_late + 0.005)/(1- occup_late + 0.005)), 
                                   log((occup_late + 0.5/n_samp)/(1- occup_late + 0.5/n_samp))),
         occup_change_logit = (occup_late_logit - occup_early_logit)/(duration_mean - 1)) %>%
  dplyr::select(-c(occup_early_logit, occup_late_logit)) %>% 
  left_join(dat_meta %>% dplyr::select(study, taxon_new, taxon_final, realm, region, climate, extent_km2, subset_abundance)) %>% 
  left_join(temp_early_change_period %>% dplyr::select(study, change_tempmean:change_tempmean_sampyr)) %>%
  # standardized temperature change by duration
  mutate(change_tempmean = change_tempmean/(duration_mean -1),
         change_tempmax = change_tempmax/(duration_mean -1),
         change_tempmin = change_tempmin/(duration_mean -1),
         change_tempmean_sampyr = change_tempmean_sampyr/(duration_mean -1)) %>% 
  inner_join(tempos_early_change_period %>% dplyr::select(study:tmeanpos_sampyr)) %>%
  inner_join(spe_thermal_range %>% dplyr::select(c(species:tempmean_median, tempmean_range))) 

# data for occupancy change between the last and first years
oc_2yr <- occupancy_change_2yr %>% 
  # calculate occupancy change as the change in the logit of occupancy between two time points (e.g. log odds ratio of occupancy; standardized by duration)
  mutate(occup_early_logit = ifelse(n_samp<100, log((occup_early + 0.005)/(1- occup_early + 0.005)), 
                                    log((occup_early + 0.5/n_samp)/(1- occup_early + 0.5/n_samp))),
         occup_late_logit = ifelse(n_samp<100, log((occup_late + 0.005)/(1- occup_late + 0.005)), 
                                   log((occup_late + 0.5/n_samp)/(1- occup_late + 0.5/n_samp))),
         occup_change_logit = (occup_late_logit - occup_early_logit)/(duration - 1)) %>%
  dplyr::select(-c(occup_early_logit, occup_late_logit)) %>% 
  left_join(dat_meta %>% dplyr::select(study, taxon_new, taxon_final, realm, region, climate, extent_km2, subset_abundance)) %>% 
  left_join(temp_first_change_2yr %>% dplyr::select(study, change_tempmean:change_tempmean_sampyr)) %>%
  # standardized temperature change by duration
  mutate(change_tempmean = change_tempmean/(duration -1),
         change_tempmax = change_tempmax/(duration -1),
         change_tempmin = change_tempmin/(duration -1),
         change_tempmean_sampyr = change_tempmean_sampyr/(duration -1)) %>% 
  inner_join(tempos_first_change_2yr %>% dplyr::select(study:tmeanpos_sampyr)) %>%
  inner_join(spe_thermal_range %>% dplyr::select(c(species:tempmean_median, tempmean_range)))


# data for abundance change between late and early periods
ac_period <- abundance_change_period %>%
  # calculate abundance change as the change in the logarithm  of abundance between two periods (e.g. log ratio of abundance; standardized by duration)
  group_by(study) %>%
  mutate(is_zero = min(c(min(N_early), min(N_late))) == 0,
         N_min = min(c(1, min(N_early[N_early>0]), min(N_late[N_late>0]))),
         N_early_no0 = ifelse(is_zero, N_early + N_min, N_early),
         N_late_no0 = ifelse(is_zero, N_late + N_min, N_late)) %>% 
  ungroup() %>%
  mutate(N_change_log = (log(N_late_no0) - log(N_early_no0))/(duration_mean - 1)) %>%
  dplyr::select(-c(is_zero, N_min, N_early_no0, N_late_no0)) %>%
  left_join(dat_meta %>% dplyr::select(study, taxon_new, taxon_final, realm, region, climate, extent_km2)) %>% 
  left_join(temp_early_change_period %>% dplyr::select(study, change_tempmean:change_tempmean_sampyr)) %>%
  # standardized temperature change by duration
  mutate(change_tempmean = change_tempmean/(duration_mean -1),
         change_tempmax = change_tempmax/(duration_mean -1),
         change_tempmin = change_tempmin/(duration_mean -1),
         change_tempmean_sampyr = change_tempmean_sampyr/(duration_mean -1)) %>% 
  inner_join(tempos_early_change_period %>% dplyr::select(study:tmeanpos_sampyr)) %>%
  inner_join(spe_thermal_range %>% dplyr::select(c(species:tempmean_median, tempmean_range))) 

# data for abundance change between last and first years
ac_2yr <- abundance_change_2yr %>%
  # calculate abundance change as the change in the logarithm  of abundance between two time points (e.g. log ratio of abundance; standardized by duration)
  group_by(study) %>%
  mutate(is_zero = min(c(min(N_early), min(N_late))) == 0,
         N_min = min(c(1, min(N_early[N_early>0]), min(N_late[N_late>0]))),
         N_early_no0 = ifelse(is_zero, N_early + N_min, N_early),
         N_late_no0 = ifelse(is_zero, N_late + N_min, N_late)) %>% 
  ungroup() %>%
  mutate(N_change_log = (log(N_late_no0) - log(N_early_no0))/(duration - 1)) %>%
  dplyr::select(-c(is_zero, N_min, N_early_no0, N_late_no0)) %>%
  left_join(dat_meta %>% dplyr::select(study, taxon_new, taxon_final, realm, region, climate, extent_km2)) %>% 
  left_join(temp_first_change_2yr %>% dplyr::select(study, change_tempmean:change_tempmean_sampyr)) %>%
  # standardized temperature change by duration
  mutate(change_tempmean = change_tempmean/(duration -1),
         change_tempmax = change_tempmax/(duration -1),
         change_tempmin = change_tempmin/(duration -1),
         change_tempmean_sampyr = change_tempmean_sampyr/(duration -1)) %>% 
  inner_join(tempos_first_change_2yr %>% dplyr::select(study:tmeanpos_sampyr)) %>%
  inner_join(spe_thermal_range %>% dplyr::select(c(species:tempmean_median, tempmean_range)))


# set factor levels of categorical groups
oc_period <- oc_period %>%
  mutate(realm = factor(realm, levels = c("Terrestrial", "Freshwater", "Marine")),
         taxon_new = factor(taxon_new, levels = c("Amphibians and reptiles", "Birds", "Fish", 
                                                  "Invertebrates", "Mammals", "Plants")),
         taxon_final = factor(taxon_final, levels = c("Amphibians and reptiles", "Birds", 
                                                      "Terrestrial invertebrates", "Terrestrial mammals", "Terrestrial plants", 
                                                      "Freshwater fish", "Freshwater invertebrates", "Freshwater plants", 
                                                      "Marine fish", "Marine invertebrates", "Marine mammals", "Marine plants")),
         region = factor(region, levels = c("Africa", "Asia", "Australia", "Europe", "North America", "South America",
                                            "Australian Ocean", "East Atlantic", "West Atlantic", "Pacific Ocean"))) %>% 
  mutate(realm_taxa = factor(paste(realm, taxon_new, sep = "_")),
         realm_region = factor(paste(realm, region, sep = "_")),
         realm_climate = factor(paste(realm, climate, sep = "_")))

oc_2yr <- oc_2yr %>%
  mutate(realm = factor(realm, levels = c("Terrestrial", "Freshwater", "Marine")),
         taxon_new = factor(taxon_new, levels = c("Amphibians and reptiles", "Birds", "Fish", 
                                                  "Invertebrates", "Mammals", "Plants")),
         taxon_final = factor(taxon_final, levels = c("Amphibians and reptiles", "Birds", 
                                                      "Terrestrial invertebrates", "Terrestrial mammals", "Terrestrial plants", 
                                                      "Freshwater fish", "Freshwater invertebrates", "Freshwater plants", 
                                                      "Marine fish", "Marine invertebrates", "Marine mammals", "Marine plants")),
         region = factor(region, levels = c("Africa", "Asia", "Australia", "Europe", "North America", "South America",
                                            "Australian Ocean", "East Atlantic", "West Atlantic", "Pacific Ocean"))) %>% 
  mutate(realm_taxa = factor(paste(realm, taxon_new, sep = "_")),
         realm_region = factor(paste(realm, region, sep = "_")),
         realm_climate = factor(paste(realm, climate, sep = "_")))

ac_period <- ac_period %>%
  mutate(realm = factor(realm, levels = c("Terrestrial", "Freshwater", "Marine")),
         taxon_new = factor(taxon_new, levels = c("Amphibians and reptiles", "Birds", "Fish", 
                                                  "Invertebrates", "Mammals", "Plants")),
         taxon_final = factor(taxon_final, levels = c("Amphibians and reptiles", "Birds", 
                                                      "Terrestrial invertebrates", "Terrestrial mammals", "Terrestrial plants", 
                                                      "Freshwater fish", "Freshwater invertebrates", "Freshwater plants", 
                                                      "Marine fish", "Marine invertebrates", "Marine mammals", "Marine plants")),
         region = factor(region, levels = c("Africa", "Asia", "Australia", "Europe", "North America", "South America",
                                            "Australian Ocean", "East Atlantic", "West Atlantic", "Pacific Ocean"))) %>% 
  mutate(realm_taxa = factor(paste(realm, taxon_new, sep = "_")),
         realm_region = factor(paste(realm, region, sep = "_")),
         realm_climate = factor(paste(realm, climate, sep = "_")))

ac_2yr <- ac_2yr %>%
  mutate(realm = factor(realm, levels = c("Terrestrial", "Freshwater", "Marine")),
         taxon_new = factor(taxon_new, levels = c("Amphibians and reptiles", "Birds", "Fish", 
                                                  "Invertebrates", "Mammals", "Plants")),
         taxon_final = factor(taxon_final, levels = c("Amphibians and reptiles", "Birds", 
                                                      "Terrestrial invertebrates", "Terrestrial mammals", "Terrestrial plants", 
                                                      "Freshwater fish", "Freshwater invertebrates", "Freshwater plants", 
                                                      "Marine fish", "Marine invertebrates", "Marine mammals", "Marine plants")),
         region = factor(region, levels = c("Africa", "Asia", "Australia", "Europe", "North America", "South America",
                                            "Australian Ocean", "East Atlantic", "West Atlantic", "Pacific Ocean"))) %>% 
  mutate(realm_taxa = factor(paste(realm, taxon_new, sep = "_")),
         realm_region = factor(paste(realm, region, sep = "_")),
         realm_climate = factor(paste(realm, climate, sep = "_")))


# save organized data for fitting models
dir.create("results", showWarnings = FALSE)
# save(oc_period, oc_2yr, ac_period, ac_2yr, dat_meta, file = "results/data_input_to_models.RDATA")
# load("results/data_input_to_models.RDATA")


########################
# fit linear mixed model
dir.create("models", showWarnings = FALSE)


## the  overall effect of thermal position on occupancy and abundance change
if(index == 1){
  t1 <- Sys.time()
  brm_oc_tmeanpos <- brm(bf(occup_change_logit ~ tmeanpos + (1 + tmeanpos|study), 
                            sigma ~ log(n_samp) + log(duration) + log(extent_km2)),
                         family = gaussian(), data = oc_period ,
                         chains = 4, cores = 4, iter = 4000, thin = 2,
                         control = list(adapt_delta = 0.9, max_treedepth = 10), 
                         file = "models/brm_oc_tmeanpos")
  print(brm_oc_tmeanpos)
  t2 <- Sys.time()
  print(t2 - t1)
  
  brm_ac_tmeanpos <- brm(bf(N_change_log ~ tmeanpos + (1 + tmeanpos|study), 
                            sigma ~ log(n_samp) + log(duration) + log(extent_km2)),
                         family = gaussian(), data = ac_period ,
                         chains = 4, cores = 4, iter = 4000, thin = 2,
                         control = list(adapt_delta = 0.9, max_treedepth = 10), 
                         file = "models/brm_ac_tmeanpos")
  print(brm_ac_tmeanpos)
  t3 <- Sys.time()
  print(t3 - t2)
}


## the effect of thermal position on occupancy and abundance change across realms
if(index == 2){
  t1 <- Sys.time()
  brm_oc_tmeanpos_realm <- brm(bf(occup_change_logit ~ 0 + realm + tmeanpos:realm + (1 + tmeanpos|study), 
                            sigma ~ log(n_samp) + log(duration) + log(extent_km2)),
                         family = gaussian(), data = oc_period ,
                         chains = 4, cores = 4, iter = 4000, thin = 2,
                         control = list(adapt_delta = 0.9, max_treedepth = 10), 
                         file = "models/brm_oc_tmeanpos_realm")
  print(brm_oc_tmeanpos_realm)
  t2 <- Sys.time()
  print(t2 - t1)
  
  brm_ac_tmeanpos_realm <- brm(bf(N_change_log ~ 0 + realm + tmeanpos:realm + (1 + tmeanpos|study), 
                                  sigma ~ log(n_samp) + log(duration) + log(extent_km2)),
                         family = gaussian(), data = ac_period ,
                         chains = 4, cores = 4, iter = 4000, thin = 2,
                         control = list(adapt_delta = 0.9, max_treedepth = 10), 
                         file = "models/brm_ac_tmeanpos_realm")
  print(brm_ac_tmeanpos_realm)
  t3 <- Sys.time()
  print(t3 - t2)
}


## the effect of thermal position across taxa: use the group with > 5 studies 
if(index == 3){
  nstudy_realm_taxa <- oc_period %>% distinct(study, realm_taxa) %>% group_by(realm_taxa) %>% count() %>% filter(n >5)
  t1 <- Sys.time()
  brm_oc_tmeanpos_taxa <- brm(bf(occup_change_logit ~ 0 + realm_taxa + tmeanpos:realm_taxa + (1 + tmeanpos|study), 
                                  sigma ~ log(n_samp) + log(duration) + log(extent_km2)),
                               family = gaussian(), data = oc_period %>% filter(realm_taxa %in% nstudy_realm_taxa$realm_taxa),
                               chains = 4, cores = 4, iter = 4000, thin = 2,
                               control = list(adapt_delta = 0.9, max_treedepth = 10), 
                               file = "models/brm_oc_tmeanpos_taxa")
  print(brm_oc_tmeanpos_taxa)
  t2 <- Sys.time()
  print(t2 - t1)
  
  nstudy_realm_taxa <- ac_period %>% distinct(study, realm_taxa) %>% group_by(realm_taxa) %>% count() %>% filter(n >5)
  brm_ac_tmeanpos_taxa <- brm(bf(N_change_log ~ 0 + realm_taxa + tmeanpos:realm_taxa + (1 + tmeanpos|study), 
                                  sigma ~ log(n_samp) + log(duration) + log(extent_km2)),
                               family = gaussian(), data = ac_period %>% filter(realm_taxa %in% nstudy_realm_taxa$realm_taxa),
                               chains = 4, cores = 4, iter = 4000, thin = 2,
                               control = list(adapt_delta = 0.9, max_treedepth = 10), 
                               file = "models/brm_ac_tmeanpos_taxa")
  print(brm_ac_tmeanpos_taxa)
  t3 <- Sys.time()
  print(t3 - t2)
}


## the effect of thermal position across regions: use the group with > 5 studies 
if(index == 4){
  nstudy_realm_region <- oc_period %>% distinct(study, realm_region) %>% group_by(realm_region) %>% count() %>% filter(n >5)
  t1 <- Sys.time()
  brm_oc_tmeanpos_region <- brm(bf(occup_change_logit ~ 0 + realm_region + tmeanpos:realm_region + (1 + tmeanpos|study), 
                                 sigma ~ log(n_samp) + log(duration) + log(extent_km2)),
                              family = gaussian(), data = oc_period %>% filter(realm_region %in% nstudy_realm_region$realm_region),
                              chains = 4, cores = 4, iter = 4000, thin = 2,
                              control = list(adapt_delta = 0.9, max_treedepth = 10), 
                              file = "models/brm_oc_tmeanpos_region")
  print(brm_oc_tmeanpos_region)
  t2 <- Sys.time()
  print(t2 - t1)
  
  nstudy_realm_region <- ac_period %>% distinct(study, realm_region) %>% group_by(realm_region) %>% count() %>% filter(n >5)
  brm_ac_tmeanpos_region <- brm(bf(N_change_log ~ 0 + realm_region + tmeanpos:realm_region + (1 + tmeanpos|study), 
                                 sigma ~ log(n_samp) + log(duration) + log(extent_km2)),
                              family = gaussian(), data = ac_period %>% filter(realm_region %in% nstudy_realm_region$realm_region),
                              chains = 4, cores = 4, iter = 4000, thin = 2,
                              control = list(adapt_delta = 0.9, max_treedepth = 10), 
                              file = "models/brm_ac_tmeanpos_region")
  print(brm_ac_tmeanpos_region)
  t3 <- Sys.time()
  print(t3 - t2)
}


## the overall effect of thermal position and interactions with temperature change
if(index == 5){
  t1 <- Sys.time()
  brm_oc_tmeanpos_tmeanchange <- brm(bf(occup_change_logit ~ tmeanpos*change_tempmean + (1 + tmeanpos|study), 
                                              sigma ~ log(n_samp) + log(duration) + log(extent_km2)),
                                           family = gaussian(), data = oc_period ,
                                           chains = 4, cores = 4, iter = 4000, thin = 2,
                                           control = list(adapt_delta = 0.9, max_treedepth = 10), 
                                           file = "models/brm_oc_tmeanpos_tmeanchange")
  print(brm_oc_tmeanpos_tmeanchange)
  t2 <- Sys.time()
  print(t2 - t1)
  
  brm_ac_tmeanpos_tmeanchange <- brm(bf(N_change_log ~ tmeanpos*change_tempmean + (1 + tmeanpos|study), 
                                              sigma ~ log(n_samp) + log(duration) + log(extent_km2)),
                                           family = gaussian(), data = ac_period ,
                                           chains = 4, cores = 4, iter = 4000, thin = 2,
                                           control = list(adapt_delta = 0.9, max_treedepth = 10), 
                                           file = "models/brm_ac_tmeanpos_tmeanchange")
  print(brm_ac_tmeanpos_tmeanchange)
  t3 <- Sys.time()
  print(t3 - t2)
}


## the effect of thermal position and interactions with temperature change across realms
if(index == 6){
  t1 <- Sys.time()
  brm_oc_tmeanpos_tmeanchange_realm <- brm(bf(occup_change_logit ~ 0 + realm +  (tmeanpos*change_tempmean):realm + (1 + tmeanpos|study), 
                                                    sigma ~ log(n_samp) + log(duration) + log(extent_km2)),
                                                 family = gaussian(), data = oc_period ,
                                                 chains = 4, cores = 4, iter = 4000, thin = 2,
                                                 control = list(adapt_delta = 0.9, max_treedepth = 10), 
                                                 file = "models/brm_oc_tmeanpos_tmeanchange_realm")
  print(brm_oc_tmeanpos_tmeanchange_realm)
  t2 <- Sys.time()
  print(t2 - t1)

  brm_ac_tmeanpos_tmeanchange_realm <- brm(bf(N_change_log ~ 0 + realm +  (tmeanpos*change_tempmean):realm + (1 + tmeanpos|study), 
                                              sigma ~ log(n_samp) + log(duration) + log(extent_km2)),
                                           family = gaussian(), data = ac_period ,
                                           chains = 4, cores = 4, iter = 4000, thin = 2,
                                           control = list(adapt_delta = 0.9, max_treedepth = 10), 
                                           file = "models/brm_ac_tmeanpos_tmeanchange_realm")
  print(brm_ac_tmeanpos_tmeanchange_realm)
  t3 <- Sys.time()
  print(t3 - t2)
}


## differences of tested effects across taxa: use the group with > 5 studies 
if(index == 7){
  nstudy_realm_taxa <- oc_period %>% distinct(study, realm_taxa) %>% group_by(realm_taxa) %>% count() %>% filter(n >5)
  t1 <- Sys.time()
  brm_oc_tmeanpos_tmeanchange_taxa <- brm(bf(occup_change_logit ~ 0 + realm_taxa +  (tmeanpos*change_tempmean):realm_taxa + (1 + tmeanpos|study), 
                                             sigma ~ log(n_samp) + log(duration) + log(extent_km2)),
                                          family = gaussian(), data = oc_period %>% filter(realm_taxa %in% nstudy_realm_taxa$realm_taxa),
                                          chains = 4, cores = 4, iter = 4000, thin = 2,
                                          control = list(adapt_delta = 0.9, max_treedepth = 10), 
                                          file = "models/brm_oc_tmeanpos_tmeanchange_taxa")
  print(brm_oc_tmeanpos_tmeanchange_taxa)
  t2 <- Sys.time()
  print(t2 - t1)
  
  nstudy_realm_taxa <- ac_period %>% distinct(study, realm_taxa) %>% group_by(realm_taxa) %>% count() %>% filter(n >5)
  brm_ac_tmeanpos_tmeanchange_taxa <- brm(bf(N_change_log ~ 0 + realm_taxa +  (tmeanpos*change_tempmean):realm_taxa + (1 + tmeanpos|study), 
                                             sigma ~ log(n_samp) + log(duration) + log(extent_km2)),
                                          family = gaussian(), data = ac_period %>% filter(realm_taxa %in% nstudy_realm_taxa$realm_taxa),
                                          chains = 4, cores = 4, iter = 4000, thin = 2,
                                          control = list(adapt_delta = 0.9, max_treedepth = 10), 
                                          file = "models/brm_ac_tmeanpos_tmeanchange_taxa")
  print(brm_ac_tmeanpos_tmeanchange_taxa)
  t3 <- Sys.time()
  print(t3 - t2)
}


## differences of tested effects across regions: use the group with > 5 studies 
if(index == 8){
  nstudy_realm_region <- oc_period %>% distinct(study, realm_region) %>% group_by(realm_region) %>% count() %>% filter(n >5)
  
  t1 <- Sys.time()
  brm_oc_tmeanpos_tmeanchange_region <- brm(bf(occup_change_logit ~ 0 + realm_region +  (tmeanpos*change_tempmean):realm_region + (1 + tmeanpos|study), 
                                               sigma ~ log(n_samp) + log(duration) + log(extent_km2)),
                                            family = gaussian(), data = oc_period %>% filter(realm_region %in% nstudy_realm_region$realm_region),
                                            chains = 4, cores = 4, iter = 4000, thin = 2,
                                            control = list(adapt_delta = 0.9, max_treedepth = 10), 
                                            file = "models/brm_oc_tmeanpos_tmeanchange_region")
  print(brm_oc_tmeanpos_tmeanchange_region)
  t2 <- Sys.time()
  print(t2 - t1)
  
  
  nstudy_realm_region <- ac_period %>% distinct(study, realm_region) %>% group_by(realm_region) %>% count() %>% filter(n >5)
  
  brm_ac_tmeanpos_tmeanchange_region <- brm(bf(N_change_log ~ 0 + realm_region +  (tmeanpos*change_tempmean):realm_region + (1 + tmeanpos|study), 
                                               sigma ~ log(n_samp) + log(duration) + log(extent_km2)),
                                            family = gaussian(), data = ac_period %>% filter(realm_region %in% nstudy_realm_region$realm_region),
                                            chains = 4, cores = 4, iter = 4000, thin = 2,
                                            control = list(adapt_delta = 0.9, max_treedepth = 10), 
                                            file = "models/brm_ac_tmeanpos_tmeanchange_region")
  print(brm_ac_tmeanpos_tmeanchange_region)
  t3 <- Sys.time()
  print(t3 - t2)
}


## sensitivity analyses: exclude species with thermal range < 5 degree
if(index == 9){
  # test the overall effects of baseline thermal position
  t1 <- Sys.time()
  brm_oc_tmeanpos_widesp <- brm(bf(occup_change_logit ~ 0 + realm + tmeanpos:realm + (1 + tmeanpos|study), 
                                  sigma ~ log(n_samp) + log(duration) + log(extent_km2)),
                               family = gaussian(), data = oc_period %>% filter(tempmean_range >= 5 & nocc_climate>=100),
                               chains = 4, cores = 4, iter = 4000, thin = 4,
                               control = list(adapt_delta = 0.9, max_treedepth = 10), 
                               file = "models/brm_oc_tmeanpos_widesp")
  print(brm_oc_tmeanpos_widesp)
  t2 <- Sys.time()
  print(t2 - t1)
  
  brm_ac_tmeanpos_widesp <- brm(bf(N_change_log ~ 0 + realm + tmeanpos:realm + (1 + tmeanpos|study), 
                                  sigma ~ log(n_samp) + log(duration) + log(extent_km2)),
                               family = gaussian(), data = ac_period %>% filter(tempmean_range >= 5 & nocc_climate>=100),
                               chains = 4, cores = 4, iter = 4000, thin = 4,
                               control = list(adapt_delta = 0.9, max_treedepth = 10), 
                               file = "models/brm_ac_tmeanpos_widesp")
  print(brm_ac_tmeanpos_widesp)
  t3 <- Sys.time()
  print(t3 - t2)
  
  # test the interaction between baseline thermal position and temperature change
  brm_oc_tmeanpos_tmeanchange_widesp <- brm(bf(occup_change_logit ~ 0 + realm +  (tmeanpos*change_tempmean):realm + (1 + tmeanpos|study), 
                                              sigma ~ log(n_samp) + log(duration) + log(extent_km2)),
                                           family = gaussian(), data = oc_period %>% filter(tempmean_range >= 5 & nocc_climate>=100),
                                           chains = 4, cores = 4, iter = 4000, thin = 4,
                                           control = list(adapt_delta = 0.9, max_treedepth = 10), 
                                           file = "models/brm_oc_tmeanpos_tmeanchange_widesp")
  print(brm_oc_tmeanpos_tmeanchange_widesp)
  t4 <- Sys.time()
  print(t4 - t3)
  
  brm_ac_tmeanpos_tmeanchange_widesp <- brm(bf(N_change_log ~ 0 + realm +  (tmeanpos*change_tempmean):realm + (1 + tmeanpos|study), 
                                               sigma ~ log(n_samp) + log(duration) + log(extent_km2)),
                                           family = gaussian(), data = ac_period %>% filter(tempmean_range >= 5 & nocc_climate>=100),
                                           chains = 4, cores = 4, iter = 4000, thin = 4,
                                           control = list(adapt_delta = 0.9, max_treedepth = 10), 
                                           file = "models/brm_ac_tmeanpos_tmeanchange_widesp")
  print(brm_ac_tmeanpos_tmeanchange_widesp)
  t5 <- Sys.time()
  print(t5 - t4)
}


##  sensitivity analyses: thermal position and temperature change based on max monthly temperature
if(index == 10){
  # test the overall effects of baseline thermal position
  t1 <- Sys.time()
  brm_oc_tmaxpos_realm <- brm(bf(occup_change_logit ~ 0 + realm + tmaxpos:realm + (1 + tmaxpos|study), 
                                   sigma ~ log(n_samp) + log(duration) + log(extent_km2)),
                                family = gaussian(), data = oc_period,
                                chains = 4, cores = 4, iter = 4000, thin = 4,
                                control = list(adapt_delta = 0.9, max_treedepth = 10), 
                                file = "models/brm_oc_tmaxpos_realm")
  print(brm_oc_tmaxpos_realm)
  t2 <- Sys.time()
  print(t2 - t1)
  
  brm_ac_tmaxpos_realm <- brm(bf(N_change_log ~ 0 + realm + tmaxpos:realm + (1 + tmaxpos|study), 
                                   sigma ~ log(n_samp) + log(duration) + log(extent_km2)),
                                family = gaussian(), data = ac_period,
                                chains = 4, cores = 4, iter = 4000, thin = 4,
                                control = list(adapt_delta = 0.9, max_treedepth = 10), 
                                file = "models/brm_ac_tmaxpos_realm")
  print(brm_ac_tmaxpos_realm)
  t3 <- Sys.time()
  print(t3 - t2)
  
  # test the interaction between baseline thermal position and temperature change
  brm_oc_tmaxpos_tmaxchange_realm <- brm(bf(occup_change_logit ~ 0 + realm +  (tmaxpos*change_tempmax):realm + (1 + tmaxpos|study), 
                                               sigma ~ log(n_samp) + log(duration) + log(extent_km2)),
                                            family = gaussian(), data = oc_period,
                                            chains = 4, cores = 4, iter = 4000, thin = 4,
                                            control = list(adapt_delta = 0.9, max_treedepth = 10), 
                                            file = "models/brm_oc_tmaxpos_tmaxchange_realm")
  print(brm_oc_tmaxpos_tmaxchange_realm)
  t4 <- Sys.time()
  print(t4 - t3)
  
  brm_ac_tmaxpos_tmaxchange_realm <- brm(bf(N_change_log ~ 0 + realm +  (tmaxpos*change_tempmax):realm + (1 + tmaxpos|study), 
                                               sigma ~ log(n_samp) + log(duration) + log(extent_km2)),
                                            family = gaussian(), data = ac_period,
                                            chains = 4, cores = 4, iter = 4000, thin = 4,
                                            control = list(adapt_delta = 0.9, max_treedepth = 10), 
                                            file = "models/brm_ac_tmaxpos_tmaxchange_realm")
  print(brm_ac_tmaxpos_tmaxchange_realm)
  t5 <- Sys.time()
  print(t5 - t4)
}


##  sensitivity analyses: thermal position and temperature change based on min monthly temperature
if(index == 11){
  # test the overall effects of baseline thermal position
  t1 <- Sys.time()
  brm_oc_tminpos_realm <- brm(bf(occup_change_logit ~ 0 + realm + tminpos:realm + (1 + tminpos|study), 
                                 sigma ~ log(n_samp) + log(duration) + log(extent_km2)),
                              family = gaussian(), data = oc_period,
                              chains = 4, cores = 4, iter = 4000, thin = 4,
                              control = list(adapt_delta = 0.9, max_treedepth = 10), 
                              file = "models/brm_oc_tminpos_realm")
  print(brm_oc_tminpos_realm)
  t2 <- Sys.time()
  print(t2 - t1)
  
  brm_ac_tminpos_realm <- brm(bf(N_change_log ~ 0 + realm + tminpos:realm + (1 + tminpos|study), 
                                 sigma ~ log(n_samp) + log(duration) + log(extent_km2)),
                              family = gaussian(), data = ac_period,
                              chains = 4, cores = 4, iter = 4000, thin = 4,
                              control = list(adapt_delta = 0.9, max_treedepth = 10), 
                              file = "models/brm_ac_tminpos_realm")
  print(brm_ac_tminpos_realm)
  t3 <- Sys.time()
  print(t3 - t2)
  
  # test the interaction between baseline thermal position and temperature change
  brm_oc_tminpos_tminchange_realm <- brm(bf(occup_change_logit ~ 0 + realm +  (tminpos*change_tempmin):realm + (1 + tminpos|study), 
                                            sigma ~ log(n_samp) + log(duration) + log(extent_km2)),
                                         family = gaussian(), data = oc_period,
                                         chains = 4, cores = 4, iter = 4000, thin = 4,
                                         control = list(adapt_delta = 0.9, max_treedepth = 10), 
                                         file = "models/brm_oc_tminpos_tminchange_realm")
  print(brm_oc_tminpos_tminchange_realm)
  t4 <- Sys.time()
  print(t4 - t3)
  
  brm_ac_tminpos_tminchange_realm <- brm(bf(N_change_log ~ 0 + realm +  (tminpos*change_tempmin):realm + (1 + tminpos|study), 
                                            sigma ~ log(n_samp) + log(duration) + log(extent_km2)),
                                         family = gaussian(), data = ac_period,
                                         chains = 4, cores = 4, iter = 4000, thin = 4,
                                         control = list(adapt_delta = 0.9, max_treedepth = 10), 
                                         file = "models/brm_ac_tminpos_tminchange_realm")
  print(brm_ac_tminpos_tminchange_realm)
  t5 <- Sys.time()
  print(t5 - t4)
}


## sensitivity analyses: use the meta-communities with small extent
if(index == 12){
  # test the overall effects of baseline thermal position
  t1 <- Sys.time()
  brm_oc_tmeanpos_smregion <- brm(bf(occup_change_logit ~ 0 + realm + tmeanpos:realm + (1 + tmeanpos|study), 
                                   sigma ~ log(n_samp) + log(duration) + log(extent_km2)),
                                family = gaussian(), data = oc_period %>% filter(extent_km2 < 10000),
                                chains = 4, cores = 4, iter = 4000, thin = 4,
                                control = list(adapt_delta = 0.9, max_treedepth = 10), 
                                file = "models/brm_oc_tmeanpos_smregion")
  print(brm_oc_tmeanpos_smregion)
  t2 <- Sys.time()
  print(t2 - t1)
  
  brm_ac_tmeanpos_smregion <- brm(bf(N_change_log ~ 0 + realm + tmeanpos:realm + (1 + tmeanpos|study), 
                                   sigma ~ log(n_samp) + log(duration) + log(extent_km2)),
                                family = gaussian(), data = ac_period %>% filter(extent_km2 < 10000),
                                chains = 4, cores = 4, iter = 4000, thin = 4,
                                control = list(adapt_delta = 0.9, max_treedepth = 10), 
                                file = "models/brm_ac_tmeanpos_smregion")
  print(brm_ac_tmeanpos_smregion)
  t3 <- Sys.time()
  print(t3 - t2)
  
  # test the interaction between baseline thermal position and temperature change
  brm_oc_tmeanpos_tmeanchange_smregion <- brm(bf(occup_change_logit ~ 0 + realm +  (tmeanpos*change_tempmean):realm + (1 + tmeanpos|study), 
                                               sigma ~ log(n_samp) + log(duration) + log(extent_km2)),
                                            family = gaussian(), data = oc_period %>% filter(extent_km2 < 10000),
                                            chains = 4, cores = 4, iter = 4000, thin = 4,
                                            control = list(adapt_delta = 0.9, max_treedepth = 10), 
                                            file = "models/brm_oc_tmeanpos_tmeanchange_smregion")
  print(brm_oc_tmeanpos_tmeanchange_smregion)
  t4 <- Sys.time()
  print(t4 - t3)
  
  brm_ac_tmeanpos_tmeanchange_smregion <- brm(bf(N_change_log ~ 0 + realm +  (tmeanpos*change_tempmean):realm + (1 + tmeanpos|study), 
                                               sigma ~ log(n_samp) + log(duration) + log(extent_km2)),
                                            family = gaussian(), data = ac_period %>% filter(extent_km2 < 10000),
                                            chains = 4, cores = 4, iter = 4000, thin = 4,
                                            control = list(adapt_delta = 0.9, max_treedepth = 10), 
                                            file = "models/brm_ac_tmeanpos_tmeanchange_smregion")
  print(brm_ac_tmeanpos_tmeanchange_smregion)
  t5 <- Sys.time()
  print(t5 - t4)
}


# sensitivity analyses: use the data in the first and last sampling years
if(index == 13){
  # test the overall effects of baseline thermal position
  t1 <- Sys.time()
  brm_oc_tmeanpos_2yr <- brm(bf(occup_change_logit ~ 0 + realm + tmeanpos:realm + (1 + tmeanpos|study), 
                                     sigma ~ log(n_samp) + log(duration) + log(extent_km2)),
                                  family = gaussian(), data = oc_2yr,
                                  chains = 4, cores = 4, iter = 4000, thin = 4,
                                  control = list(adapt_delta = 0.9, max_treedepth = 10), 
                                  file = "models/brm_oc_tmeanpos_2yr")
  print(brm_oc_tmeanpos_2yr)
  t2 <- Sys.time()
  print(t2 - t1)
  
  brm_ac_tmeanpos_2yr <- brm(bf(N_change_log ~ 0 + realm + tmeanpos:realm + (1 + tmeanpos|study), 
                                     sigma ~ log(n_samp) + log(duration) + log(extent_km2)),
                                  family = gaussian(), data = ac_2yr,
                                  chains = 4, cores = 4, iter = 4000, thin = 4,
                                  control = list(adapt_delta = 0.9, max_treedepth = 10), 
                                  file = "models/brm_ac_tmeanpos_2yr")
  print(brm_ac_tmeanpos_2yr)
  t3 <- Sys.time()
  print(t3 - t2)
  
  # test the interaction between baseline thermal position and temperature change
  brm_oc_tmeanpos_tmeanchange_2yr <- brm(bf(occup_change_logit ~ 0 + realm +  (tmeanpos*change_tempmean):realm + (1 + tmeanpos|study), 
                                                 sigma ~ log(n_samp) + log(duration) + log(extent_km2)),
                                              family = gaussian(), data = oc_2yr,
                                              chains = 4, cores = 4, iter = 4000, thin = 4,
                                              control = list(adapt_delta = 0.9, max_treedepth = 10), 
                                              file = "models/brm_oc_tmeanpos_tmeanchange_2yr")
  print(brm_oc_tmeanpos_tmeanchange_2yr)
  t4 <- Sys.time()
  print(t4 - t3)
  
  brm_ac_tmeanpos_tmeanchange_2yr <- brm(bf(N_change_log ~ 0 + realm +  (tmeanpos*change_tempmean):realm + (1 + tmeanpos|study), 
                                                 sigma ~ log(n_samp) + log(duration) + log(extent_km2)),
                                              family = gaussian(), data = ac_2yr,
                                              chains = 4, cores = 4, iter = 4000, thin = 4,
                                              control = list(adapt_delta = 0.9, max_treedepth = 10), 
                                              file = "models/brm_ac_tmeanpos_tmeanchange_2yr")
  print(brm_ac_tmeanpos_tmeanchange_2yr)
  t5 <- Sys.time()
  print(t5 - t4)
}


## sensitivity analyses: use the temperature at the sampling years
if(index == 14){
  # test the overall effects of baseline thermal position
  t1 <- Sys.time()
  brm_oc_tmeanpos_sampyr <- brm(bf(occup_change_logit ~ 0 + realm + tmeanpos_sampyr:realm + (1 + tmeanpos_sampyr|study), 
                                     sigma ~ log(n_samp) + log(duration) + log(extent_km2)),
                                  family = gaussian(), data = oc_period,
                                  chains = 4, cores = 4, iter = 4000, thin = 4,
                                  control = list(adapt_delta = 0.9, max_treedepth = 10), 
                                  file = "models/brm_oc_tmeanpos_sampyr")
  print(brm_oc_tmeanpos_sampyr)
  t2 <- Sys.time()
  print(t2 - t1)
  
  brm_ac_tmeanpos_sampyr <- brm(bf(N_change_log ~ 0 + realm + tmeanpos_sampyr:realm + (1 + tmeanpos_sampyr|study), 
                                     sigma ~ log(n_samp) + log(duration) + log(extent_km2)),
                                  family = gaussian(), data = ac_period,
                                  chains = 4, cores = 4, iter = 4000, thin = 4,
                                  control = list(adapt_delta = 0.9, max_treedepth = 10), 
                                  file = "models/brm_ac_tmeanpos_sampyr")
  print(brm_ac_tmeanpos_sampyr)
  t3 <- Sys.time()
  print(t3 - t2)
  
  # test the interaction between baseline thermal position and temperature change
  brm_oc_tmeanpos_tmeanchange_sampyr <- brm(bf(occup_change_logit ~ 0 + realm +  (tmeanpos_sampyr*change_tempmean_sampyr):realm + (1 + tmeanpos_sampyr|study), 
                                                 sigma ~ log(n_samp) + log(duration) + log(extent_km2)),
                                              family = gaussian(), data = oc_period,
                                              chains = 4, cores = 4, iter = 4000, thin = 4,
                                              control = list(adapt_delta = 0.9, max_treedepth = 10), 
                                              file = "models/brm_oc_tmeanpos_tmeanchange_sampyr")
  print(brm_oc_tmeanpos_tmeanchange_sampyr)
  t4 <- Sys.time()
  print(t4 - t3)
  
  brm_ac_tmeanpos_tmeanchange_sampyr <- brm(bf(N_change_log ~ 0 + realm +  (tmeanpos_sampyr*change_tempmean_sampyr):realm + (1 + tmeanpos_sampyr|study), 
                                                 sigma ~ log(n_samp) + log(duration) + log(extent_km2)),
                                              family = gaussian(), data = ac_period,
                                              chains = 4, cores = 4, iter = 4000, thin = 4,
                                              control = list(adapt_delta = 0.9, max_treedepth = 10), 
                                              file = "models/brm_ac_tmeanpos_tmeanchange_sampyr")
  print(brm_ac_tmeanpos_tmeanchange_sampyr)
  t5 <- Sys.time()
  print(t5 - t4)
}


## sensitivity analyses: use the thermal  position estimated as 2.5% and 97.5% quantile of observed temperatures
if(index == 15){
  # test the overall effects of baseline thermal position
  t1 <- Sys.time()
  brm_oc_tmeanpos98_realm <- brm(bf(occup_change_logit ~ 0 + realm + tmeanpos98:realm + (1 + tmeanpos98|study), 
                                  sigma ~ log(n_samp) + log(duration) + log(extent_km2)),
                               family = gaussian(), data = oc_period ,
                               chains = 4, cores = 4, iter = 4000, thin = 4,
                               control = list(adapt_delta = 0.9, max_treedepth = 10), 
                               file = "models/brm_oc_tmeanpos98_realm")
  print(brm_oc_tmeanpos98_realm)
  t2 <- Sys.time()
  print(t2 - t1)
  
  brm_ac_tmeanpos98_realm <- brm(bf(N_change_log ~ 0 + realm + tmeanpos98:realm + (1 + tmeanpos98|study), 
                                    sigma ~ log(n_samp) + log(duration) + log(extent_km2)),
                               family = gaussian(), data = ac_period ,
                               chains = 4, cores = 4, iter = 4000, thin = 4,
                               control = list(adapt_delta = 0.9, max_treedepth = 10), 
                               file = "models/brm_ac_tmeanpos98_realm")
  print(brm_ac_tmeanpos98_realm)
  t3 <- Sys.time()
  print(t3 - t2)
  

  # test the interaction between baseline thermal position and temperature change
  brm_oc_tmeanpos98_tmeanchange_realm <- brm(bf(occup_change_logit ~ 0 + realm +  (tmeanpos98*change_tempmean):realm + (1 + tmeanpos98|study), 
                                              sigma ~ log(n_samp) + log(duration) + log(extent_km2)),
                                           family = gaussian(), data = oc_period ,
                                           chains = 4, cores = 4, iter = 4000, thin = 4,
                                           control = list(adapt_delta = 0.9, max_treedepth = 10), 
                                           file = "models/brm_oc_tmeanpos98_tmeanchange_realm")
  print(brm_oc_tmeanpos98_tmeanchange_realm)
  t4 <- Sys.time()
  print(t4 - t3)
  
  brm_ac_tmeanpos98_tmeanchange_realm <- brm(bf(N_change_log ~ 0 + realm +  (tmeanpos98*change_tempmean):realm + (1 + tmeanpos98|study), 
                                                sigma ~ log(n_samp) + log(duration) + log(extent_km2)),
                                           family = gaussian(), data = ac_period ,
                                           chains = 4, cores = 4, iter = 4000, thin = 4,
                                           control = list(adapt_delta = 0.9, max_treedepth = 10), 
                                           file = "models/brm_ac_tmeanpos98_tmeanchange_realm")
  print(brm_ac_tmeanpos98_tmeanchange_realm)
  t5 <- Sys.time()
  print(t5 - t4)
}


## sensitivity analyses: use the data with abundance data for occupancy change analyses
if(index == 16){
  # test the overall effects of baseline thermal position
  t1 <- Sys.time()
  brm_oc_tmeanpos_abundancedata <- brm(bf(occup_change_logit ~ 0 + realm + tmeanpos:realm + (1 + tmeanpos|study), 
                                     sigma ~ log(n_samp) + log(duration) + log(extent_km2)),
                                  family = gaussian(), data = oc_period %>% filter(subset_abundance),
                                  chains = 4, cores = 4, iter = 4000, thin = 4,
                                  control = list(adapt_delta = 0.9, max_treedepth = 10), 
                                  file = "models/brm_oc_tmeanpos_abundancedata")
  print(brm_oc_tmeanpos_abundancedata)
  t2 <- Sys.time()
  print(t2 - t1)
  
  
  # test the interaction between baseline thermal position and temperature change
  brm_oc_tmeanpos_tmeanchange_abundancedata <- brm(bf(occup_change_logit ~ 0 + realm +  (tmeanpos*change_tempmean):realm + (1 + tmeanpos|study), 
                                                 sigma ~ log(n_samp) + log(duration) + log(extent_km2)),
                                              family = gaussian(), data = oc_period %>% filter(subset_abundance),
                                              chains = 4, cores = 4, iter = 4000, thin = 4,
                                              control = list(adapt_delta = 0.9, max_treedepth = 10), 
                                              file = "models/brm_oc_tmeanpos_tmeanchange_abundancedata")
  print(brm_oc_tmeanpos_tmeanchange_abundancedata)
  t3 <- Sys.time()
  print(t3 - t2)
}


## sensitivity analyses: use the meta-communities with duration at least 20 years
if(index == 17){
  # test the overall effects of baseline thermal position
  t1 <- Sys.time()
  brm_oc_tmeanpos_longduration <- brm(bf(occup_change_logit ~ 0 + realm + tmeanpos:realm + (1 + tmeanpos|study), 
                                     sigma ~ log(n_samp) + log(duration) + log(extent_km2)),
                                  family = gaussian(), data = oc_period %>% filter(duration >= 20),
                                  chains = 4, cores = 4, iter = 4000, thin = 4,
                                  control = list(adapt_delta = 0.9, max_treedepth = 10), 
                                  file = "models/brm_oc_tmeanpos_longduration")
  print(brm_oc_tmeanpos_longduration)
  t2 <- Sys.time()
  print(t2 - t1)
  
  brm_ac_tmeanpos_longduration <- brm(bf(N_change_log ~ 0 + realm + tmeanpos:realm + (1 + tmeanpos|study), 
                                     sigma ~ log(n_samp) + log(duration) + log(extent_km2)),
                                  family = gaussian(), data = ac_period %>% filter(duration >= 20),
                                  chains = 4, cores = 4, iter = 4000, thin = 4,
                                  control = list(adapt_delta = 0.9, max_treedepth = 10), 
                                  file = "models/brm_ac_tmeanpos_longduration")
  print(brm_ac_tmeanpos_longduration)
  t3 <- Sys.time()
  print(t3 - t2)
  
  # test the interaction between baseline thermal position and temperature change
  brm_oc_tmeanpos_tmeanchange_longduration <- brm(bf(occup_change_logit ~ 0 + realm +  (tmeanpos*change_tempmean):realm + (1 + tmeanpos|study), 
                                                 sigma ~ log(n_samp) + log(duration) + log(extent_km2)),
                                              family = gaussian(), data = oc_period %>% filter(duration >= 20),
                                              chains = 4, cores = 4, iter = 4000, thin = 4,
                                              control = list(adapt_delta = 0.9, max_treedepth = 10), 
                                              file = "models/brm_oc_tmeanpos_tmeanchange_longduration")
  print(brm_oc_tmeanpos_tmeanchange_longduration)
  t4 <- Sys.time()
  print(t4 - t3)
  
  brm_ac_tmeanpos_tmeanchange_longduration <- brm(bf(N_change_log ~ 0 + realm +  (tmeanpos*change_tempmean):realm + (1 + tmeanpos|study), 
                                                 sigma ~ log(n_samp) + log(duration) + log(extent_km2)),
                                              family = gaussian(), data = ac_period %>% filter(duration >= 20),
                                              chains = 4, cores = 4, iter = 4000, thin = 4,
                                              control = list(adapt_delta = 0.9, max_treedepth = 10), 
                                              file = "models/brm_ac_tmeanpos_tmeanchange_longduration")
  print(brm_ac_tmeanpos_tmeanchange_longduration)
  t5 <- Sys.time()
  print(t5 - t4)
}
