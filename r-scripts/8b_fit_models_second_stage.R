######### extract study-level slopes of thermal position from first-stage models and fit models assessing the effects of temperature change and other variables on the slopes

rm(list = ls())

# load packages
needed_libs <- c("tidyverse","brms")

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


## input the observation data and and first-stage models 
load("results/data_input_to_models.RDATA")
brm_oc_tmeanpos <- readRDS("models/brm_oc_tmeanpos.rds")
brm_ac_tmeanpos <- readRDS("models/brm_ac_tmeanpos.rds")


#######################
## extract study-level slopes of thermal position from first-stage models

## study-level effects of thermal position on occupancy change
coef_oc_tmeanpos <- coef(brm_oc_tmeanpos, robust = TRUE, probs = c(0.025, 0.975))[[1]][,,1:2] 
coef_oc_tmeanpos <- as_tibble(coef_oc_tmeanpos) %>%
  mutate(study = rownames(coef_oc_tmeanpos)) 
colnames(coef_oc_tmeanpos)[1:8] <- c("estimate_intercept", "se_intercept", "Q2.5_intercept", "Q97.5_intercept", 
                                     "estimate_slope", "se_slope", "Q2.5_slope", "Q97.5_slope")

# add study-level meta data
coef_oc_tmeanpos <- coef_oc_tmeanpos %>% 
  # indicate significance of slopes
  mutate(sig_slope = ifelse(Q2.5_slope < 0 & Q97.5_slope >0, "neutral", ifelse(Q97.5_slope <= 0, "negative", "positive")),
         sig_slope = factor(sig_slope, levels = c("negative", "positive", "neutral"))) %>%
  left_join(oc_period %>% group_by(study) %>%
              summarise(xmin = min(tmeanpos),
                        xmax = max(tmeanpos))) %>%
  left_join(oc_period %>% select(c(study:n_samp, duration_mean, taxon_new:change_tempmean, realm_taxa, realm_climate)) %>% distinct()) %>% 
  left_join(dat_meta %>% dplyr::select(study, cent_lat, cent_long)) %>%
  relocate(study, database, studyID, study_name)

## fixed effect of thermal position
fixef_oc_tmeanpos <- fixef(brm_oc_tmeanpos) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "term") %>%
  filter(grepl("tmeanpos", term))


## study-level effects of thermal position on abundance change
coef_ac_tmeanpos <- coef(brm_ac_tmeanpos, robust = TRUE, probs = c(0.025, 0.975))[[1]][,,1:2] 
coef_ac_tmeanpos <- as_tibble(coef_ac_tmeanpos) %>%
  mutate(study = rownames(coef_ac_tmeanpos)) 
colnames(coef_ac_tmeanpos)[1:8] <- c("estimate_intercept", "se_intercept", "Q2.5_intercept", "Q97.5_intercept", 
                                     "estimate_slope", "se_slope", "Q2.5_slope", "Q97.5_slope")

# add study-level meta data
coef_ac_tmeanpos <- coef_ac_tmeanpos %>% 
  # indicate significance of slopes
  mutate(sig_slope = ifelse(Q2.5_slope < 0 & Q97.5_slope >0, "neutral", ifelse(Q97.5_slope <= 0, "negative", "positive")),
         sig_slope = factor(sig_slope, levels = c("negative", "positive", "neutral"))) %>%
  left_join(ac_period %>% group_by(study) %>%
              summarise(xmin = min(tmeanpos),
                        xmax = max(tmeanpos))) %>%
  left_join(ac_period %>% select(c(study:study_name, sprich:n_samp, duration_mean, taxon_new:change_tempmean, realm_taxa, realm_climate)) %>% distinct()) %>% 
  left_join(dat_meta %>% dplyr::select(study, cent_lat, cent_long)) %>%
  relocate(study, database, studyID, study_name)

## fixed effect of thermal position
fixef_ac_tmeanpos <- fixef(brm_ac_tmeanpos) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "term") %>%
  filter(grepl("tmeanpos", term))


save(coef_oc_tmeanpos, coef_ac_tmeanpos, fixef_oc_tmeanpos, fixef_ac_tmeanpos, file = "results/coef_oc_ac_tmeanpos.RDATA")


########################################
## fit models assessing the effects of temperature change and other variables on the slopes
load("results/coef_oc_ac_tmeanpos.RDATA")

## test the relationship between regional-level slopes of thermal position and temperature change
t1 <- Sys.time()
brm_slope_oc_tmeanchange_realm <- brm(bf(estimate_slope|se(se_slope) ~ 0 + realm + change_tempmean:realm + (1|study)),
                                      data = coef_oc_tmeanpos,
                                      control = list(adapt_delta = 0.9, max_treedepth = 10),
                                      cores = 4, chains = 4, iter = 40000, thin = 10,
                                      file = "models/brm_slope_oc_tmeanchange_realm")
print(brm_slope_oc_tmeanchange_realm)
t2 <- Sys.time()
print(t2 - t1)

brm_slope_ac_tmeanchange_realm <- brm(bf(estimate_slope|se(se_slope) ~ 0 + realm + change_tempmean:realm + (1|study)),
                                      data = coef_ac_tmeanpos,
                                      control = list(adapt_delta = 0.9, max_treedepth = 10),
                                      cores = 4, chains = 4, iter = 40000, thin = 10,
                                      file = "models/brm_slope_ac_tmeanchange_realm")
print(brm_slope_ac_tmeanchange_realm)
t3 <- Sys.time()
print(t3 - t2)


## test the relationship between regional-level slopes of thermal position and temperature change and other covariables
coef_oc_tmeanpos <- coef_oc_tmeanpos %>% 
  mutate(lat_cent = abs(cent_lat) - mean(abs(cent_lat)),
         sprich_cent = log10(sprich) - mean(log10(sprich)),       
         duration_cent = log10(duration_mean) - mean(log10(duration_mean)),
         startyr_cent = start_year - mean(start_year),
         extent_cent = log10(extent_km2) - mean(log10(extent_km2)),
         nsamp_cent = log10(n_samp) - mean(log10(n_samp)))

coef_ac_tmeanpos <- coef_ac_tmeanpos %>% 
  mutate(lat_cent = abs(cent_lat) - mean(abs(cent_lat)),
         sprich_cent = log10(sprich) - mean(log10(sprich)),       
         duration_cent = log10(duration_mean) - mean(log10(duration_mean)),
         startyr_cent = start_year - mean(start_year),
         extent_cent = log10(extent_km2) - mean(log10(extent_km2)),
         nsamp_cent = log10(n_samp) - mean(log10(n_samp)))

t1 <- Sys.time()
brm_slope_oc_covariates_realm <- brm(bf(estimate_slope|se(se_slope) ~ 
                                          0 + realm + (change_tempmean + lat_cent + sprich_cent + duration_cent + startyr_cent + 
                                                         extent_cent + nsamp_cent):realm + (1|study)),
                                     data = coef_oc_tmeanpos,
                                     control = list(adapt_delta = 0.9, max_treedepth = 10),
                                     cores = 4, chains = 4, iter = 40000, thin = 10,
                                     file = "models/brm_slope_oc_covariates_realm")
print(brm_slope_oc_covariates_realm) # nsamp_cent singnificant
t2 <- Sys.time()
print(t2 - t1)

brm_slope_ac_covariates_realm <- brm(bf(estimate_slope|se(se_slope) ~ 
                                          0 + realm + (change_tempmean + lat_cent + sprich_cent + duration_cent + startyr_cent + 
                                                         extent_cent + nsamp_cent):realm + (1|study)),
                                     data = coef_ac_tmeanpos,
                                     control = list(adapt_delta = 0.9, max_treedepth = 10),
                                     cores = 4, chains = 4, iter = 40000, thin = 10,
                                     file = "models/brm_slope_ac_covariates_realm")
print(brm_slope_ac_covariates_realm)
t3 <- Sys.time()
print(t3 - t2)
