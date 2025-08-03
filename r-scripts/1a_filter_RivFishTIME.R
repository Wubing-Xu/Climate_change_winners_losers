## filter RivFishTIME database to choose subsets of datasets with number of locations in each year >=4, duration >=10, 
## and matching locations across years to keep keep sites in same locations

rm(list = ls())
setwd("C:/Dropbox/iDiv/Populations_warming/analyses_202406")

library(tidyverse)
library(vegan)
library(reshape2)
library(magrittr)

ft <- read_csv("data/RivFishTIME/1873_10_1873_2_RivFishTIME_SurveyTable.csv")
meta <- read_csv('data/RivFishTIME/1873_10_1873_2_RivFishTIME_TimeseriesTable.csv')[,-13]

# studies with at least 4 time series
study_4loc <- meta %>%
  group_by(SourceID) %>%
  summarise(n_loc = n_distinct(TimeSeriesID),
            n_protocol = n_distinct(Protocol)) %>% 
  ungroup() %>% 
  filter(n_loc > 3)

# 10 studies with multiple protocols
study_mp <- study_4loc  %>% 
  filter(n_protocol > 1)

# Choose the protocol with the maximum number of time series
study_mp_selected <- meta %>% 
  filter(SourceID %in% pull(study_mp, SourceID)) %>%
  group_by(SourceID, Protocol) %>%
  summarise(n_loc = n_distinct(TimeSeriesID)) %>%
  filter(n_loc > 3) %>%
  group_by(SourceID) %>%
  filter(n_loc == max(n_loc)) %>%
  ungroup() %>% 
  distinct(SourceID, .keep_all = TRUE) #one study have 2 protocol with same number of time series. Keep one of them.

# the studies with at least 4 time series in the same protocol 
meta_4loc <- bind_rows(meta %>% filter(SourceID %in% 
                                         (study_4loc %>% filter(n_protocol == 1) %>% pull(SourceID))),
                       meta %>% inner_join(study_mp_selected %>% dplyr::select(-n_loc)))

# check number of sites and time series for each study, and whether they are equal
meta_4loc_nsites <- meta_4loc %>% 
  group_by(SourceID) %>% 
  summarise(n_site = n_distinct(SiteID),
            n_ts = n_distinct(TimeSeriesID))
table(meta_4loc_nsites$n_site == meta_4loc_nsites$n_ts) # all 39 studies is TRUE, meaning each time series has unique site

# communities for the selected time series
ft_4loc <- inner_join(ft, meta_4loc, by = "TimeSeriesID") %>% 
  relocate(SourceID)


###################
## Check abundance  information
table(ft_4loc$UnitAbundance)

ft_4loc %>% filter(UnitAbundance == "Leslie_index") %$% hist(Abundance)
ft_4loc %>% filter(UnitAbundance == "Leslie_index") %>% distinct(SourceID)

ft_4loc %>% filter(UnitAbundance == "Moyle_class") %$% hist(Abundance)
ft_4loc %>% filter(UnitAbundance == "Moyle_class") %>% distinct(SourceID)
ft_4loc %>% filter(UnitAbundance == "Moyle_class") %$% table(Abundance)

# check datasets with small min and max abundance values 
ft_study_abundance_range <- ft_4loc %>% 
  group_by(SourceID, UnitAbundance) %>%
  summarise(max_abundance = max(range(Abundance)),
            min_abundance = min(range(Abundance)),
            range_abundance = diff(range(Abundance))) 

ft_study_abundance_range %>% filter(max_abundance <= 10 | min_abundance < 1)
ft_4loc %>% filter(SourceID == 16) %>% ggplot(aes(x = Abundance)) + geom_histogram()

# check whether some studies have multiple abundance units within each study
ft_4loc %>%
  group_by(SourceID) %>%
  summarise(n_unitAbundance = n_distinct(UnitAbundance)) %>% 
  filter(n_unitAbundance > 1)

# 2 studies have 2 units within each study. Choose the one with the maximum number of time series
ft_4loc %>%
  filter(SourceID %in% c(27, 32)) %>%
  distinct(SourceID, TimeSeriesID, UnitAbundance) %>%
  group_by(SourceID, UnitAbundance) %>%
  summarise(n_distinct(TimeSeriesID))

ft_4loc %>% filter(SourceID == 27) %$% table(UnitAbundance)
ft_4loc %>% filter(SourceID == 32) %$% table(UnitAbundance)

ft_4loc <- ft_4loc %>% 
  filter(!(SourceID == 27 & UnitAbundance == "Count")) %>%
  filter(!(SourceID == 32 & UnitAbundance == "Ind.100m2"))

# note: use UnitAbundance values as "Count" and "CPUE", "Ind.100m2", "Leslie_index" (delete "Moyle_class") 


########################
## the number and locations of sites are somewhat different across years
# filter the dataset to maximum the number of same samples for at least two years within a region

ft_4loc_filtered <- NULL
for(i in 1:length(unique(ft_4loc$SourceID))){
  # perform loop for each study
  study <- ft_4loc %>% 
    filter(SourceID == unique(ft_4loc$SourceID)[i])
  
  # get a matrix with row are years and colunmns is samples
  year_sample <- as.matrix(xtabs(~ Year + TimeSeriesID, data = unique(study[,c("Year","TimeSeriesID")]), sparse=TRUE))
  
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
    filter(TimeSeriesID %in% colnames(year_sample) & Year %in% rownames(year_sample))
  
  ft_4loc_filtered <- bind_rows(ft_4loc_filtered, study_filtered)
}


ggplot(data = study_filtered, aes(Longitude, Latitude)) + 
  facet_wrap(~Year) +
  geom_point(size = 1, alpha = 0.7) + 
  coord_fixed()


# keep only years with at least 4 sites,
# and keep studies with at least 2 time points and duration >10 years
ft_4loc_filtered <- ft_4loc_filtered %>%
  group_by(SourceID, Year) %>%
  mutate(n_site = n_distinct(TimeSeriesID)) %>%
  filter(n_site > 3) %>%
  group_by(SourceID) %>%
  mutate(n_samples = n_distinct(TimeSeriesID ),
         n_years = n_distinct(Year),
         duration = max(Year) - min(Year) +1) %>%
  ungroup() %>%
  filter(n_years >= 2 & duration >= 10)


# check how many studies and their attributes
ft_studies <- ft_4loc_filtered %>% distinct(SourceID, n_samples, n_years, duration) #35 studies
table(ft_studies$n_samples) # 24 studies >= 10
table(ft_studies$n_years)  # 21 studies <= 3; 8 studies >=10 
table(ft_studies$duration)


# add the column "keep" to distinguish the records that should be kept or removed
ft_4loc <- ft_4loc %>% 
  left_join(ft_4loc_filtered %>% 
              dplyr::select(- (n_site:duration)) %>%
              mutate(keep = "yes")) %>%
  mutate(keep = ifelse(is.na(keep), "no", keep))


# plot distributions of sites and indicate which records will be removed
pdf('data/RivFishTIME/RivFishTIME_4locations_filtered.pdf', width = 12, height = 10)
id_study <- unique(ft_4loc$SourceID)
for(i in 1:length(id_study)){
  study <- ft_4loc %>% 
    filter(SourceID %in% id_study[i]) %>% 
    distinct(SourceID, TimeSeriesID, Year, Latitude, Longitude, keep)
  
  p <- ggplot(data = study, aes(Longitude, Latitude)) + 
    facet_wrap(~Year) +
    geom_point( aes(colour = keep), size = 1, alpha = 0.7) + 
    labs(title = id_study[i]) +
    theme(legend.position = "top", legend.text = element_text(size=12)) + 
    coord_fixed() +
    scale_color_manual(values=c("yes" = "deepskyblue", "no" = "coral"))
  
  print(p)
}
dev.off()


# check how many time series have different number of surveys across years
n_survey <- ft_4loc_filtered %>% 
  dplyr::select(-(n_site:duration)) %>%
  group_by(SourceID, TimeSeriesID, Year) %>%
  summarise(n_survey = n_distinct(SurveyID)) %>%
  group_by(SourceID, TimeSeriesID) %>% 
  summarise(min_survey = min(n_survey),
            max_survey = max(n_survey))

# most time series have same number of surveys between years
table(n_survey$min_survey == n_survey$max_survey)

# only 34 time series have multiple surveys at all year
table(n_survey$min_survey > 1) 

# no studies have more than 1 surveys through sites and years
n_survey %>% group_by(SourceID) %>% 
  summarise(min_survey_years = min(min_survey),
            max_survey_years = min(max_survey)) %>%
  filter(min_survey_years > 1)

# keep one survey for each year of each time series
one_survey_perYear <- ft_4loc_filtered %>% 
  dplyr::select(SourceID, TimeSeriesID, Year, SurveyID) %>%
  distinct(SourceID, TimeSeriesID, Year, .keep_all=TRUE)

ft_4loc_filtered <-  ft_4loc_filtered %>% inner_join(one_survey_perYear)

# save the filtered data
ft_filtered <- ft_4loc_filtered %>% dplyr::select(-(n_site:duration)) 
save(ft_filtered, file = "data/RivFishTIME/RivFishTIME_4sites_10years_filtered.RDATA")

