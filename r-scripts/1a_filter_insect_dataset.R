## filter InsectChange database to choose subsets of datasets with number of locations in each year >=4, duration >=10, 
## and matching locations across years to keep keep sites in same locations

rm(list = ls())
setwd("C:/Dropbox/iDiv/Populations_warming/analyses_202406")

library(tidyverse)
library(vegan)
library(reshape2)
library(magrittr)

ic <- readRDS("data/InsectChange/Raw data insect metacommunities.Rdata")

# the additional insect studies after the publication of InsectChange
load("data/InsectChange/mosquito_enc_standardized.RDATA")


# to combine these new studies with the main insect studies, change the column names
mosquito_enc_std <- mosquito_enc_std %>% 
  rename(Datasource_name = study, Plot_ID = plot, Latitude = Latitudes, Longitude = Longitudes, 
         Year = year, species = Species, Reference = Citations, Datasource_nameREDUNDANT = Tag , Number = N, Realm = realm) %>%
  mutate(Abundance.Biomass = "A",
         Original_unit = "total number per year",
         Unit.y = "abundance") %>%
  relocate(Datasource_ID)

mosquito_enc_std %>% distinct(Datasource_ID, Datasource_name)


# combine insect datasets, and and filtered out datasets with < 4 locations
ic_4loc <- ic %>% 
  as_tibble() %>% 
  filter(Level == "Species") %>% 
  mutate(species = paste(Genus, Species, sep = " ")) %>% 
  # add the new studies
  mutate(Datasource_ID = as.character(Datasource_ID),
         Plot_ID = as.character(Plot_ID)) %>%
  bind_rows(mosquito_enc_std %>% 
              mutate(Datasource_ID = as.character(Datasource_ID))) %>%
  group_by(Datasource_ID) %>%
  mutate(n_loc = n_distinct(Plot_ID)) %>% 
  ungroup() %>%
  filter(n_loc > 3) %>%
  dplyr::select(!n_loc)

# the study 1367 survey in 1980 and 1981, and resurvey in 2008 and 2009; change the years as two periods
ic_4loc %>% filter(Datasource_ID == "1367") %>% dplyr::select(Datasource_ID, Plot_ID, Year)
ic_4loc %>% filter(Datasource_ID == "1367") %>% pull(Year) %>% table()
id1 <- ic_4loc$Datasource_ID == "1367" & ic_4loc$Year == 1981
id2 <- ic_4loc$Datasource_ID == "1367" & ic_4loc$Year == 2009
ic_4loc$Year[id1] <- 1980
ic_4loc$Year[id2] <- 2008


# check whether different plots have same coordinates.
# only few plots have same coordinates with the others; use the Plot_ID to indicate locations (sites)
ic_4loc %>% 
  distinct(Datasource_ID, Plot_ID, Longitude, Latitude) %>%
  filter(duplicated(.[, c("Longitude", "Latitude")]))

# 2 studies (ID = 1367, 1408) provide only regional central coordinates
ic_4loc %>% 
  group_by(Datasource_ID) %>%
  summarise(n_loc = n_distinct(Longitude, Latitude)) %>%
  filter(n_loc == 1)


###########################
## Check abundance and biomass information
x <- ic_4loc %>% select(Datasource_ID, Plot_ID, Number, Biomass, Abundance.Biomass, Taxon, species, Level, Rank, Original_unit, Unit.y)
table(x$Original_unit, useNA = "always")
table(x$Unit.y, useNA = "always") # abundance,   density
table(x$Abundance.Biomass, useNA = "always") # A, AB

x %>% 
  group_by(Datasource_ID, Plot_ID) %>%
  summarise(range_abundance = diff(range(Number))) %>%
  filter(range_abundance == 0)

x %>% 
  group_by(Datasource_ID) %>%
  summarise(range_abundance = diff(range(Number))) %>%
  filter(range_abundance == 0)

range(ic_4loc$Number, na.rm = TRUE)
table(ic_4loc$Number == 0, useNA = "always")

# check datasets with zero abundance values 
ic_4loc %>% filter(Number == 0) %>% pull(Datasource_ID) %>% table()
x %>% filter(Datasource_ID == 1533) %$% table(Number == 0)
x %>% filter(Datasource_ID == 1533 & Number == 0) %$% table(Biomass, useNA = "always")
x %>% filter(Datasource_ID == 502) %$% table(Number == 0)
x %>% filter(Datasource_ID == 502 & Number == 0) %$% table(Biomass, useNA = "always")
x %>% filter(Datasource_ID == 79) %$% table(Number == 0)
x %>% filter(Datasource_ID == 79 & Number == 0) %$% table(Biomass, useNA = "always")

y <- x %>% filter(Datasource_ID == 1533); 
View(y)
ic_4loc %>% filter(Datasource_ID == 1533) %>% 
  mutate(is.zero.abundance = Number == 0) %>% 
  group_by(Plot_ID, Year, is.zero.abundance) %>% 
  summarise(sprich = n_distinct(species)) # most records are zero; exclude this dataset

y <- x %>% filter(Datasource_ID == 502); 
View(y)
ic_4loc %>% filter(Datasource_ID == 502) %>% 
  mutate(is.zero.abundance = Number == 0) %>% 
  group_by(Plot_ID, is.zero.abundance) %>% 
  summarise(sprich = n_distinct(species))

y <- x %>% filter(Datasource_ID == 79); 
View(y)
ic_4loc %>% filter(Datasource_ID == 79) %>% 
  mutate(is.zero.abundance = Number == 0) %>% 
  group_by(Plot_ID, is.zero.abundance) %>% 
  summarise(sprich = n_distinct(species))

# check dataset with abundance values as NA 
ic_4loc %>% filter(is.na(Number)) %>% pull(Datasource_ID) %>% table()
x %>% filter(Datasource_ID == 1542) %$% table(is.na(Number))
y <- x %>% filter(Datasource_ID == 1542); 
View(y)

# note: use three columns for abundance information: Number, Abundance.Biomass, Unit.y
# some "Number" values are 0, meaning true 0 count (absence).For dataset 1533, most records are zero; exclude this dataset 
# some "Number" values are NA


############
## the number and locations of sites are somewhat different across years
# filter the dataset to maxiamize the number of same samples for at least two years within a region

ic_4loc_filtered <- NULL
for(i in 1:length(unique(ic_4loc$Datasource_ID))){
  # perform loop for each study
  study <- ic_4loc %>% 
    filter(Datasource_ID == unique(ic_4loc$Datasource_ID)[i])
  
  # get a matrix with row are years and colunmns is samples
  year_sample <- as.matrix(xtabs(~ Year + Plot_ID, data = unique(study[,c("Year","Plot_ID")]), sparse=TRUE))
  
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
    filter(Plot_ID %in% colnames(year_sample) & Year %in% rownames(year_sample))
  
  ic_4loc_filtered <- bind_rows(ic_4loc_filtered, study_filtered)
}

ggplot(data = study_filtered, aes(Longitude, Latitude)) + 
  facet_wrap(~Year) +
  geom_point(size = 1, alpha = 0.7) + 
  coord_fixed()


# keep only years with at least 4 plots,
# and keep studies with at least 2 time points and duration >10 years
ic_4loc_filtered <- ic_4loc_filtered %>%
  group_by(Datasource_ID, Year) %>%
  mutate(n_samp = n_distinct(Plot_ID)) %>%
  filter(n_samp > 3) %>%
  group_by(Datasource_ID) %>%
  mutate(n_samples = n_distinct(Plot_ID),
         n_years = n_distinct(Year),
         duration = max(Year) - min(Year) +1) %>%
  ungroup() %>%
  filter(n_years >= 2 & duration >= 10)

# check how many studies and their attributes
ic_studies <- ic_4loc_filtered %>% distinct(Datasource_ID, n_samples, n_years, duration) # 23 studies
table(ic_studies$n_samples)
table(ic_studies$n_years)  
table(ic_studies$duration)

# save the filtered data
ic_filtered <- ic_4loc_filtered %>% dplyr::select(-(n_samp:duration)) 
save(ic_filtered, file = "data/InsectChange/InsectChange_filtered.RDATA")


# locations of the filtered dataset
ic_loc_filtered <- ic_4loc_filtered %>% 
  distinct(Datasource_ID, Plot_ID, Year, Latitude, Longitude)

# locations of the full dataset
ic_loc <- ic_4loc %>% 
  distinct(Datasource_ID, Plot_ID, Year, Latitude, Longitude ) %>% 
  left_join(ic_loc_filtered %>%
              mutate(keep = "yes")) %>%
  mutate(keep = ifelse(is.na(keep), "no", keep))

# plot distributions of samples and indicate which records will be removed
pdf('data/InsectChange/InsectChange_filtered.pdf', width = 12, height = 10)
id_study <- unique(ic_loc$Datasource_ID)
for(i in 1:length(id_study)){
  study <- ic_loc %>% 
    filter(Datasource_ID %in% id_study[i])
  
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

