## filter BioTIME database to choose subsets of datasets with number of locations in each year >=4, duration >=10, 
## and matching locations across years to keep keep sites in same locations

rm(list = ls())
setwd("C:/Dropbox/iDiv/Populations_warming/analyses_202406")

library(tidyverse)
library(vegan)
library(reshape2)
library(magrittr)


###############
## Choose subsets with number of locations in each year >=4, duration >=10, and number of years >2

load("data/BioTIME/biotime_2021.RDATA")

bt_loc <- bt %>% dplyr::select(STUDY_ID, YEAR, LATITUDE, LONGITUDE) %>% distinct()

# remove years with number of locations < 4.
bt_loc <-  bt_loc %>% 
  group_by(STUDY_ID, YEAR) %>%
  mutate(n_loc = n_distinct(LATITUDE, LONGITUDE)) %>%
  ungroup() %>%
  filter(n_loc >= 4) %>%
  dplyr::select(-n_loc)

# calculate number of locations, number of years and duration
# and remove studies with number of years < 2 and duration <10 years
meta_year <- bt_loc %>%
  group_by(STUDY_ID, YEAR) %>%
  summarise(n_loc = n_distinct(LATITUDE, LONGITUDE, na.rm = TRUE)) %>%
  group_by(STUDY_ID) %>%
  mutate(t_loc = sum(n_loc), 
         mean_loc = round(mean(n_loc),1), 
         min_loc = min(n_loc), 
         max_loc = max(n_loc),
         n_years = n_distinct(YEAR, na.rm = TRUE),
         duration = max(YEAR) - min(YEAR) +1) %>% 
  ungroup() %>%
  filter(n_years >=2 & duration >=10)

meta <- meta_year %>% 
  dplyr::select(-c(YEAR, n_loc)) %>% 
  distinct()


# frequency distribution of number of time points, duration, locations in a study
ggplot(meta) + geom_histogram(aes(n_years)) +  theme(axis.text=element_text(size=14), axis.title=element_text(size=20))

ggplot(meta) + geom_histogram(aes(duration)) + theme(axis.text=element_text(size=14), axis.title=element_text(size=20))

ggplot(meta) + geom_histogram(aes(mean_loc)) + scale_x_log10() + theme(axis.text=element_text(size=14), axis.title=element_text(size=20))


# remove studies with number of years < 2 and duration <10 years
bt_4loc_10yr <- bt_loc %>% 
  unite(col = study_year, STUDY_ID, YEAR, remove = FALSE) %>%
  filter(study_year %in% (meta_year %>% unite(col = study_year, STUDY_ID, YEAR) %>% pull(study_year))) %>%
  dplyr::select(-study_year)


# BioTIME have no plot id to indicate resurvey sites for most studies. We will use coordinates to define sites. 
# As sites surveyed in different years may be note exactly same, particularly in the marine, we will round coordinates to 2 decimals (~1 km) 
# to match sites through years for studies with no same resurveyed locations (coordinates).
# We individually check site coordinates for each study. If some of their coordinates are same through years, we will not round them.

btstudy <- bt_4loc_10yr %>% 
  group_by(STUDY_ID) %>%
  summarise(lat_span = max(LATITUDE) - min(LATITUDE),
            long_span = max(LONGITUDE) - min(LONGITUDE))

bt1 <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[1]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =67, not to round 
View(bt1)
bt2 <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[2]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =68, to round, 1 decimals 
View(bt2)
bt3 <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[3]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =71, to round
View(bt3)
bt4 <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[4]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =78, not to round
View(bt4)
bt5 <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[5]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =81, to round
View(bt5)
bt6 <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[6]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =84, not to round
View(bt6)
bt7 <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[7]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =85, not to round
View(bt7)
bt8 <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[8]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =86, to round
View(bt8)
bt9 <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[9]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =90, to round
View(bt9)
bt10 <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[10]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =92, not to round?
View(bt10)
bt11 <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[11]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =97, not to round?
View(bt11)
bt12 <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[12]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =99, not to round?
View(bt12)``
bt13 <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[13]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =108, in 2 decimals
View(bt13)
bt14 <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[14]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =110, not to round
View(bt14)
bt15 <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[15]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =112, in 2 decimals
View(bt15)
bt16 <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[16]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =112, not to round?
View(bt16)
bt17 <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[17]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =117, in 2 decimals
View(bt17)
bt18 <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[18]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =119, to round
View(bt18)
bt19 <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[19]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =121, not to round
View(bt19)
bt20 <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[20]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =123, to round, 1 decimals
View(bt20)
bt21 <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[21]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =125, in 2 decimals
View(bt21)
bt22 <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[22]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =127, in 2 decimals
View(bt22)
bt23 <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[23]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =135, not to round?
View(bt23)
bt24 <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[24]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =147, to round
View(bt24)
bt25 <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[25]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =148, to round
View(bt25)
bt26 <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[26]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =152, in 1 decimals
View(bt26)
bt27 <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[27]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =162, to round
View(bt27)
bt28 <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[28]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =163, in 2 decimals
View(bt28)
bt29 <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[29]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =166, to round
View(bt29)
bt30<- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[30]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =169, to round
View(bt30)
bt31 <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[31]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =171, to round
View(bt31)
bt32 <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[32]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =172, to round
View(bt32)
bt33 <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[33]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =173, to round
View(bt33)
bt34 <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[34]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =176, to round
View(bt34)
bt35 <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[35]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =178, to round
View(bt35)
bt36 <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[36]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =180, in 2 decimals
View(bt36)
bt37 <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[37]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =182, to round
View(bt37)
bt38 <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[38]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =183, to round, 1 decimals
View(bt38)
bt39 <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[39]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =187, to round
View(bt39)
bt40 <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[40]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =189, to round
View(bt40)
bt41 <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[41]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =190, to round
View(bt41)
bt42 <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[42]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =191, to round
View(bt42)
bt43 <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[43]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =192, in 1 decimals
View(bt43)
bt44 <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[44]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =195, not to round
View(bt44)
bt45 <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[45]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =196, not to round
View(bt45)
bt46 <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[46]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =197, to round
View(bt46)
bt47 <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[47]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =198, to round
View(bt47)
bt48 <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[48]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =200, not to round
View(bt48)
bt49 <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[49]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =204, to round
View(bt49)
bt50 <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[50]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =206, to round
View(bt50)
bt51 <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[51]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =208, to round
View(bt51)
bt52 <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[52]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =209, to round
View(bt52)
bt53 <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[53]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =210, not to round
View(bt53)
bt54 <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[54]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =213,  to round
View(bt54)
bt55 <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[55]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =214,  not to round
View(bt55)
bt56 <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[56]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =215,  not to round
View(bt56)
bt57 <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[57]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =217,  not to round
View(bt57)
bt58 <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[58]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =219,  not to round
View(bt58)
bt59 <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[59]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =220,  not to round
View(bt59)
bt60 <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[60]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =221,  not to round
View(bt60)
bt61 <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[61]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =229,  not to round
View(bt61)
bt62 <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[62]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =232,  not to round?
View(bt62)
bt63 <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[63]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =237,  not to round
View(bt63)
bt64 <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[64]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =244,  not to round
View(bt64)
bt65 <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[65]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =247,  not to round
View(bt65)
bt66 <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[66]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =252,  in 2 decimals
View(bt66)
bt67 <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[67]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =253,  not to round
View(bt67)
bt68 <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[68]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =256,  to round, 1 decimals
View(bt68)
bt69 <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[69]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =271,  not to round
View(bt69)
bt70 <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[70]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =272,  not to round
View(bt70)
bt_ <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[71]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =273,  not to round
View(bt_)
bt_ <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[72]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =274,  not to round
View(bt_)
bt_ <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[73]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =278,  to round
View(bt_)
bt_ <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[74]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =285,  not to round?
View(bt_)
bt_ <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[75]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =286, to round, 1 decimals
View(bt_)
bt_ <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[76]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =288, to round, 1 decimals
View(bt_)
bt_ <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[77]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =289, not to round?
View(bt_)
bt_ <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[78]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =290, not to round?
View(bt_)
bt_ <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[79]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =291, not to round?
View(bt_)
bt_ <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[80]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =292, not to round?
View(bt_)
bt_ <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[81]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =297, not to round
View(bt_)
bt_ <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[82]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =300, not to round
View(bt_)
bt_ <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[83]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =304, not to round
View(bt_)
bt_ <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[84]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =309, not to round
View(bt_)
bt_ <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[85]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =321, not to round
View(bt_)
bt_ <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[86]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =330, not to round?
View(bt_)
bt_ <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[87]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =331, not to round
View(bt_)
bt_ <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[88]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =354, not to round
View(bt_)
bt_ <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[89]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =356, not to round
View(bt_)
bt_ <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[90]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =359, not to round
View(bt_)
bt_ <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[91]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =374, not to round
View(bt_)
bt_ <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[92]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =375, not to round
View(bt_)
bt_ <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[93]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =418, not to round
View(bt_)
bt_ <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[94]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =428, not to round
View(bt_)
bt_ <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[95]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =430, to round
View(bt_)
bt_ <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[96]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =431, to round
View(bt_)
bt_ <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[97]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =432, to round
View(bt_)
bt_ <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[98]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =433, not to round
View(bt_)
bt_ <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[99]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =466, to round
View(bt_)
bt_ <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[100]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =467, not to round
View(bt_)
bt_ <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[101]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =468, not to round
View(bt_)
bt_ <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[102]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =469, not to round
View(bt_)
bt_ <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[103]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =499, not to round
View(bt_)
bt_ <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[104]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =500, not to round
View(bt_)
bt_ <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[105]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =505, not to round
View(bt_)
bt_ <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[106]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =507, not to round
View(bt_)
bt_ <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[107]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =512, not to round
View(bt_)
bt_ <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[108]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =516, not to round
View(bt_)
bt_ <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[109]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =525, not to round
View(bt_)
bt_ <- bt_4loc_10yr %>% distinct() %>% filter(STUDY_ID == btstudy$STUDY_ID[110]) %>% arrange(LATITUDE, LONGITUDE, YEAR) %>% as.data.frame() # id =527, in 2 decimals
View(bt_)


# round coordinates to 2 decimals and define locations for these 37 datasets
studyid_round <- c(68, 71, 81, 86, 90, 119, 147, 148, 162, 166, 169, 171, 172, 173, 176, 178, 182, 183, 
                   187, 189, 190, 191, 197, 198, 204, 206, 208, 209, 213, 256, 278, 286, 288, 430, 431, 432, 466)

 # realm of these studies: 34 marine, 3 freshwater; of them, 23 marine and 1 freshwater were finally kept and used in analyses 
bt_meta %>% filter(STUDY_ID %in% studyid_round) %>% pull(REALM) %>% table()

bt_4loc_10yr <- bt_4loc_10yr %>% 
  mutate(LATITUDE = ifelse(STUDY_ID %in% studyid_round, 
                           round(LATITUDE, 2), LATITUDE),
         LONGITUDE = ifelse(STUDY_ID %in% studyid_round, 
                            round(LONGITUDE, 2), LONGITUDE)) %>%
  unite(col = location, LATITUDE, LONGITUDE, remove = FALSE) %>% 
  distinct()

bt <- bt %>% 
  filter(STUDY_ID %in% bt_4loc_10yr$STUDY_ID) %>%
  mutate(LATITUDE = ifelse(STUDY_ID %in% studyid_round, 
                           round(LATITUDE, 2), LATITUDE),
         LONGITUDE = ifelse(STUDY_ID %in% studyid_round, 
                            round(LONGITUDE, 2), LONGITUDE))
  

############
## the number and locations of sites are somewhat different across years
# filter the dataset to maximum the number of same samples for at least two years within a region

bt_4loc_10yr_filtered <- NULL
for(i in 1:length(unique(bt_4loc_10yr$STUDY_ID))){
  # perform loop for each study
  study <- bt_4loc_10yr %>% 
    filter(STUDY_ID == unique(bt_4loc_10yr$STUDY_ID)[i])
  
  # get a matrix with row are years and colunmns is samples
  year_sample <- as.matrix(xtabs(~ YEAR + location, data = unique(study[,c("YEAR", "location")]), sparse=TRUE))
  
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
    filter(location %in% colnames(year_sample) & YEAR %in% rownames(year_sample))
  
  bt_4loc_10yr_filtered <- bind_rows(bt_4loc_10yr_filtered, study_filtered)
}

ggplot(data = study_filtered, aes(LONGITUDE, LATITUDE)) + 
  facet_wrap(~YEAR) +
  geom_point(size = 1, alpha = 0.7) + 
  coord_fixed()

# remove years with < 4 locations and studies with number of years < 2 and duration < 10 years
bt_4loc_10yr_filtered <- bt_4loc_10yr_filtered %>% 
  group_by(STUDY_ID, YEAR) %>%
  mutate(n_loc = n_distinct(LATITUDE, LONGITUDE)) %>% 
  filter(n_loc >= 4) %>%
  group_by(STUDY_ID) %>%
  mutate(n_years = n_distinct(YEAR, na.rm = TRUE),
         duration = max(YEAR) - min(YEAR) + 1) %>%
  filter(n_years >= 2 & duration >= 10) %>%
  ungroup() %>%
  dplyr::select(-c(n_loc, n_years, duration))


# check how many studies left
bt_4loc_10yr_filtered %>% distinct(STUDY_ID) # 87 studies

# add the column "keep" to distinguish the locations that should be kept or removed
bt_4loc_10yr <- bt_4loc_10yr %>% 
  left_join(bt_4loc_10yr_filtered %>% mutate(keep = "yes")) %>%
  mutate(keep = ifelse(is.na(keep), "no", keep))

# plot distributions of locations and indicate which locations will be removed
pdf('data/BioTIME/BioTIME_4locations_10years_filtered.pdf', width = 12, height = 10)
id_study <- unique(bt_4loc_10yr$STUDY_ID)
for(i in 1:length(id_study)){
  study <- bt_4loc_10yr %>% filter(STUDY_ID %in% id_study[i])
  
  p <- ggplot(data = study, aes(LONGITUDE, LATITUDE)) + 
    facet_wrap(~YEAR) +
    geom_point( aes(colour = keep), size = 1, alpha = 0.7) + 
    labs(title = id_study[i]) +
    theme(legend.position = "top", legend.text = element_text(size=12)) + 
    coord_fixed() +
    scale_color_manual(values = c("yes" = "deepskyblue", "no" = "coral"))

  print(p)
}
dev.off()


# occurrence data for filtered locations
bt_filtered <- bt_4loc_10yr_filtered %>% 
  left_join(bt %>% dplyr::select(-X1), by = c("STUDY_ID",  "YEAR", "LATITUDE", "LONGITUDE")) %>%
  left_join(bt_meta, by = "STUDY_ID")

# Save data
save(bt_filtered, file = "data/BioTIME/BioTIME_4locations_10years_filtered.RDATA")


# check how many samples for each location
location_nsamples <- bt_filtered %>% 
  dplyr::select(STUDY_ID, YEAR, location, SAMPLE_DESC) %>% 
  distinct() %>%
  group_by(STUDY_ID, YEAR, location) %>%
  summarise(n_samp = n_distinct(SAMPLE_DESC)) %>%
  ungroup()

# ~18.5% of locations have more than one sample
table(location_nsamples$n_samp > 1)

# the minimum number of samples across locations within studies
study_location_nsamples <- location_nsamples %>% 
  group_by(STUDY_ID) %>% 
  summarise(n_location = n_distinct(location),
            min_samp = min(n_samp),
            max_samp = max(n_samp),
            median_samp = median(n_samp),
            mean_samp = mean(n_samp)) 

# 12 studies have multiple samples at all locations
table(study_location_nsamples$min_samp > 1)

# keep minimum number of samples across locations within studies for each location
set.seed(1)
bt_filtered_location <- bt_filtered %>% 
  dplyr::select(STUDY_ID, YEAR, location, SAMPLE_DESC) %>% 
  distinct() %>%
  group_by(STUDY_ID, YEAR, location) %>% 
  mutate(n_samp = n_distinct(SAMPLE_DESC)) %>%
  group_by(STUDY_ID) %>% 
  mutate(min_samp = min(n_samp)) %>%
  dplyr::select(-n_samp) %>%
  group_by(STUDY_ID, YEAR, location) %>% 
  sample_n(size = unique(min_samp))

bt_filtered <- bt_filtered %>% inner_join(bt_filtered_location %>% dplyr::select(-min_samp))

save(bt_filtered, file = "data/BioTIME/BioTIME_4locations_10years_filtered_final.RDATA")


## Check abundance and biomass information
table(bt_filtered$AB_BIO)
table(bt_filtered$ABUNDANCE_TYPE, useNA = "always")
table(bt_filtered$BIOMASS_TYPE, useNA = "always")

bt_filtered %>% filter(is.na(ABUNDANCE_TYPE)) %>% 
  distinct(STUDY_ID, AB_BIO, ABUNDANCE_TYPE)
# note: records with NA abundance type only have biomass estimates 

# check data set with ABUNDANCE_TYPE as "Presence/Absence"
x <- bt_filtered %>% filter(ABUNDANCE_TYPE == "Presence/Absence") %>% 
  select(STUDY_ID, AB_BIO, sum.allrawdata.ABUNDANCE, sum.allrawdata.BIOMASS, ABUNDANCE_TYPE) %>%
  group_by(STUDY_ID) %>%
  summarise(max_abundance = max(sum.allrawdata.ABUNDANCE))
print(x) # study 84, 110 and 173 have max abundance greater than  1. 

# check the above three study individually
y <- bt_filtered %>% filter(ABUNDANCE_TYPE == "Presence/Absence") %>% 
  select(STUDY_ID, AB_BIO, sum.allrawdata.ABUNDANCE, sum.allrawdata.BIOMASS, ABUNDANCE_TYPE) %>%
  filter(STUDY_ID == 173) # 84, 173, 110
view(y)
# note: most abundance information based on presence/absence are 1, so they are not really abundance 

# check data set with ABUNDANCE_TYPE as "Density"
x <- bt_filtered %>% filter(ABUNDANCE_TYPE == "Density") %>% 
  select(STUDY_ID, AB_BIO, sum.allrawdata.ABUNDANCE, sum.allrawdata.BIOMASS, ABUNDANCE_TYPE) %>%
  group_by(STUDY_ID) %>%
  summarise(min_abundance = min(sum.allrawdata.ABUNDANCE),
            mean_abundance = mean(sum.allrawdata.ABUNDANCE),
            max_abundance = max(sum.allrawdata.ABUNDANCE))
print(x) # note: some values are decimals (very small, close to zero)


# check data set with ABUNDANCE_TYPE as "Count"
x <- bt_filtered %>% filter(ABUNDANCE_TYPE == "Count") %>% 
  select(STUDY_ID, AB_BIO, sum.allrawdata.ABUNDANCE, sum.allrawdata.BIOMASS, ABUNDANCE_TYPE) %>%
  group_by(STUDY_ID) %>%
  summarise(min_abundance = min(sum.allrawdata.ABUNDANCE),
            mean_abundance = mean(sum.allrawdata.ABUNDANCE),
            max_abundance = max(sum.allrawdata.ABUNDANCE),
            range_abundance = max_abundance - min_abundance)
print(x, n = 66)
# note: even recorded as count, some abundance values are decimals and smaller than 1.

x %>% filter(max_abundance <= 10 | min_abundance < 1)
# note: all abundance values in study 67 are 1, so it is not really abundance 
# note: study 67 and 467 have min abundance valuaes as zero.

# check datasets with small min and max abundance values 
bt_filtered %>% filter(STUDY_ID == 176) %>% ggplot(aes(x =sum.allrawdata.ABUNDANCE)) + geom_histogram() + scale_x_log10()

bt_filtered %>% filter(STUDY_ID == 253) %>% ggplot(aes(x =sum.allrawdata.ABUNDANCE)) + geom_histogram() + scale_x_log10()

bt_filtered %>% filter(STUDY_ID == 304) %>% ggplot(aes(x =sum.allrawdata.ABUNDANCE)) + geom_histogram() 
# note: 5 abundance calss, 1: 1-4, 2:5-10, 3: 11-50, 4: 51-100, 5: >100; exclude this dataset in abundance analyses

bt_filtered %>% filter(STUDY_ID == 321) %>% ggplot(aes(x =sum.allrawdata.ABUNDANCE)) + geom_histogram() 
# note: most abundance values are one; exclude this dataset in abundance analyses

bt_filtered %>% filter(STUDY_ID == 433) %>% ggplot(aes(x =sum.allrawdata.ABUNDANCE)) + geom_histogram() + scale_x_log10()

bt_filtered %>% filter(STUDY_ID == 467) %>% ggplot(aes(x =sum.allrawdata.ABUNDANCE)) + geom_histogram() + scale_x_log10()

bt_filtered %>% filter(STUDY_ID == 176) %$% table(sum.allrawdata.ABUNDANCE == 0)
bt_filtered %>% filter(STUDY_ID == 176) %>% filter(sum.allrawdata.ABUNDANCE == 0) %$% table(sum.allrawdata.BIOMASS > 0) 
#note: 200 of 1214 abundance values are zero, but corresponding biomass are not zero; replace abundance zero with 1

bt_filtered %>% filter(STUDY_ID == 467) %$% table(sum.allrawdata.ABUNDANCE == 0)
bt_filtered %>% filter(STUDY_ID == 467) %>% filter(sum.allrawdata.ABUNDANCE == 0) %$% table(sum.allrawdata.BIOMASS > 0) 
#note: 6 of 819 abundance values are zero, but corresponding biomass are not zero; replace abundance zero with 1


# check data set with ABUNDANCE_TYPE as "Weight"
x <- bt_filtered %>% filter(BIOMASS_TYPE == "Weight") %>% 
  select(STUDY_ID, AB_BIO, sum.allrawdata.ABUNDANCE, sum.allrawdata.BIOMASS, ABUNDANCE_TYPE, BIOMASS_TYPE) %>%
  group_by(STUDY_ID, AB_BIO, ABUNDANCE_TYPE, BIOMASS_TYPE) %>%
  summarise(min_abundance = min(sum.allrawdata.ABUNDANCE),
            mean_abundance = mean(sum.allrawdata.ABUNDANCE),
            max_abundance = max(sum.allrawdata.ABUNDANCE),
            min_biomass = min(sum.allrawdata.BIOMASS),
            mean_biomass = mean(sum.allrawdata.BIOMASS),
            max_biomass = max(sum.allrawdata.BIOMASS))
print(x) # note: AB_BIO values is "AB" or "B". some biomass values are zero and very small values. 

# check data set with ABUNDANCE_TYPE as "Cover"
x <- bt_filtered %>% filter(BIOMASS_TYPE == "Cover") %>% 
  select(STUDY_ID, AB_BIO, sum.allrawdata.ABUNDANCE, sum.allrawdata.BIOMASS, ABUNDANCE_TYPE, BIOMASS_TYPE) %>%
  group_by(STUDY_ID, AB_BIO, ABUNDANCE_TYPE, BIOMASS_TYPE) %>%
  summarise(min_abundance = min(sum.allrawdata.ABUNDANCE),
            mean_abundance = mean(sum.allrawdata.ABUNDANCE),
            max_abundance = max(sum.allrawdata.ABUNDANCE),
            min_biomass = min(sum.allrawdata.BIOMASS),
            mean_biomass = mean(sum.allrawdata.BIOMASS),
            max_biomass = max(sum.allrawdata.BIOMASS))
print(x) # note: values from 1 to 91. All AB_BIO values is "B".


## Solution: only use abundance data measured with "count" and "density" (delete "Presence/Absence")
# for datasets 176 and 467, replace abundance zero with 1
# for datasets 67, 304, and 321, exclude them in abundance analyses
bt_filtered <- bt_filtered %>% 
  mutate(sum.allrawdata.ABUNDANCE = ifelse(sum.allrawdata.ABUNDANCE==0 & STUDY_ID == 176, 1, sum.allrawdata.ABUNDANCE)) %>% 
  mutate(sum.allrawdata.ABUNDANCE = ifelse(sum.allrawdata.ABUNDANCE==0 & STUDY_ID == 467, 1, sum.allrawdata.ABUNDANCE)) %>% 
  mutate(ABUNDANCE_TYPE = ifelse(STUDY_ID %in% c(67, 304, 321), "Count_delete", ABUNDANCE_TYPE)) 
  
save(bt_filtered, file = "data/BioTIME/BioTIME_4locations_10years_filtered_final.RDATA")


# the StudyID used
bt_studyID_filtered <- bt_filtered %>% distinct(STUDY_ID, HABITAT, TAXA, ORGANISMS, TITLE )
write.csv(bt_studyID_filtered, file = "data/BioTIME/bt_studyID_filtered.csv", row.names = FALSE)
