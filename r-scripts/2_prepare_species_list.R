## Combine species from four assemblage database, add and unify realm, kingdom, and match to GBIF taxonomy backbone
# the matched GBIF species names were used to extract GBIF occurrences

rm(list = ls())
setwd("C:/Dropbox/iDiv/Populations_warming/analyses_202406")

library(tidyverse)
library(rgbif)

####################
#### combine all species from four database


#### BioTIME species list
load("data/BioTIME/BioTIME_4locations_10years_filtered_final.RDATA")

# species list in BioTIME
bt_species <- bt_filtered %>% 
  distinct(STUDY_ID, GENUS_SPECIES, REALM, HABITAT, TAXA, ORGANISMS)

# add kingdom according to descriptions in TAXA and ORGANISMS
id_plantae <- bt_species$TAXA %in% c("Terrestrial plants", "Marine plants", "Freshwater plants") &
  !bt_species$ORGANISMS %in% c("Macroalgae", "fungi", "moss", "lichen", "tropical algae")
table(id_plantae)

id_Animalia <- (bt_species$TAXA %in% c("Birds", "Fish", " Mammals", "Terrestrial invertebrates", "Marine invertebrates", 
                                       "Freshwater invertebrates", "Reptiles", "Amphibians") &
                  !bt_species$ORGANISMS %in% c("intertidal", "marine invertebrates & marine plants", "Microplankton", "Phytoplankton",
                                               "Phytoplankton + zooplankton (and some plants - changed taxa from /plants)", "turf algae")) |
  (bt_species$TAXA %in% c("Benthos") & bt_species$ORGANISMS %in% c("benthic animals", "Echinodermata",  "Epibenthic megafauna", 
                                                                   "macrozoobenthos", "Polychaeta","zoobenthos")) |
  (bt_species$TAXA %in% c("All") & bt_species$ORGANISMS %in% c("All animals observed around waterholes", "birds and some marine mammals", "Birds carnivores and ungulates",                                                             
                                                               "Cetaceans. seabirds + turtles", "fish and marine invertebrates", "Fish and marine invertebrates",  
                                                               "Herbivores carnivores and eagles", "Herpetofauna", "mostly seabirds + some marine mammals", "pelagic seabirds"))
table(id_Animalia)

bt_species[id_plantae, "kingdom"] <- "Plantae"
bt_species[id_Animalia, "kingdom"] <- "Animalia"

# the BioTIME species list
btsplist <- bt_species %>% 
  distinct(GENUS_SPECIES, REALM, TAXA, kingdom) %>%
  setNames(tolower(names(.))) %>% 
  rename(species = genus_species, taxon = taxa) %>% 
  mutate(realm = tolower(realm))


#### checklist of metacommunity resurvey database
load("data/Metacommunity_Resurvey/metacommunityResurvey_filtered.RDATA")

# species list
mr_species <- mr_filtered %>% distinct(dataset_id, species, realm, taxon)

# check whether all species have realm and taxon information
mr_species %>% filter(is.na(realm))

# add kingdom to species list
unique(mr_species$taxon)
taxon_kingdom <- data.frame(taxon = c("Fish", "Invertebrates", "Marine plants", "Herpetofauna",  "Birds", "Plants", "Mammals", "Phytoplankton"),
                            kingdom = c("Animalia", "Animalia", "Plantae", "Animalia", "Animalia", "Plantae", "Animalia", "Chromista"))

# check whether all names can be matched
unique(mr_species$taxon) %in% taxon_kingdom$taxon

mr_species <- mr_species %>% left_join(taxon_kingdom)

# the prepared species list for metacommunity resurvey database
mrsplist <- mr_species %>% 
  dplyr::select(species, realm, taxon, kingdom) %>% 
  distinct()


#### RivFishTIME species checklist
load("data/RivFishTIME/RivFishTIME_4sites_10years_filtered.RDATA")

ftsplist <- ft_filtered %>% distinct(species = Species) %>%
  mutate(realm = "freshwater", taxon = "fish", kingdom = "Animalia")


#### species checklist of insect dataset
load("data/InsectChange/InsectChange_filtered.RDATA")

icsplist <- ic_filtered %>% 
  distinct(species, realm = Realm) %>%
  mutate(taxon = "invertebrate", kingdom = "Animalia")


##combine checklist
splist <- bind_rows(mrsplist, btsplist, ftsplist, icsplist) %>% distinct()

# update "realm" for those species with inconsistent information of the same species 
splist <- splist %>%
  mutate(realm = tolower(realm))

check.realm <- splist %>% 
  distinct(species, realm) %>% 
  filter(species %in% pull(.[duplicated(species),], species)) %>%
  mutate(realm.new = NA)

# combine realms for the same species
check.realm.splist <- distinct(check.realm, species) %>% unlist()
for(i in 1:length(check.realm.splist)){
  id <- check.realm[,1] == check.realm.splist[i]
  check.realm$realm.new[id]<- paste(sort(unlist(check.realm[id, 2])), collapse="_")
}
check.realm <- check.realm %>% distinct(species, realm.new)

# use the updated to replace the raw
splist <- splist %>% left_join(check.realm) %>% 
  mutate(realm = ifelse(is.na(realm.new), realm, realm.new)) %>% 
  dplyr::select(!realm.new) %>% distinct()


# update "kingdom" for those species with inconsistent information of the same species 
check.kingdom <- splist %>% 
  distinct(species, kingdom) %>% 
  filter(species %in% pull(.[duplicated(species),], species)) %>%
  mutate(kingdom.new = NA)

check.kingdom.splist <- distinct(check.kingdom, species) %>% unlist()
for(i in 1:length(check.kingdom.splist)){
  id <- check.kingdom[,1] == check.kingdom.splist[i]
  kingdom <- unlist(check.kingdom[id, 2])
  kingdom <- kingdom[!is.na(kingdom)]
  check.kingdom$kingdom.new[id] <- ifelse(length(kingdom) ==1, kingdom, NA)
}

check.kingdom <- check.kingdom %>% 
  distinct(species, kingdom.new) %>%
  mutate(kingdom.new = ifelse(is.na(kingdom.new), "unknown", kingdom.new))

# use the updated to replace the raw
splist <- splist %>% left_join(check.kingdom) %>% 
  mutate(kingdom = ifelse(is.na(kingdom.new), kingdom, kingdom.new)) %>% 
  mutate(kingdom = ifelse(kingdom == "unknown", NA, kingdom)) %>% 
  dplyr::select(!kingdom.new) %>% distinct()


# the final splist
splist <- distinct(splist, species, realm, kingdom)



#########
# find names in gbif backbone
spgbif <- lapply(1:nrow(splist),function(x) name_backbone(name=splist[x,1], kingdom=splist[x, 3]))
spgbif <- bind_rows(spgbif)

# Combine species checklist with information from gbif backbone
splist.gbif.all <- splist %>% 
  rename(input.species = species, input.kingdom = kingdom) %>% 
  bind_cols(spgbif %>% setNames(tolower(names(.))))

# the simplified species list with names in gbif backbone
splist.gbif <- splist.gbif.all %>% select(c(input.species:input.kingdom, species, specieskey)) # 15469 species

splist.gbif.specieskey <- splist.gbif %>% filter(!is.na(specieskey)) %>% distinct(specieskey) %>% pull()
length(splist.gbif.specieskey) # 14,987 species keys


# checklist with only names from GBIF. remove duplicated
spgbif <- splist.gbif %>% 
  filter(!is.na(specieskey)) %>% 
  dplyr::select(species, specieskey, realm) %>% 
  distinct() 

## update "realm" for those species with inconsistent information of the same species 
check.realm <- spgbif %>% 
  distinct(specieskey, realm) %>% 
  filter(specieskey %in% pull(.[duplicated(specieskey),], specieskey)) %>%
  mutate(realm.new = NA)

# combine realms for the same species
check.realm.splist <- distinct(check.realm, specieskey) %>% unlist()
for(i in 1:length(check.realm.splist)){
  id <- check.realm[ ,1] == check.realm.splist[i]
  realm.new <- paste(sort(unique(unlist(strsplit(unlist(check.realm[id, 2]), "_")))),collapse="_")
  check.realm$realm.new[id]<- realm.new
}
check.realm <- check.realm %>% distinct(specieskey, realm.new)

# use the updated to replace the raw
spgbif <- spgbif %>% 
  left_join(check.realm) %>% 
  mutate(realm = ifelse(is.na(realm.new), realm, realm.new)) %>% 
  dplyr::select(!realm.new) %>% distinct()

# determine species' distribution habitat type: land or ocean. It used to filter GBIF occurrences
spgbif <- spgbif %>% 
  mutate(keep = ifelse(realm %in% c("terrestrial", "freshwater", "freshwater_terrestrial"), "land", NA)) %>%
  mutate(keep = ifelse(realm == "marine", "ocean", keep)) 

table(spgbif$realm, useNA = "always")
table(spgbif$keep, useNA = "always")

# output
dir.create("data/combined_checklists/")
save(splist.gbif.all, file="data/combined_checklists/splist_gbif_all.columns.RDATA")
save(splist.gbif, splist.gbif.specieskey, spgbif, file="data/combined_checklists/splist_gbif.RDATA")

splist.gbif.all <- splist.gbif.all %>% select(-c(verbatim_name, verbatim_kingdom))
write_csv(splist.gbif.all, file="data/combined_checklists/splist_gbif_all.columns.csv") #save as a csv to inspect

