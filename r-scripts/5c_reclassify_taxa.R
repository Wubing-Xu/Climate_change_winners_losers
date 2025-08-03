##################################
## reclassify the taxonomic groups for those that are listed as "multiple taxa" and "Benthos".
# check the taxonomic rank of species in my filtered data using the GBIF backbone
# the taxonomic group with most species within a study will be used  

rm(list = ls())

# Set user dependent working directories
user <- Sys.info()["nodename"]
path2wd <- switch(user,
                  "IDIVNB341" = "C:/Dropbox/iDiv/Populations_warming/analyses_202406",
                  "IDIVTS01" = "H:/wubing/iDiv/Populations_warming/analyses_202406")
setwd(path2wd)

# load packages
needed_libs <- c("tidyverse","vegan", "reshape2", "sf")

usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1])) {   
    install.packages(p, dep = TRUE)
  }
  require(p, character.only = TRUE)
}

sapply(needed_libs, usePackage)
rm(usePackage)


load("data/Combined_assemblages.RDATA")
load("data/combined_checklists/splist_gbif_all.columns.RDATA")


#####################
## the studies are listed as "multiple taxa".
meta_mtaxa <-  dat_meta %>% 
    filter(taxon_new == "Multiple taxa")
meta_mtaxa$study # bt_169" "bt_172" "bt_182" "bt_274" "bt_505"

species_dat_mtaxa <- dat %>%
  filter(study %in% meta_mtaxa$study) %>%
  distinct(study, species, specieskey) %>%
  left_join(splist.gbif.all %>% 
              dplyr::select(specieskey, kingdom:family) %>%
              distinct())

# study bt_169: 7 mammals and 5 birds. some migratory birds have been removed
species_dat_mtaxa %>% filter(study == "bt_169") %>% 
  count(kingdom)

species_dat_mtaxa %>% filter(study == "bt_169") %>% 
  count(order)


# study bt_172: mammals
species_dat_mtaxa %>% filter(study == "bt_172") %>% 
  count(kingdom)

species_dat_mtaxa %>% filter(study == "bt_172") %>% 
  count(order)
# order: 13 mammals (13 Cetacea), 
# 1 birds (1 Procellariiformes); 
# 2 turtles (Testudines)


# study bt_182: fishes
species_dat_mtaxa %>% filter(study == "bt_182") %>% 
  count(kingdom)

species_dat_mtaxa %>% filter(study == "bt_182") %>% 
  count(phylum)

species_dat_mtaxa %>% filter(study == "bt_182") %>% 
  count(order)
# order: fishes: 16 fishes (5 Pleuronectiformes,5 Scorpaeniformes, 3 Perciformes, 2 Gadiformes, Rajiformes)
# invertebrates: 11 invertebrates (3 Decapoda, Amphilepidida, Camarodonta, Dendrochirotida, Euryalida, Mytilida, Myxiniformes, Pectinida, Phyllodocida)
# Chordariales: 1;  1 Echinolampadacea; 


# study bt_274: invertebrates
species_dat_mtaxa %>% filter(study == "bt_274") %>% 
  count(kingdom)

species_dat_mtaxa %>% filter(study == "bt_274") %>% 
  count(phylum)
# phylum: invertebrate: 37 species (10 Cnidaria, 7 Bryozoa, 6 Annelida,  2 Chordata, 4 Echinodermata, 4 Mollusca, 4 Porifera)
# plants (alga): 29 species (17 Rhodophyta, 9 Ochrophyta, 2 Chlorophyta, 2 Tracheophyta)


# study bt_505: most are fishes
species_dat_mtaxa %>% filter(study == "bt_505") %>% 
  count(kingdom)

species_dat_mtaxa %>% filter(study == "bt_505") %>% 
  count(phylum)

species_dat_mtaxa %>% filter(study == "bt_505") %>% 
  count(order) %>% as.data.frame()
# 89 Chordata (most are fishes, 48 Perciformes,5 Scorpaeniformes, 5 Clupeiformes, 5 Pleuronectiformes, 5 Tetraodontiformes), 
# 9 Decapoda, 8 Mollusca, 1 Echinodermata, 1 Cnidaria


#########################
## the studies are listed as "Benthos".

meta_benthos <-  dat_meta %>% 
  filter(taxon_new == "Benthos")
meta_benthos$study # "bt_78"  "bt_110" "bt_162" "bt_163" "bt_187" "bt_196" "bt_204" "bt_213" "bt_468" "bt_469"


species_dat_benthos <- dat %>%
  filter(study %in% meta_benthos$study) %>%
  distinct(study, species, specieskey) %>%
  left_join(splist.gbif.all %>% 
              dplyr::select(specieskey, kingdom:family) %>%
              distinct())


# study bt_78: most are invertebrates
species_dat_benthos %>% filter(study == "bt_78") %>% 
  count(kingdom)

species_dat_benthos %>% filter(study == "bt_78") %>% 
  count(phylum)
# 58 Annelida, 31 Arthropoda, 25 Mollusca, 10 Cnidaria, 8 Bryozoa


# study bt_110: most are invertebrates
species_dat_benthos %>% filter(study == "bt_110") %>% 
  count(kingdom)

species_dat_benthos %>% filter(study == "bt_110") %>% 
  count(phylum)
# invertebrates: 139 Mollusca, 77 Annelida, 53 Arthropoda, 35 Cnidaria, 34 Chordata, 29 Echinodermata, 15 Bryozoa, 11 Porifera
# Plantae: 64 Rhodophyta, 27 Ochrophyta, 
# Fungi 10: Ascomycota


# study bt_162: all are invertebrates
species_dat_benthos %>% filter(study == "bt_162") %>% 
  count(kingdom)

species_dat_benthos %>% filter(study == "bt_162") %>% 
  count(phylum)
# 43 Annelida, 29 Mollusca, 12 Arthropoda, 2 Echinodermata


# study bt_163: most are fish
species_dat_benthos %>% filter(study == "bt_163") %>% 
  count(kingdom)

species_dat_benthos %>% filter(study == "bt_163") %>% 
  count(phylum)

species_dat_benthos %>% filter(study == "bt_163") %>% 
  count(order) %>% as.data.frame()
# fish: 54 Scorpaeniformes, 22 Pleuronectiformes, 15 Perciformes, 12 Rajiformes
# invertebrates: 19 Arthropoda, 2 Mollusca, Annelida, Cnidaria


# study bt_187: all are invertebrates
species_dat_benthos %>% filter(study == "bt_187") %>% 
  count(kingdom)

species_dat_benthos %>% filter(study == "bt_187") %>% 
  count(phylum)
# 17 Annelida, 12 Mollusca, 6 Arthropoda


# study bt_196: about half are invertebrates
species_dat_benthos %>% filter(study == "bt_196") %>% 
  count(kingdom)

species_dat_benthos %>% filter(study == "bt_196") %>% 
  count(phylum)
# invertebrates: 21 Mollusca, 9 Arthropoda, 6 Porifera, 5 Chordata, 4 Bryozoa, 4 Cnidaria, 4 Echinodermata
# Fungi: 11
# Plantae: 35 + 19 Chromista


# study bt_204: all are invertebrates
species_dat_benthos %>% filter(study == "bt_204") %>% 
  count(kingdom)

species_dat_benthos %>% filter(study == "bt_204") %>% 
  count(phylum)
# 32 Annelida, 24 Arthropoda, 19 Mollusca, 4 Echinodermata


# study bt_213: most are fishes
species_dat_benthos %>% filter(study == "bt_213") %>% 
  count(kingdom)

species_dat_benthos %>% filter(study == "bt_213") %>% 
  count(phylum)

species_dat_benthos %>% filter(study == "bt_213") %>% 
  count(order)  %>% as.data.frame()
# fish: 27 Perciformes, 12 Pleuronectiformes, 11 Gadiformes, 10 Scorpaeniformes, 7 Clupeiformes, 4 Rajiformes
# invertebrates: 14 Arthropoda, 5 Mollusca


# study bt_468: about half are plants
species_dat_benthos %>% filter(study == "bt_468") %>% 
  count(kingdom)

species_dat_benthos %>% filter(study == "bt_468") %>% 
  count(phylum)


# study bt_469: most are invertebrates
species_dat_benthos %>% filter(study == "bt_469") %>% 
  count(kingdom)

species_dat_benthos %>% filter(study == "bt_469") %>% 
  count(phylum)
# invertebrates: 12 Mollusca, 8 Echinodermata, 4 Cnidaria, 4 Chordata, Porifera, Arthropoda
# Chromista: 6


#############
# to update taxonomic groups
dat_meta_mtaxa_manul <- data.frame(study= c("bt_169", "bt_172", "bt_182", "bt_274", "bt_505"),
                                  taxon_new = c("Mammals", "Mammals", "Fish", "Invertebrates", "Fish"))

dat_meta_benthos_manual <- data.frame(study= c("bt_78","bt_110", "bt_162", "bt_163",  
                                             "bt_196", "bt_204", "bt_213", "bt_468", "bt_469"),
                                  taxon_new = c("Invertebrates", "Invertebrates", "Invertebrates", "Fish", 
                                                "Invertebrates", "Invertebrates", "Fish", "Plants", "Invertebrates"))

# taxon information
dat_meta_taxa <- dat_meta %>%
  dplyr::select(study, database, studyID, study_name, realm, taxon, taxon_new)

# update the "multiple taxa"
id <- match(dat_meta_mtaxa_manul[, "study"], dat_meta_taxa$study)
dat_meta_taxa[id, c("taxon_new")] <- dat_meta_mtaxa_manul[, c("taxon_new")]

# update the "benthos"
id <- match(dat_meta_benthos_manual[, "study"], dat_meta_taxa$study)
dat_meta_taxa[id, c("taxon_new")] <- dat_meta_benthos_manual[, c("taxon_new")]


# distinguish plants, invertebrates and fish in different realms (terrestrial, freshwater and marine)
dat_meta_taxa <- dat_meta_taxa %>% 
  mutate(taxon_final = ifelse(taxon_new %in% c("Fish", "Invertebrates", "Plants", "Mammals"), 
                              paste(realm, tolower(taxon_new), sep = " "), taxon_new))
table(dat_meta_taxa$taxon_final)

save(dat_meta_taxa, file = "data/Assemblages_taxa.RDATA")
