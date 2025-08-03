#### create all supplementary plots through this r-script

rm(list = ls())
gc()

# load packages
needed_libs <- c("tidyverse","ggplot2", "brms", "tidybayes", "abind", "cowplot", "ggridges")

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


######################
## figures showing regional variation in the relationships between occupancy/abundance change and baseline thermal position
rm(list = ls())
gc()

load("results/data_input_to_models.RDATA") # observation data
load("results/coef_oc_ac_tmeanpos.RDATA") # the regional-level slopes of baseline thermal position


# maps show regional variation in the slopes of baseline thermal position on occupancy and abundance changes
## draw figure
# world map 
worldmap <- ne_countries(scale = 'medium', type = 'map_units',  returnclass = 'sf')
worldmap <- worldmap %>% filter(worldmap$admin != "Antarctica")

# generate a spatial point file containing temperature changes
coef_oc_tmeanpos_points <- coef_oc_tmeanpos %>% 
  dplyr::select(cent_long, cent_lat, realm, change_tempmean, estimate_slope) %>% 
  mutate(change_tempmean = ifelse(change_tempmean <= -0.1, -0.1, ifelse(change_tempmean >= 0.1, 0.1, change_tempmean)),
         estimate_slope = ifelse(estimate_slope <= -0.1, -0.1, ifelse(estimate_slope >= 0.1, 0.1, estimate_slope))) %>%
  st_as_sf(coords = c('cent_long', 'cent_lat'))  %>% 
  st_set_crs(4326) %>%
  st_transform(crs = "+proj=moll")

coef_ac_tmeanpos_points <- coef_ac_tmeanpos %>% 
  dplyr::select(cent_long, cent_lat, realm, change_tempmean, estimate_slope) %>% 
  mutate(estimate_slope = ifelse(estimate_slope <= -0.07, -0.07, ifelse(estimate_slope >= 0.07, 0.07, estimate_slope))) %>%
  st_as_sf(coords = c('cent_long', 'cent_lat'))  %>% 
  st_set_crs(4326) %>%
  st_transform(crs = "+proj=moll")

# bounding of north america
namerica <- st_sfc(st_point(c(-115, 21)), st_point(c(-65, 55)), crs = 4326) %>%
  st_transform(crs = "+proj=moll") %>% 
  st_bbox() %>%
  st_as_sfc()

# bounding of Europe
europe <- st_sfc(st_point(c(-17, 38)), st_point(c(30, 60)), crs = 4326) %>%
  st_transform(crs = "+proj=moll") %>% 
  st_bbox() %>%
  st_as_sfc()

# generate the legend of realm
plot_legend_realm <- ggplot() + 
  geom_point(data = coef_oc_tmeanpos, aes(x = cent_long, y = cent_lat, shape = realm), fill = "gray50") + 
  scale_shape_manual(values = c("Terrestrial" = 22, "Freshwater" = 21, "Marine" = 24)) + 
  labs(shape = NULL) +
  theme_bw() +
  theme(legend.direction = "horizontal",
        legend.key.size = unit(5, 'mm'),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10, hjust = 0.2)) + 
  guides(fill = "none") + 
  guides(shape = guide_legend(nrow = 1))

legend_realm <- get_legend(plot_legend_realm)


###### map of occupancy change
map_occupchange_global <- ggplot() + 
  geom_sf(data = worldmap, fill = "gray95", colour = "gray50", size = 0.2) + 
  geom_sf(data = namerica, fill = NA, colour = "black") + 
  geom_sf(data = europe, fill = NA, colour = "black") + 
  geom_sf(data = coef_oc_tmeanpos_points, aes(fill = estimate_slope, shape = realm), color = "black",
          size = 1.5, alpha = 0.8) + 
  coord_sf(crs = st_crs('+proj=moll'), ylim = c(-8172663, 9020048), expand = FALSE) + 
  scale_fill_gradient2(low = "blue", mid = "darkgray", high = "red") + #breaks = c(-2, -1, 0, 1, 2), 
  scale_shape_manual(values = c("Terrestrial" = 22, "Freshwater" = 21, "Marine" = 24)) + 
  #  annotate("text", x =1.125*-11046300, y= 0.80*6386580, label = "North\nAmerica", size = 2.46) +
  #  annotate("text", x =1.75*-1469518, y= 0.90*6876759, label = "Europe", size = 2.46) +
  theme_bw() +
  theme(plot.margin = unit(c(0,-0.5,0,-1), "cm"),
        legend.position = c(0.55, 0.10), 
        legend.direction = 'horizontal',
        legend.key.size = unit(3, 'mm'),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 7),
        legend.background = element_rect(fill = 'transparent'),
        panel.border = element_blank(),
        title = element_text(size = 8),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  guides(shape = FALSE) + 
  guides(fill = guide_colourbar(title = "Slope of baseline thermal position", title.position = "top", barwidth = 5))

histgram_occupchange <- ggplot(data = coef_oc_tmeanpos) +
  geom_histogram(aes(estimate_slope), fill = "gray60") +
  geom_vline(xintercept = 0, lty = 2) + 
  geom_vline(data = fixef_oc_tmeanpos, 
             aes(xintercept = Estimate), linetype =1) + 
  geom_rect(data = fixef_oc_tmeanpos, 
            aes(xmin = Q2.5, xmax = Q97.5, ymin = -Inf, ymax = Inf), alpha = 0.6, fill = "gray30") + 
  labs(x= "Slope", y = "Number of regions" )+ 
  theme_classic() +  
  theme(axis.title = element_text(size = 7),
        axis.text = element_text(size = 6),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank())

map_occupchange_global <- ggdraw(map_occupchange_global) +
  draw_plot(histgram_occupchange, x = 0.01, y = 0.01, width = 0.3, height = 0.5)

# map of north america
map_occupchange_namerica <- ggplot() + 
  geom_sf(data = worldmap, fill = "gray95", colour = "gray50", size = 0.2) + 
  geom_sf(data = coef_oc_tmeanpos_points, aes(fill = estimate_slope, shape = realm), color = "black",
          size = 2, alpha = 0.8) + 
  coord_sf(crs = st_crs('+proj=moll'), expand = FALSE, 
           xlim = st_coordinates(namerica)[1:2, "X"], ylim = st_coordinates(namerica)[2:3, "Y"]) + 
  annotate("text", x =0.945*-11046300, y= 0.92*6386580, label = "North\nAmerica", size = 2.46) +
  scale_fill_gradient2(low = "blue", mid = "gray", high = "red") +
  scale_shape_manual(values = c("Terrestrial" = 22, "Freshwater" = 21, "Marine" = 24)) +
  theme_bw() +
  theme(plot.margin = unit(c(1,1,0,0), "mm"),
        legend.position = "n", 
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# map of Europe
map_occupchange_europe <- ggplot() + 
  geom_sf(data = worldmap, fill = "gray95", colour = "gray50", size = 0.2) + 
  geom_sf(data = coef_oc_tmeanpos_points, aes(fill = estimate_slope, shape = realm), color = "black",
          size = 2, alpha = 0.8) + 
  coord_sf(crs = st_crs('+proj=moll'), expand = TRUE, 
           xlim = st_coordinates(europe)[1:2, "X"], ylim = st_coordinates(europe)[2:3, "Y"]) + 
  annotate("text", x =0.85*-1469518, y= 0.99*6876759, label = "Europe", size = 2.46) +
  scale_fill_gradient2(low = "blue", mid = "gray", high = "red") +
  scale_shape_manual(values = c("Terrestrial" = 22, "Freshwater" = 21, "Marine" = 24)) +
  theme_bw() +
  theme(plot.margin = unit(c(1,1,0,0), "mm"),
        legend.position = "n", 
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

right <- cowplot::plot_grid(map_occupchange_namerica, map_occupchange_europe, ncol =1, nrow = 2, rel_heights = c(1, 1.135))
map_occupchange <- cowplot::plot_grid(map_occupchange_global, right,  ncol =2, nrow = 1, rel_widths = c(2.52, 1),
                                      labels = "a Occupancy change", label_size = 10, hjust = 0, vjust = 1) +
  theme(plot.background = element_rect(fill = "white", colour = NA))

# add the legend of realm
map_occupchange <- cowplot::plot_grid(legend_realm, map_occupchange, ncol =1, nrow = 2, rel_heights = c(0.1, 1)) +
  theme(plot.background = element_rect(fill = "white", colour = NA))


###### map of abundance change
map_Nchange_global <- ggplot() + 
  geom_sf(data = worldmap, fill = "gray95", colour = "gray50", size = 0.2) + 
  geom_sf(data = namerica, fill = NA, colour = "black") + 
  geom_sf(data = europe, fill = NA, colour = "black") + 
  geom_sf(data = coef_ac_tmeanpos_points, aes(fill = estimate_slope, shape = realm), color = "black",
          size = 1.5, alpha = 0.8) + 
  coord_sf(crs = st_crs('+proj=moll'), ylim = c(-8172663, 9020048), expand = FALSE) + 
  scale_fill_gradient2(low = "blue", mid = "darkgray", high = "red") + #breaks = c(-2, -1, 0, 1, 2), 
  scale_shape_manual(values = c("Terrestrial" = 22, "Freshwater" = 21, "Marine" = 24)) + 
  #  annotate("text", x =1.125*-11046300, y= 0.80*6386580, label = "North\nAmerica", size = 2.46) +
  #  annotate("text", x =1.75*-1469518, y= 0.90*6876759, label = "Europe", size = 2.46) +
  theme_bw() +
  theme(plot.margin = unit(c(0,-0.5,0,-1), "cm"),
        legend.position = c(0.55, 0.10), 
        legend.direction = 'horizontal',
        legend.key.size = unit(3, 'mm'),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 7),
        legend.background = element_rect(fill = 'transparent'),
        panel.border = element_blank(),
        title = element_text(size = 8),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  guides(shape = FALSE) + 
  guides(fill = guide_colourbar(title = "Slope of baseline thermal position", title.position = "top", hjust =0.5, barwidth = 5))

histgram_Nchange <- ggplot(data = coef_ac_tmeanpos) +
  geom_histogram(aes(estimate_slope), fill = "gray60") +
  geom_vline(xintercept = 0, lty = 2) + 
  geom_vline(data = fixef_ac_tmeanpos, 
             aes(xintercept = Estimate), linetype =1) + 
  geom_rect(data = fixef_ac_tmeanpos, 
            aes(xmin = Q2.5, xmax = Q97.5, ymin = -Inf, ymax = Inf), alpha = 0.6, fill = "gray30") + 
  labs(x= "Slope", y = "Number of regions" )+ 
  theme_classic() +  
  theme(axis.title = element_text(size = 7),
        axis.text = element_text(size = 6),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank())

map_Nchange_global <- ggdraw(map_Nchange_global) +
  draw_plot(histgram_Nchange, x = 0.01, y = 0.01, width = 0.3, height = 0.5)

# map of north america
map_Nchange_namerica <- ggplot() + 
  geom_sf(data = worldmap, fill = "gray95", colour = "gray50", size = 0.2) + 
  geom_sf(data = coef_ac_tmeanpos_points, aes(fill = estimate_slope, shape = realm), color = "black",
          size = 2, alpha = 0.8) + 
  coord_sf(crs = st_crs('+proj=moll'), expand = FALSE, 
           xlim = st_coordinates(namerica)[1:2, "X"], ylim = st_coordinates(namerica)[2:3, "Y"]) + 
  #  annotate("text", x =0.945*-11046300, y= 0.92*6386580, label = "North\nAmerica", size = 2.46) +
  scale_fill_gradient2(low = "blue", mid = "gray", high = "red") +
  scale_shape_manual(values = c("Terrestrial" = 22, "Freshwater" = 21, "Marine" = 24)) +
  theme_bw() +
  theme(plot.margin = unit(c(1,1,0,0), "mm"),
        legend.position = "n", 
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# map of Europe
map_Nchange_europe <- ggplot() + 
  geom_sf(data = worldmap, fill = "gray95", colour = "gray50", size = 0.2) + 
  geom_sf(data = coef_ac_tmeanpos_points, aes(fill = estimate_slope, shape = realm), color = "black",
          size = 2, alpha = 0.8) + 
  coord_sf(crs = st_crs('+proj=moll'), expand = TRUE, 
           xlim = st_coordinates(europe)[1:2, "X"], ylim = st_coordinates(europe)[2:3, "Y"]) + 
  #  annotate("text", x =0.85*-1469518, y= 0.99*6876759, label = "Europe", size = 2.46) +
  scale_fill_gradient2(low = "blue", mid = "gray", high = "red") +
  scale_shape_manual(values = c("Terrestrial" = 22, "Freshwater" = 21, "Marine" = 24)) +
  theme_bw() +
  theme(plot.margin = unit(c(1,1,0,0), "mm"),
        legend.position = "n", 
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

right <- cowplot::plot_grid(map_Nchange_namerica, map_Nchange_europe, ncol =1, nrow = 2, rel_heights = c(1, 1.135))
map_Nchange <- cowplot::plot_grid(map_Nchange_global, right,  ncol =2, nrow = 1, rel_widths = c(2.52, 1),
                                  labels = "b Abundance change", label_size = 10, hjust = 0, vjust = 1) +
  theme(plot.background = element_rect(fill = "white", colour = NA))


## combine maps of occupancy change and abundance change
maps_population_change <- cowplot::plot_grid(map_occupchange, map_Nchange, 
                                             ncol = 1, nrow = 2, rel_heights = c(1.1, 1)) +
  theme(plot.background = element_rect(fill = "white", colour = NA))

ggsave(maps_population_change, file = "results/supplement_Fig.S3_maps_of_regional_levelp_slopes.png", 
       width = 180, height = 139, units = 'mm', dpi = 600) # 1:0.7722
# ggsave(maps_population_change, file = "results/supplement_Fig.2_maps_of_regional_levelp_slopes.pdf", width = 180, height = 139, units = 'mm')


## plots showing regional-level regression lines of relationships between occupancy/abundance change and baseline thermal position
# show colour as the realm
realm_col <- c("Terrestrial" = "#D95F02", "Freshwater" = "#44bbff", "Marine" = "#7570B3")
realm_col <- c("All" = "#000000", "Terrestrial" = "#A65628", "Freshwater" = "#33A02C", "Marine" = "#1F78B4")

plot_oc_tp <- ggplot() +
  facet_wrap(~ realm) + 
  geom_segment(data = coef_oc_tmeanpos, 
               aes(x = xmin, xend = xmax, 
                   y = estimate_intercept + estimate_slope * xmin, 
                   yend = estimate_intercept + estimate_slope * xmax,
                   group = study, colour = realm), linewidth = 0.3, alpha = 0.3) + 
  geom_hline(yintercept = 0, lty = 2, linewidth = 0.3) + 
  labs(tag = "c", x = "Baseline thermal position", y = "Occupancy change") + 
  scale_colour_manual(name = NULL, values = realm_col) + 
  theme_bw()  +
  theme(legend.position = "no", 
        plot.margin = unit(c(0.2,0.2,0,0.2), "cm"),
        text = element_text(size =10),
        strip.text = element_text(size = 9),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 10), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

plot_ac_tp <- ggplot() +
  facet_wrap(~ realm) + 
  geom_segment(data = coef_ac_tmeanpos, 
               aes(x = xmin, xend = xmax, 
                   y = estimate_intercept + estimate_slope * xmin, 
                   yend = estimate_intercept + estimate_slope * xmax,
                   group = study, colour = realm), linewidth = 0.3, alpha = 0.3) + 
  geom_hline(yintercept = 0, lty = 2, linewidth = 0.3) + 
  labs(tag = "d", x = "Baseline thermal position", y = "Abundance change") + 
  scale_colour_manual(name = NULL, values = realm_col) + 
  theme_bw()  +
  theme(legend.position = "no", 
        plot.margin = unit(c(0,0.2,0.2,0.2), "cm"),
        text = element_text(size =10),
        strip.text = element_text(size = 9),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 10), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

plot_oc_ac_tp <- cowplot::plot_grid(plot_oc_tp, plot_ac_tp, nrow = 2, rel_heights = c(1, 1))

ggsave(plot_oc_ac_tp, file = "results/supplement_Fig.S3_occupancy_abundance_change_thermal_position.png", 
       width = 150, height = 115, units = 'mm', dpi = 300)


## combine maps with line plots
maps_plots <- cowplot::plot_grid(maps_population_change, plot_oc_ac_tp, nrow = 2, rel_heights = c(1.007, 1))

ggsave(maps_plots, file = "results/supplement_Fig.S3_maps_occupancy_abundance_change_thermal_position.png", 
       width = 180, height = 277, units = 'mm', dpi = 300) # 1:1.538867



######################
## figure showing relationships between occupancy/abundance change and baseline thermal position across taxa groups
rm(list = ls())
gc()

load("results/data_input_to_models.RDATA") # observation data
brm_oc_tmeanpos_taxa <- readRDS("models/brm_oc_tmeanpos_taxa.rds")
brm_ac_tmeanpos_taxa <- readRDS("models/brm_ac_tmeanpos_taxa.rds")

# sequence of thermal position
position_seq <- seq(0, 1, length = 101)

# new data for predicting relations
newdata_oc_taxa <- oc_period %>% 
  distinct(study, realm, realm_taxa, n_samp, duration, extent_km2) %>% 
  group_by(realm, realm_taxa) %>% 
  summarise(nstudy = n_distinct(study),
            n_samp = round(exp(mean(log(n_samp)))), # use the geometric mean 
            duration = round(exp(mean(log(duration)))),
            extent_km2 = round(exp(mean(log(extent_km2))), 1)) %>% 
  filter(nstudy >5) %>% 
  mutate(realm_taxa = as.character(realm_taxa)) %>% 
  cross_join(tibble(tmeanpos = position_seq)) 

newdata_ac_taxa <- ac_period %>% 
  distinct(study, realm, realm_taxa, n_samp, duration, extent_km2) %>% 
  group_by(realm, realm_taxa) %>% 
  summarise(nstudy = n_distinct(study),
            n_samp = round(exp(mean(log(n_samp)))), # use the geometric mean 
            duration = round(exp(mean(log(duration)))),
            extent_km2 = round(exp(mean(log(extent_km2))), 1)) %>% 
  filter(nstudy >5) %>% 
  mutate(realm_taxa = as.character(realm_taxa)) %>% 
  cross_join(tibble(tmeanpos = position_seq)) 

## text showing fixed effect
## for occupancy change 
beta_oc_tmeanpos_taxa <- fixef(brm_oc_tmeanpos_taxa, robust = TRUE, probs = c(0.025, 0.975)) %>%
  as.data.frame() %>% 
  rownames_to_column(var = "term") %>%
  filter(grepl(":tmeanpos", term)) %>%
  mutate(realm_taxa = gsub(":tmeanpos", "", term),
         realm_taxa = gsub("realm_taxa", "", realm_taxa)) %>%
  mutate(label = paste('beta == ', round(Estimate, 3),'~(', round(Q2.5, 3), "*','~", round(Q97.5, 3),')'))

beta_oc_tmeanpos_taxa <- beta_oc_tmeanpos_taxa %>% 
  mutate(realm_taxa = factor(realm_taxa, levels = c("Terrestrial_Birds", "Terrestrial_Invertebrates", "Terrestrial_Plants", "Freshwater_Fish",
                                        "Freshwater_Invertebrates", "Marine_Fish", "Marine_Invertebrates", "Marine_Plants")))

## for abundance change
beta_ac_tmeanpos_taxa <- fixef(brm_ac_tmeanpos_taxa, robust = TRUE, probs = c(0.025, 0.975)) %>%
  as.data.frame() %>% 
  rownames_to_column(var = "term") %>%
  filter(grepl(":tmeanpos", term)) %>%
  mutate(realm_taxa = gsub(":tmeanpos", "", term),
         realm_taxa = gsub("realm_taxa", "", realm_taxa)) %>%
  mutate(label = paste('beta == ', round(Estimate, 3),'~(', round(Q2.5, 3), "*','~", round(Q97.5, 3),')'))

beta_ac_tmeanpos_taxa <- beta_ac_tmeanpos_taxa %>% 
  mutate(realm_taxa = factor(realm_taxa, levels = c("Terrestrial_Birds", "Terrestrial_Invertebrates", "Terrestrial_Plants", "Freshwater_Fish",
                                        "Freshwater_Invertebrates", "Marine_Fish", "Marine_Invertebrates")))

## fitted values along thermal position
## for occupancy change
fitted_oc_tmeanpos_taxa <- fitted(brm_oc_tmeanpos_taxa, newdata = newdata_oc_taxa, re_formula = NA, ndraws = 1000, probs = c(0.025, 0.1, 0.9, 0.975)) %>%
  as_tibble()  %>% 
  bind_cols(newdata_oc_taxa)

fitted_oc_tmeanpos_taxa <- fitted_oc_tmeanpos_taxa %>% 
  mutate(realm_taxa = factor(realm_taxa, levels = c("Terrestrial_Birds", "Terrestrial_Invertebrates", "Terrestrial_Plants", "Freshwater_Fish",
                                        "Freshwater_Invertebrates", "Marine_Fish", "Marine_Invertebrates", "Marine_Plants")))

## for abundance change
fitted_ac_tmeanpos_taxa <- fitted(brm_ac_tmeanpos_taxa, newdata = newdata_ac_taxa, re_formula = NA, ndraws = 1000, probs = c(0.025, 0.1, 0.9, 0.975)) %>%
  as_tibble()  %>% 
  bind_cols(newdata_ac_taxa)

fitted_ac_tmeanpos_taxa <- fitted_ac_tmeanpos_taxa %>% 
  mutate(realm_taxa = factor(realm_taxa, levels = c("Terrestrial_Birds", "Terrestrial_Invertebrates", "Terrestrial_Plants", "Freshwater_Fish",
                                        "Freshwater_Invertebrates", "Marine_Fish", "Marine_Invertebrates")))

## draw figures
# colors for realm
realm_col <- c("Terrestrial" = "#D95F02", "Freshwater" = "#44bbff", "Marine" = "#7570B3")
realm_col <- c("All" = "#000000", "Terrestrial" = "#A65628", "Freshwater" = "#33A02C", "Marine" = "#1F78B4")

# label showing number of studies
nstudy_oc <- newdata_oc_taxa %>% 
  distinct(realm_taxa, nstudy) %>% 
  mutate(label = paste0("italic(n) == ", nstudy)) %>% 
  mutate(realm_taxa = factor(realm_taxa, levels = c("Terrestrial_Birds", "Terrestrial_Invertebrates", "Terrestrial_Plants", "Freshwater_Fish",
                                                    "Freshwater_Invertebrates", "Marine_Fish", "Marine_Invertebrates", "Marine_Plants")))

nstudy_ac <- newdata_ac_taxa %>% 
  distinct(realm_taxa, nstudy) %>% 
  mutate(label = paste0("italic(n) == ", nstudy)) %>% 
  mutate(realm_taxa = factor(realm_taxa, levels = c("Terrestrial_Birds", "Terrestrial_Invertebrates", "Terrestrial_Plants", "Freshwater_Fish",
                                                    "Freshwater_Invertebrates", "Marine_Fish", "Marine_Invertebrates")))

plot_oc <- ggplot(fitted_oc_tmeanpos_taxa) +
  facet_wrap(~ realm_taxa, nrow = 2) + 
  geom_hline(yintercept = 0, lty = 2, linewidth = 0.3) + 
  geom_ribbon(aes(x =  tmeanpos, ymin = Q2.5, ymax = Q97.5, fill = realm), alpha = 0.3) + 
  geom_line(aes(x = tmeanpos, y = Estimate, color = realm), lwd = 1) + 
  geom_text(data = beta_oc_tmeanpos_taxa, mapping = aes(x= 0.12, y = Inf, label = label), 
            hjust = 0, vjust = 1.5,  parse = T, size = 2.46) + 
  geom_text(data = nstudy_oc, mapping = aes(x= 0.12, y = Inf, label = label), 
            hjust = 0, vjust = 4,  parse = T, size = 2.46) + 
  # scale_y_continuous(limits = c(-0.28,0.4)) +
  labs(x = "Baseline thermal position", y = "Occupancy change") + 
  scale_colour_manual(name = NULL, values = realm_col) + 
  scale_fill_manual(name = NULL, values = realm_col) + 
  theme_bw()  +
  theme(legend.position = "no", 
        text = element_text(size =9),
        strip.text = element_text(size = 9),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 9), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

plot_ac <- ggplot(fitted_ac_tmeanpos_taxa) +
  facet_wrap(~ realm_taxa, nrow = 2) + 
  geom_hline(yintercept = 0, lty = 2, linewidth = 0.3) + 
  geom_ribbon(aes(x =  tmeanpos, ymin = Q2.5, ymax = Q97.5, fill = realm), alpha = 0.3) + 
  geom_line(aes(x = tmeanpos, y = Estimate, color = realm), lwd = 1) + 
  geom_text(data = beta_ac_tmeanpos_taxa, mapping = aes(x= 0.12, y = Inf, label = label), 
            hjust = 0, vjust = 1.5,  parse = T, size = 2.46) + 
  geom_text(data = nstudy_ac, mapping = aes(x= 0.12, y = Inf, label = label), 
            hjust = 0, vjust = 4,  parse = T, size = 2.46) + 
  # scale_y_continuous(limits = c(-0.28,0.4)) +
  labs(x = "Baseline thermal position", y = "Abundance change") + 
  scale_colour_manual(name = NULL, values = realm_col) + 
  scale_fill_manual(name = NULL, values = realm_col) + 
  theme_bw()  +
  theme(legend.position = "no", 
        text = element_text(size =9),
        strip.text = element_text(size = 9),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 9), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

plot_oc_ac <- cowplot::plot_grid(plot_oc, plot_ac,
                                 nrow = 2, labels = c("a", "b"), label_size = 10)

ggsave(plot_oc_ac, file = "results/supplement_Fig.S_occupancy_abundance_change_thermal_position_taxa.png", 
       width = 180, height = 200, units = 'mm', dpi = 600)



######################
## figure showing relationships between occupancy/abundance change and baseline thermal position across regions
rm(list = ls())
gc()

load("results/data_input_to_models.RDATA") # observation data
brm_oc_tmeanpos_region <- readRDS("models/brm_oc_tmeanpos_region.rds")
brm_ac_tmeanpos_region <- readRDS("models/brm_ac_tmeanpos_region.rds")

# sequence of thermal position
position_seq <- seq(0, 1, length = 101)

# new data for predicting relations
newdata_oc_region <- oc_period %>% 
  distinct(study, realm, realm_region, n_samp, duration, extent_km2) %>% 
  group_by(realm, realm_region) %>% 
  summarise(nstudy = n_distinct(study),
            n_samp = round(exp(mean(log(n_samp)))), # use the geometric mean 
            duration = round(exp(mean(log(duration)))),
            extent_km2 = round(exp(mean(log(extent_km2))), 1)) %>% 
  filter(nstudy >5) %>% 
  mutate(realm_region = as.character(realm_region),
         realm_region1 = gsub(" ", "", realm_region)) %>% 
  cross_join(tibble(tmeanpos = position_seq)) 

newdata_ac_region <- ac_period %>% 
  distinct(study, realm, realm_region, n_samp, duration, extent_km2) %>% 
  group_by(realm, realm_region) %>% 
  summarise(nstudy = n_distinct(study),
            n_samp = round(exp(mean(log(n_samp)))), # use the geometric mean 
            duration = round(exp(mean(log(duration)))),
            extent_km2 = round(exp(mean(log(extent_km2))), 1)) %>% 
  filter(nstudy >5) %>% 
  mutate(realm_region = as.character(realm_region),
         realm_region1 = gsub(" ", "", realm_region)) %>% 
  cross_join(tibble(tmeanpos = position_seq)) 

## text showing fixed effect
## for occupancy change 
beta_oc_tmeanpos_region <- fixef(brm_oc_tmeanpos_region, robust = TRUE, probs = c(0.025, 0.975)) %>%
  as.data.frame() %>% 
  rownames_to_column(var = "term") %>%
  filter(grepl(":tmeanpos", term)) %>%
  mutate(realm_region1 = gsub(":tmeanpos", "", term),
         realm_region1 = gsub("realm_region", "", realm_region1)) %>%
  mutate(label = paste('beta == ', round(Estimate, 3),'~(', round(Q2.5, 3), "*','~", round(Q97.5, 3),')'))

beta_oc_tmeanpos_region <- beta_oc_tmeanpos_region %>% 
  left_join(newdata_oc_region %>% distinct(realm_region, realm_region1)) %>% 
  mutate(realm_region = factor(realm_region, levels = c("Terrestrial_Africa", "Terrestrial_Australia", "Terrestrial_Europe", "Terrestrial_North America",
                                                    "Freshwater_Europe", "Freshwater_North America", "Marine_Australian Ocean","Marine_East Atlantic", 
                                                    "Marine_West Atlantic",  "Marine_Pacific Ocean")))

## for abundance change
beta_ac_tmeanpos_region <- fixef(brm_ac_tmeanpos_region, robust = TRUE, probs = c(0.025, 0.975)) %>%
  as.data.frame() %>% 
  rownames_to_column(var = "term") %>%
  filter(grepl(":tmeanpos", term)) %>%
  mutate(realm_region1 = gsub(":tmeanpos", "", term),
         realm_region1 = gsub("realm_region", "", realm_region1)) %>%
  mutate(label = paste('beta == ', round(Estimate, 3),'~(', round(Q2.5, 3), "*','~", round(Q97.5, 3),')'))

beta_ac_tmeanpos_region <- beta_ac_tmeanpos_region %>% 
  left_join(newdata_oc_region %>% distinct(realm_region, realm_region1)) %>% 
  mutate(realm_region = factor(realm_region, levels = c("Terrestrial_Australia", "Terrestrial_Europe", "Terrestrial_North America",
                                                        "Freshwater_Europe", "Freshwater_North America", "Marine_Australian Ocean",
                                                        "Marine_East Atlantic",  "Marine_West Atlantic",  "Marine_Pacific Ocean")))
## fitted values along thermal position
## for occupancy change
fitted_oc_tmeanpos_region <- fitted(brm_oc_tmeanpos_region, newdata = newdata_oc_region, re_formula = NA, ndraws = 1000, probs = c(0.025, 0.1, 0.9, 0.975)) %>%
  as_tibble()  %>% 
  bind_cols(newdata_oc_region)

fitted_oc_tmeanpos_region <- fitted_oc_tmeanpos_region %>% 
  mutate(realm_region = factor(realm_region, levels = c("Terrestrial_Africa", "Terrestrial_Australia", "Terrestrial_Europe", "Terrestrial_North America",
                                                          "Freshwater_Europe", "Freshwater_North America", "Marine_Australian Ocean","Marine_East Atlantic", 
                                                          "Marine_West Atlantic",  "Marine_Pacific Ocean")))

## for abundance change
fitted_ac_tmeanpos_region <- fitted(brm_ac_tmeanpos_region, newdata = newdata_ac_region, re_formula = NA, ndraws = 1000, probs = c(0.025, 0.1, 0.9, 0.975)) %>%
  as_tibble()  %>% 
  bind_cols(newdata_ac_region)

fitted_ac_tmeanpos_region <- fitted_ac_tmeanpos_region %>% 
  mutate(realm_region = factor(realm_region, levels = c("Terrestrial_Australia", "Terrestrial_Europe", "Terrestrial_North America",
                                                          "Freshwater_Europe", "Freshwater_North America", "Marine_Australian Ocean",
                                                        "Marine_East Atlantic",  "Marine_West Atlantic",  "Marine_Pacific Ocean")))

## draw figures
# colors for realm
realm_col <- c("Terrestrial" = "#D95F02", "Freshwater" = "#44bbff", "Marine" = "#7570B3")
realm_col <- c("All" = "#000000", "Terrestrial" = "#A65628", "Freshwater" = "#33A02C", "Marine" = "#1F78B4")

# label showing number of studies
nstudy_oc <- newdata_oc_region %>% 
  distinct(realm_region, nstudy) %>% 
  mutate(label = paste0("italic(n) == ", nstudy)) %>% 
  mutate(realm_region = factor(realm_region, levels = c("Terrestrial_Africa", "Terrestrial_Australia", "Terrestrial_Europe", "Terrestrial_North America",
                                                        "Freshwater_Europe", "Freshwater_North America", "Marine_Australian Ocean","Marine_East Atlantic", 
                                                        "Marine_West Atlantic",  "Marine_Pacific Ocean")))

nstudy_ac <- newdata_ac_region %>% 
  distinct(realm_region, nstudy) %>% 
  mutate(label = paste0("italic(n) == ", nstudy)) %>% 
  mutate(realm_region = factor(realm_region, levels = c("Terrestrial_Australia", "Terrestrial_Europe", "Terrestrial_North America",
                                                        "Freshwater_Europe", "Freshwater_North America", "Marine_Australian Ocean",
                                                        "Marine_East Atlantic",  "Marine_West Atlantic",  "Marine_Pacific Ocean")))

plot_oc <- ggplot(fitted_oc_tmeanpos_region) +
  facet_wrap(~ realm_region, nrow = 2) + 
  geom_hline(yintercept = 0, lty = 2, linewidth = 0.3) + 
  geom_ribbon(aes(x =  tmeanpos, ymin = Q2.5, ymax = Q97.5, fill = realm), alpha = 0.3) + 
  geom_line(aes(x = tmeanpos, y = Estimate, color = realm), lwd = 1) + 
  geom_text(data = beta_oc_tmeanpos_region, mapping = aes(x= 0.12, y = Inf, label = label), 
            hjust = 0, vjust = 1.5,  parse = T, size = 2.46) + 
  geom_text(data = nstudy_oc, mapping = aes(x= 0.12, y = Inf, label = label), 
            hjust = 0, vjust = 4,  parse = T, size = 2.46) + 
  # scale_y_continuous(limits = c(-0.28,0.4)) +
  labs(x = "Baseline thermal position", y = "Occupancy change") + 
  scale_colour_manual(name = NULL, values = realm_col) + 
  scale_fill_manual(name = NULL, values = realm_col) + 
  theme_bw()  +
  theme(legend.position = "no", 
        text = element_text(size =9),
        strip.text = element_text(size = 9),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 9), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

plot_ac <- ggplot(fitted_ac_tmeanpos_region) +
  facet_wrap(~ realm_region, nrow = 2) + 
  geom_hline(yintercept = 0, lty = 2, linewidth = 0.3) + 
  geom_ribbon(aes(x =  tmeanpos, ymin = Q2.5, ymax = Q97.5, fill = realm), alpha = 0.3) + 
  geom_line(aes(x = tmeanpos, y = Estimate, color = realm), lwd = 1) + 
  geom_text(data = beta_ac_tmeanpos_region, mapping = aes(x= 0.12, y = Inf, label = label), 
            hjust = 0, vjust = 1.5,  parse = T, size = 2.46) + 
  geom_text(data = nstudy_ac, mapping = aes(x= 0.12, y = Inf, label = label), 
            hjust = 0, vjust = 4,  parse = T, size = 2.46) + 
  # scale_y_continuous(limits = c(-0.28,0.4)) +
  labs(x = "Baseline thermal position", y = "Abundance change") + 
  scale_colour_manual(name = NULL, values = realm_col) + 
  scale_fill_manual(name = NULL, values = realm_col) + 
  theme_bw()  +
  theme(legend.position = "no", 
        text = element_text(size =9),
        strip.text = element_text(size = 9),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 9), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

plot_oc_ac <- cowplot::plot_grid(plot_oc, plot_ac,
                                 nrow = 2, labels = c("a", "b"), label_size = 10)

ggsave(plot_oc_ac, file = "results/supplement_Fig.S_occupancy_abundance_change_thermal_position_regions.png", 
       width = 220, height = 200, units = 'mm', dpi = 600)



######################
## figure showing the interaction between baseline thermal position and temperature change across taxa groups
rm(list = ls())
gc()

load("results/data_input_to_models.RDATA") # observation data
brm_oc_tmeanpos_tmeanchange_taxa <- readRDS("models/brm_oc_tmeanpos_tmeanchange_taxa.rds")
brm_ac_tmeanpos_tmeanchange_taxa <- readRDS("models/brm_ac_tmeanpos_tmeanchange_taxa.rds")

## prepare new data to predict occupancy or abundance change
# sequence of thermal position
position_seq <- seq(0, 1, length = 101)

# temperature change 
tchange_ranges <- dat_meta %>%
  summarise(Mean = mean(change_tempmean/(duration_mean-1)),
            sd = sd(change_tempmean/(duration_mean-1)),
            low = round((Mean - 2*sd), 3),
            high = round((Mean + 2*sd), 3)) %>%
  dplyr::select(-sd)

tchange_seq <- seq(tchange_ranges$low, tchange_ranges$high, length = 101)

# new data for prediction
newdata_oc_taxa <- oc_period %>% 
  distinct(study, realm, realm_taxa, n_samp, duration, extent_km2) %>% 
  group_by(realm, realm_taxa) %>% 
  summarise(nstudy = n_distinct(study),
            n_samp = round(exp(mean(log(n_samp)))), # use the geometric mean 
            duration = round(exp(mean(log(duration)))),
            extent_km2 = round(exp(mean(log(extent_km2))), 1)) %>% 
  filter(nstudy >5) %>% 
  mutate(realm_taxa = as.character(realm_taxa)) %>% 
  cross_join(tibble(tmeanpos = position_seq)) %>% 
  cross_join(tibble(change_tempmean = tchange_seq)) 

newdata_ac_taxa <- ac_period %>% 
  distinct(study, realm, realm_taxa, n_samp, duration, extent_km2) %>% 
  group_by(realm, realm_taxa) %>% 
  summarise(nstudy = n_distinct(study),
            n_samp = round(exp(mean(log(n_samp)))), # use the geometric mean 
            duration = round(exp(mean(log(duration)))),
            extent_km2 = round(exp(mean(log(extent_km2))), 1)) %>% 
  filter(nstudy >5) %>% 
  mutate(realm_taxa = as.character(realm_taxa)) %>% 
  cross_join(tibble(tmeanpos = position_seq)) %>% 
  cross_join(tibble(change_tempmean = tchange_seq)) 

## predicted occupancy and abundance change 
fitted_oc_tmeanpos_taxa <- fitted(brm_oc_tmeanpos_tmeanchange_taxa, newdata = newdata_oc_taxa, re_formula = NA, ndraws = 1000, probs = c(0.025, 0.975)) %>%
  as_tibble()  %>% 
  bind_cols(newdata_oc_taxa)

fitted_oc_tmeanpos_taxa <- fitted_oc_tmeanpos_taxa %>% 
  mutate(realm_taxa = factor(realm_taxa, levels = c("Terrestrial_Birds", "Terrestrial_Invertebrates", "Terrestrial_Plants", "Freshwater_Fish",
                                                    "Freshwater_Invertebrates", "Marine_Fish", "Marine_Invertebrates", "Marine_Plants")))

fitted_ac_tmeanpos_taxa <- fitted(brm_ac_tmeanpos_tmeanchange_taxa, newdata = newdata_ac_taxa, re_formula = NA, ndraws = 1000, probs = c(0.025, 0.975)) %>%
  as_tibble()  %>% 
  bind_cols(newdata_ac_taxa)

fitted_ac_tmeanpos_taxa <- fitted_ac_tmeanpos_taxa %>% 
  mutate(realm_taxa = factor(realm_taxa, levels = c("Terrestrial_Birds", "Terrestrial_Invertebrates", "Terrestrial_Plants", "Freshwater_Fish",
                                                    "Freshwater_Invertebrates", "Marine_Fish", "Marine_Invertebrates")))

## text showing interaction effect
beta_oc_interaction_taxa <- fixef(brm_oc_tmeanpos_tmeanchange_taxa) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "term") %>%
  filter(grepl("tmeanpos:change_tempmean", term)) %>%
  mutate(realm_taxa = gsub(":tmeanpos:change_tempmean", "", term),
         realm_taxa = gsub("realm_taxa", "", realm_taxa)) %>%
  mutate(realm_taxa = factor(realm_taxa, levels = c("Terrestrial_Birds", "Terrestrial_Invertebrates", "Terrestrial_Plants", "Freshwater_Fish",
                                                    "Freshwater_Invertebrates", "Marine_Fish", "Marine_Invertebrates", "Marine_Plants"))) %>% 
  mutate(label = paste('beta == ', round(Estimate, 2),'~(', round(Q2.5, 2), "*','~", round(Q97.5, 2),')'))

beta_ac_interaction_taxa <- fixef(brm_ac_tmeanpos_tmeanchange_taxa) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "term") %>%
  filter(grepl("tmeanpos:change_tempmean", term)) %>%
  mutate(realm_taxa = gsub(":tmeanpos:change_tempmean", "", term),
         realm_taxa = gsub("realm_taxa", "", realm_taxa)) %>% 
  mutate(realm_taxa = factor(realm_taxa, levels = c("Terrestrial_Birds", "Terrestrial_Invertebrates", "Terrestrial_Plants", "Freshwater_Fish",
                                                    "Freshwater_Invertebrates", "Marine_Fish", "Marine_Invertebrates"))) %>%
  mutate(label = paste('beta == ', round(Estimate, 2),'~(', round(Q2.5, 2), "*','~", round(Q97.5, 2),')'))

## label showing number of studies
nstudy_oc <- newdata_oc_taxa %>% 
  distinct(realm_taxa, nstudy) %>% 
  mutate(label = paste0("italic(n) == ", nstudy)) %>% 
  mutate(realm_taxa = factor(realm_taxa, levels = c("Terrestrial_Birds", "Terrestrial_Invertebrates", "Terrestrial_Plants", "Freshwater_Fish",
                                                    "Freshwater_Invertebrates", "Marine_Fish", "Marine_Invertebrates", "Marine_Plants")))
nstudy_ac <- newdata_ac_taxa %>% 
  distinct(realm_taxa, nstudy) %>% 
  mutate(label = paste0("italic(n) == ", nstudy)) %>% 
  mutate(realm_taxa = factor(realm_taxa, levels = c("Terrestrial_Birds", "Terrestrial_Invertebrates", "Terrestrial_Plants", "Freshwater_Fish",
                                                    "Freshwater_Invertebrates", "Marine_Fish", "Marine_Invertebrates")))

## draw figures
plot_oc <- ggplot(fitted_oc_tmeanpos_taxa) + 
  facet_wrap(~ realm_taxa, nrow = 2) +
  geom_line(aes(x = tmeanpos, y = Estimate, color = change_tempmean, group = change_tempmean)) +
  geom_hline(yintercept = 0, lty =2, lwd = 0.4) + 
  geom_text(data = beta_oc_interaction_taxa, mapping = aes(x= 0.1, y = Inf, label = label), 
            hjust = 0, vjust = 1.5,  parse = T, size = 2.46) + 
  geom_text(data = nstudy_oc, mapping = aes(x= 0.1, y = Inf, label = label), 
            hjust = 0, vjust = 4,  parse = T, size = 2.46) + 
  labs(tag = "a", x = "Baseline thermal position", y = "Occupancy change", 
       color = expression('Temperature change ('*degree*'C/year)')) + 
  scale_color_gradient2(low = "blue", mid = "darkgray", high="red") +
  theme_bw() + 
  theme(legend.position = c(0.75, 1.15), 
        plot.margin = unit(c(0.8,0.2,0,0.2), "cm"),
        legend.direction = 'horizontal',
        legend.key.size = unit(4, 'mm'),
        legend.background = element_blank(),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 9),
        text = element_text(size = 10),
        strip.text = element_text(size = 9),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  guides(color = guide_colourbar(barwidth = 6, barheight = 0.8))

plot_ac <- ggplot(fitted_ac_tmeanpos_taxa) + 
  facet_wrap(~ realm_taxa, nrow = 2) +
  geom_line(aes(x = tmeanpos, y = Estimate, color = change_tempmean, group = change_tempmean)) +
  geom_hline(yintercept = 0, lty =2, lwd = 0.4) + 
  geom_text(data = beta_ac_interaction_taxa, mapping = aes(x= 0.1, y = Inf, label = label), 
            hjust = 0, vjust = 1.5,  parse = T, size = 2.46) + 
  geom_text(data = nstudy_ac, mapping = aes(x= 0.1, y = Inf, label = label), 
            hjust = 0, vjust = 4,  parse = T, size = 2.46) + 
  labs(tag = "b", x = "Baseline thermal position", y = "Abundance change") + 
  scale_color_gradient2(low = "blue", mid = "darkgray", high="red") +
  theme_bw() + 
  theme(legend.position = "no",
        plot.margin = unit(c(0,0.2,0.2,0.2), "cm"),
        text = element_text(size =10),
        strip.text = element_text(size = 9),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

cowplot::plot_grid(plot_oc, plot_ac,
                   nrow = 2, rel_heights = c(1.08, 1), align = "v")

ggsave(file = "results/supplement_Fig.S_interaction_taxa.png", 
       width = 180, height = 200, units = 'mm', dpi = 300)



######################
## figure showing the interaction between baseline thermal position and temperature change across regions
rm(list = ls())
gc()

load("results/data_input_to_models.RDATA") # observation data
brm_oc_tmeanpos_tmeanchange_region <- readRDS("models/brm_oc_tmeanpos_tmeanchange_region.rds")
brm_ac_tmeanpos_tmeanchange_region <- readRDS("models/brm_ac_tmeanpos_tmeanchange_region.rds")

## prepare new data to predict occupancy or abundance change
# sequence of thermal position
position_seq <- seq(0, 1, length = 101)

# temperature change 
tchange_ranges <- dat_meta %>%
  summarise(Mean = mean(change_tempmean/(duration_mean-1)),
            sd = sd(change_tempmean/(duration_mean-1)),
            low = round((Mean - 2*sd), 3),
            high = round((Mean + 2*sd), 3)) %>%
  dplyr::select(-sd)

tchange_seq <- seq(tchange_ranges$low, tchange_ranges$high, length = 101)

# new data for prediction
newdata_oc_region <- oc_period %>% 
  distinct(study, realm, realm_region, n_samp, duration, extent_km2) %>% 
  group_by(realm, realm_region) %>% 
  summarise(nstudy = n_distinct(study),
            n_samp = round(exp(mean(log(n_samp)))), # use the geometric mean 
            duration = round(exp(mean(log(duration)))),
            extent_km2 = round(exp(mean(log(extent_km2))), 1)) %>% 
  filter(nstudy >5) %>% 
  mutate(realm_region = as.character(realm_region),
         realm_region1 = gsub(" ", "", realm_region)) %>% 
  cross_join(tibble(tmeanpos = position_seq)) %>% 
  cross_join(tibble(change_tempmean = tchange_seq)) 

newdata_ac_region <- ac_period %>% 
  distinct(study, realm, realm_region, n_samp, duration, extent_km2) %>% 
  group_by(realm, realm_region) %>% 
  summarise(nstudy = n_distinct(study),
            n_samp = round(exp(mean(log(n_samp)))), # use the geometric mean 
            duration = round(exp(mean(log(duration)))),
            extent_km2 = round(exp(mean(log(extent_km2))), 1)) %>% 
  filter(nstudy >5) %>% 
  mutate(realm_region = as.character(realm_region),
         realm_region1 = gsub(" ", "", realm_region)) %>% 
  cross_join(tibble(tmeanpos = position_seq)) %>% 
  cross_join(tibble(change_tempmean = tchange_seq)) 

## predicted occupancy and abundance change 
fitted_oc_tmeanpos_region <- fitted(brm_oc_tmeanpos_tmeanchange_region, newdata = newdata_oc_region, re_formula = NA, ndraws = 1000, probs = c(0.025, 0.975)) %>%
  as_tibble()  %>% 
  bind_cols(newdata_oc_region)

fitted_oc_tmeanpos_region <- fitted_oc_tmeanpos_region %>% 
  mutate(realm_region = factor(realm_region, levels = c("Terrestrial_Africa", "Terrestrial_Australia", "Terrestrial_Europe", "Terrestrial_North America",
                                                        "Freshwater_Europe", "Freshwater_North America", "Marine_Australian Ocean","Marine_East Atlantic", 
                                                        "Marine_West Atlantic",  "Marine_Pacific Ocean")))

fitted_ac_tmeanpos_region <- fitted(brm_ac_tmeanpos_tmeanchange_region, newdata = newdata_ac_region, re_formula = NA, ndraws = 1000, probs = c(0.025, 0.975)) %>%
  as_tibble()  %>% 
  bind_cols(newdata_ac_region)

fitted_ac_tmeanpos_region <- fitted_ac_tmeanpos_region %>% 
  mutate(realm_region = factor(realm_region, levels = c("Terrestrial_Australia", "Terrestrial_Europe", "Terrestrial_North America",
                                                        "Freshwater_Europe", "Freshwater_North America", "Marine_Australian Ocean",
                                                        "Marine_East Atlantic",  "Marine_West Atlantic",  "Marine_Pacific Ocean")))

## text showing interaction effect
beta_oc_interaction_region <- fixef(brm_oc_tmeanpos_tmeanchange_region) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "term") %>%
  filter(grepl("tmeanpos:change_tempmean", term)) %>%
  mutate(realm_region1 = gsub(":tmeanpos:change_tempmean", "", term),
         realm_region1 = gsub("realm_region", "", realm_region1)) %>%
  left_join(newdata_oc_region %>% distinct(realm_region, realm_region1)) %>% 
  mutate(realm_region = factor(realm_region, levels = c("Terrestrial_Africa", "Terrestrial_Australia", "Terrestrial_Europe", "Terrestrial_North America",
                                                        "Freshwater_Europe", "Freshwater_North America", "Marine_Australian Ocean","Marine_East Atlantic", 
                                                        "Marine_West Atlantic",  "Marine_Pacific Ocean"))) %>% 
  mutate(label = paste('beta == ', round(Estimate, 2),'~(', round(Q2.5, 2), "*','~", round(Q97.5, 2),')'))

beta_ac_interaction_region <- fixef(brm_ac_tmeanpos_tmeanchange_region) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "term") %>%
  filter(grepl("tmeanpos:change_tempmean", term)) %>%
  mutate(realm_region1 = gsub(":tmeanpos:change_tempmean", "", term),
         realm_region1 = gsub("realm_region", "", realm_region1)) %>% 
  left_join(newdata_oc_region %>% distinct(realm_region, realm_region1)) %>% 
  mutate(realm_region = factor(realm_region, levels = c("Terrestrial_Australia", "Terrestrial_Europe", "Terrestrial_North America",
                                                        "Freshwater_Europe", "Freshwater_North America", "Marine_Australian Ocean",
                                                        "Marine_East Atlantic",  "Marine_West Atlantic",  "Marine_Pacific Ocean"))) %>% 
  mutate(label = paste('beta == ', round(Estimate, 2),'~(', round(Q2.5, 2), "*','~", round(Q97.5, 2),')'))

## label showing number of studies
nstudy_oc <- newdata_oc_region %>% 
  distinct(realm_region, nstudy) %>% 
  mutate(label = paste0("italic(n) == ", nstudy)) %>% 
  mutate(realm_region = factor(realm_region, levels = c("Terrestrial_Africa", "Terrestrial_Australia", "Terrestrial_Europe", "Terrestrial_North America",
                                                        "Freshwater_Europe", "Freshwater_North America", "Marine_Australian Ocean","Marine_East Atlantic", 
                                                        "Marine_West Atlantic",  "Marine_Pacific Ocean")))
nstudy_ac <- newdata_ac_region %>% 
  distinct(realm_region, nstudy) %>% 
  mutate(label = paste0("italic(n) == ", nstudy)) %>% 
  mutate(realm_region = factor(realm_region, levels = c("Terrestrial_Australia", "Terrestrial_Europe", "Terrestrial_North America",
                                                        "Freshwater_Europe", "Freshwater_North America", "Marine_Australian Ocean",
                                                        "Marine_East Atlantic",  "Marine_West Atlantic",  "Marine_Pacific Ocean")))

## draw figures
#  for visualization, remove extreme occupancy change values at high baseline thermal positions in Terrestrial-Africa
fitted_oc_tmeanpos_region1 <- fitted_oc_tmeanpos_region %>%  
  filter(!(realm_region == "Terrestrial_Africa" & tmeanpos > 0.4))

plot_oc <- ggplot(fitted_oc_tmeanpos_region1) + 
  facet_wrap(~ realm_region, nrow = 2) +
  geom_line(aes(x = tmeanpos, y = Estimate, color = change_tempmean, group = change_tempmean)) +
  geom_hline(yintercept = 0, lty =2, lwd = 0.4) + 
  geom_text(data = beta_oc_interaction_region, mapping = aes(x= 0.1, y = Inf, label = label), 
            hjust = 0, vjust = 1.5,  parse = T, size = 2.46) + 
  geom_text(data = nstudy_oc, mapping = aes(x= 0.1, y = Inf, label = label), 
            hjust = 0, vjust = 4,  parse = T, size = 2.46) + 
  labs(tag = "a", x = "Baseline thermal position", y = "Occupancy change", 
       color = expression('Temperature change ('*degree*'C/year)')) + 
  scale_color_gradient2(low = "blue", mid = "darkgray", high="red") +
  theme_bw() + 
  theme(legend.position = c(0.8, 1.15), 
        plot.margin = unit(c(0.8,0.2,0,0.2), "cm"),
        legend.direction = 'horizontal',
        legend.key.size = unit(4, 'mm'),
        legend.background = element_blank(),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 9),
        text = element_text(size = 10),
        strip.text = element_text(size = 9),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  guides(color = guide_colourbar(barwidth = 6, barheight = 0.8))

plot_ac <- ggplot(fitted_ac_tmeanpos_region) + 
  facet_wrap(~ realm_region, nrow = 2) +
  geom_line(aes(x = tmeanpos, y = Estimate, color = change_tempmean, group = change_tempmean)) +
  geom_hline(yintercept = 0, lty =2, lwd = 0.4) + 
  geom_text(data = beta_ac_interaction_region, mapping = aes(x= 0.1, y = Inf, label = label), 
            hjust = 0, vjust = 1.5,  parse = T, size = 2.46) + 
  geom_text(data = nstudy_ac, mapping = aes(x= 0.1, y = Inf, label = label), 
            hjust = 0, vjust = 4,  parse = T, size = 2.46) + 
  labs(tag = "b", x = "Baseline thermal position", y = "Abundance change") + 
  scale_color_gradient2(low = "blue", mid = "darkgray", high="red") +
  theme_bw() + 
  theme(legend.position = "no",
        plot.margin = unit(c(0,0.2,0.2,0.2), "cm"),
        text = element_text(size =10),
        strip.text = element_text(size = 9),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

cowplot::plot_grid(plot_oc, plot_ac,
                   nrow = 2, rel_heights = c(1.08, 1), align = "v")

ggsave(file = "results/supplement_Fig.S_interaction_region.png", 
       width = 220, height = 200, units = 'mm', dpi = 300)



######################
## sensitivity analyses based on species with thermal range > 5-degrees and at least 20 occurrences
## figure showing the relationship between occupancy/abundance change and thermal position 

rm(list = ls())
gc()

load("results/data_input_to_models.RDATA") # observation data
brm_oc_tmeanpos_widesp <- readRDS("models/brm_oc_tmeanpos_widesp.rds")
brm_ac_tmeanpos_widesp <- readRDS("models/brm_ac_tmeanpos_widesp.rds")


## new data for predicting relations across realms
position_seq <- seq(0, 1, length = 101) # sequence of thermal position

newdata <- oc_period %>% distinct(study, realm, n_samp, duration, extent_km2) %>% group_by(realm) %>% 
  summarise(n_samp = round(exp(mean(log(n_samp)))), # use the geometric mean 
            duration = round(exp(mean(log(duration)))),
            extent_km2 = round(exp(mean(log(extent_km2))), 1)) %>% 
  left_join(tibble(tmeanpos = rep(position_seq, 3), 
                   realm = rep(c("Terrestrial", "Freshwater", "Marine"), each = length(position_seq))))

## text showing fixed effect
beta_oc_tmeanpos <- fixef(brm_oc_tmeanpos_widesp, robust = TRUE, probs = c(0.025, 0.975)) %>%
  as.data.frame() %>% 
  rownames_to_column(var = "term") %>%
  filter(grepl(":tmeanpos", term)) %>%
  mutate(realm = gsub(":tmeanpos", "", term),
         realm = gsub("realm", "", realm),
         realm = factor(realm, levels = c("Terrestrial", "Freshwater", "Marine"))) %>%
  mutate(label = paste('beta == ', round(Estimate, 3),'~(', round(Q2.5, 3), "*','~", round(Q97.5, 3),')'))

beta_ac_tmeanpos <- fixef(brm_ac_tmeanpos_widesp, robust = TRUE, probs = c(0.025, 0.975)) %>%
  as.data.frame() %>% 
  rownames_to_column(var = "term") %>%
  filter(grepl(":tmeanpos", term)) %>%
  mutate(realm = gsub(":tmeanpos", "", term),
         realm = gsub("realm", "", realm),
         realm = factor(realm, levels = c("Terrestrial", "Freshwater", "Marine"))) %>%
  mutate(label = paste('beta == ', round(Estimate, 3),'~(', round(Q2.5, 3), "*','~", round(Q97.5, 3),')'))

## fitted values along thermal position
fitted_oc_tmeanpos <- fitted(brm_oc_tmeanpos_widesp, newdata = newdata, re_formula = NA, ndraws = 1000, probs = c(0.025, 0.1, 0.9, 0.975)) %>%
  as_tibble()  %>% 
  bind_cols(newdata) %>% 
  mutate(realm = factor(realm, levels = c("Terrestrial", "Freshwater", "Marine")))

fitted_ac_tmeanpos <- fitted(brm_ac_tmeanpos_widesp, newdata = newdata, re_formula = NA, ndraws = 1000, probs = c(0.025, 0.1, 0.9, 0.975)) %>%
  as_tibble()  %>% 
  bind_cols(newdata) %>% 
  mutate(realm = factor(realm, levels = c("Terrestrial", "Freshwater", "Marine")))


## draw figures
# colors for realm
realm_col <- c("Terrestrial" = "#D95F02", "Freshwater" = "#44bbff", "Marine" = "#7570B3")
realm_col <- c("All" = "#000000", "Terrestrial" = "#A65628", "Freshwater" = "#33A02C", "Marine" = "#1F78B4")

plot_oc <- ggplot(fitted_oc_tmeanpos) +
  facet_wrap(~ realm, nrow = 1) + 
  geom_hline(yintercept = 0, lty = 2, linewidth = 0.3) + 
  geom_ribbon(aes(x =  tmeanpos, ymin = Q2.5, ymax = Q97.5, fill = realm), alpha = 0.3) + 
  geom_line(aes(x = tmeanpos, y = Estimate, color = realm), lwd = 1) + 
  geom_text(data = beta_oc_tmeanpos, mapping = aes(x= 0.12, y = Inf, label = label), 
            hjust = 0, vjust = 1.5,  parse = T, size = 2.46) + 
  labs(x = "Baseline thermal position", y = "Occupancy change") + 
  scale_colour_manual(name = NULL, values = realm_col) + 
  scale_fill_manual(name = NULL, values = realm_col) + 
  theme_bw()  +
  theme(legend.position = "no", 
        text = element_text(size =9),
        strip.text = element_text(size = 9),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 9), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

plot_ac <- ggplot(fitted_ac_tmeanpos) +
  facet_wrap(~ realm, nrow = 1) + 
  geom_hline(yintercept = 0, lty = 2, linewidth = 0.3) + 
  geom_ribbon(aes(x =  tmeanpos, ymin = Q2.5, ymax = Q97.5, fill = realm), alpha = 0.3) + 
  geom_line(aes(x = tmeanpos, y = Estimate,  color = realm), lwd = 1) + 
  geom_text(data = beta_ac_tmeanpos, mapping = aes(x= 0.12, y = Inf, label = label), 
            hjust = 0, vjust = 1.5,  parse = T, size = 2.46) + 
  # scale_y_continuous(limits = c(-2.5,2.5)) +
  labs(x = "Baseline thermal position", y = "Abundance change") + 
  scale_colour_manual(name = NULL, values = realm_col) + 
  scale_fill_manual(name = NULL, values = realm_col) + 
  theme_bw()  +
  theme(legend.position = "no", 
        text = element_text(size =9),
        strip.text = element_text(size = 9),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 9), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

plot_oc_ac <- cowplot::plot_grid(plot_oc, plot_ac, nrow = 2, labels = c("a", "b"), label_size = 10)

ggsave(plot_oc_ac, file = "results/Supplement_Fig.S_effect_thermal_position_wideRangeSpecies.png", 
       width = 140, height = 110, units = 'mm', dpi = 600) 



## figure showing the interaction between baseline thermal position and temperature change
rm(list = ls())
gc()

load("results/data_input_to_models.RDATA") # observation data
brm_oc_tmeanpos_tmeanchange_widesp <- readRDS("models/brm_oc_tmeanpos_tmeanchange_widesp.rds")
brm_ac_tmeanpos_tmeanchange_widesp <- readRDS("models/brm_ac_tmeanpos_tmeanchange_widesp.rds")

## prepare new data to predict occupancy or abundance change
# sequence of thermal position
position_seq <- seq(0, 1, length = 101)

# temperature change 
tchange_ranges <- dat_meta %>%
  summarise(Mean = mean(change_tempmean/(duration_mean-1)),
            sd = sd(change_tempmean/(duration_mean-1)),
            low = round((Mean - 2*sd), 3),
            high = round((Mean + 2*sd), 3)) %>%
  dplyr::select(-sd)

tchange_seq <- seq(tchange_ranges$low, tchange_ranges$high, length = 101)
tchange_realm <- tibble(expand.grid(change_tempmean = tchange_seq, realm = c("Terrestrial", "Freshwater", "Marine")))

# new data for prediction
newdata <- tchange_realm %>% 
  left_join(expand.grid(change_tempmean = unique(tchange_realm$change_tempmean), tmeanpos = position_seq), relationship = "many-to-many") %>%
  left_join(oc_period %>% distinct(study, realm, n_samp, duration, extent_km2) %>% group_by(realm) %>% 
              summarise(n_samp = round(exp(mean(log(n_samp)))), # use the geometric mean 
                        duration = round(exp(mean(log(duration)))),
                        extent_km2 = round(exp(mean(log(extent_km2))), 1))) %>% 
  mutate(realm = factor(realm, levels = c("Terrestrial", "Freshwater", "Marine")))

## predicted occupancy and abundance change 
fitted_oc_tmeanpos <- fitted(brm_oc_tmeanpos_tmeanchange_widesp, newdata = newdata, re_formula = NA, ndraws = 1000, probs = c(0.025, 0.975)) %>%
  as_tibble()  %>% 
  bind_cols(newdata)

fitted_ac_tmeanpos <- fitted(brm_ac_tmeanpos_tmeanchange_widesp, newdata = newdata, re_formula = NA, ndraws = 1000, probs = c(0.025, 0.975)) %>%
  as_tibble()  %>% 
  bind_cols(newdata)

# text showing interaction effect
beta_oc_interaction_text <- fixef(brm_oc_tmeanpos_tmeanchange_widesp) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "term") %>%
  filter(grepl("tmeanpos:change_tempmean", term)) %>%
  mutate(realm = gsub(":tmeanpos:change_tempmean", "", term),
         realm = gsub("realm", "", realm),
         realm = factor(realm, levels = c("Terrestrial", "Freshwater", "Marine"))) %>%
  mutate(label = paste('beta == ', round(Estimate, 2),'~(', round(Q2.5, 2), "*','~", round(Q97.5, 2),')'))

beta_ac_interaction_text <- fixef(brm_ac_tmeanpos_tmeanchange_widesp) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "term") %>%
  filter(grepl("tmeanpos:change_tempmean", term)) %>%
  mutate(realm = gsub(":tmeanpos:change_tempmean", "", term),
         realm = gsub("realm", "", realm),
         realm = factor(realm, levels = c("Terrestrial", "Freshwater", "Marine"))) %>%
  mutate(label = paste('beta == ', round(Estimate, 2),'~(', round(Q2.5, 2), "*','~", round(Q97.5, 2),')'))


## draw figures
plot_oc <- ggplot(fitted_oc_tmeanpos) + 
  facet_wrap(~ realm, scales = "fixed") +
  geom_line(aes(x = tmeanpos, y = Estimate, color = change_tempmean, group = change_tempmean)) +
  geom_hline(yintercept = 0, lty =2, lwd = 0.4) + 
  geom_text(data = beta_oc_interaction_text, mapping = aes(x= 0.1, y = Inf, label = label), 
            hjust = 0, vjust = 1.5,  parse = T, size = 2.46) + 
  labs(tag = "a", x = "Baseline thermal position", y = "Occupancy change", 
       color = expression('Temperature change ('*degree*'C/year)')) + 
  scale_color_gradient2(low = "blue", mid = "darkgray", high="red") +
  theme_bw() + 
  theme(legend.position = c(0.7, 1.3), 
        plot.margin = unit(c(0.6,0.2,0,0.2), "cm"),
        legend.direction = 'horizontal',
        legend.key.size = unit(4, 'mm'),
        legend.background = element_blank(),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 9),
        text = element_text(size = 10),
        strip.text = element_text(size = 9),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  guides(color = guide_colourbar(barwidth = 4, barheight = 0.5))

plot_ac <- ggplot(fitted_ac_tmeanpos) + 
  facet_wrap(~ realm, scales = "fixed") +
  geom_line(aes(x = tmeanpos, y = Estimate, color = change_tempmean, group = change_tempmean)) +
  geom_hline(yintercept = 0, lty =2, lwd = 0.4) + 
  geom_text(data = beta_ac_interaction_text, mapping = aes(x= 0.1, y = Inf, label = label), 
            hjust = 0, vjust = 1.5,  parse = T, size = 2.46) + 
  labs(tag = "b", x = "Baseline thermal position", y = "Abundance change") + 
  scale_color_gradient2(low = "blue", mid = "darkgray", high="red") +
  theme_bw() + 
  theme(legend.position = "no",
        plot.margin = unit(c(0,0.2,0.2,0.2), "cm"),
        text = element_text(size =10),
        strip.text = element_text(size = 9),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

cowplot::plot_grid(plot_oc, plot_ac,
                   nrow = 2, rel_heights = c(1.06, 1), align = "v")

ggsave(file = "results/supplement_Fig.S_interaction_wideRangeSpecies.png", 
       width = 150, height = 120, units = 'mm', dpi = 300)



######################
## sensitivity analyses based on monthly maximum temperature
## figure showing the relationship between occupancy/abundance change and thermal position 

rm(list = ls())
gc()

load("results/data_input_to_models.RDATA") # observation data
brm_oc_tmaxpos_realm <- readRDS("models/brm_oc_tmaxpos_realm.rds")
brm_ac_tmaxpos_realm <- readRDS("models/brm_ac_tmaxpos_realm.rds")


## new data for predicting relations across realms
position_seq <- seq(0, 1, length = 101) # sequence of thermal position

newdata <- oc_period %>% distinct(study, realm, n_samp, duration, extent_km2) %>% group_by(realm) %>% 
  summarise(n_samp = round(exp(mean(log(n_samp)))), # use the geometric mean 
            duration = round(exp(mean(log(duration)))),
            extent_km2 = round(exp(mean(log(extent_km2))), 1)) %>% 
  left_join(tibble(tmaxpos = rep(position_seq, 3), 
                   realm = rep(c("Terrestrial", "Freshwater", "Marine"), each = length(position_seq))))

## text showing fixed effect
beta_oc_tmaxpos <- fixef(brm_oc_tmaxpos_realm, robust = TRUE, probs = c(0.025, 0.975)) %>%
  as.data.frame() %>% 
  rownames_to_column(var = "term") %>%
  filter(grepl(":tmaxpos", term)) %>%
  mutate(realm = gsub(":tmaxpos", "", term),
         realm = gsub("realm", "", realm),
         realm = factor(realm, levels = c("Terrestrial", "Freshwater", "Marine"))) %>%
  mutate(label = paste('beta == ', round(Estimate, 3),'~(', round(Q2.5, 3), "*','~", round(Q97.5, 3),')'))

beta_ac_tmaxpos <- fixef(brm_ac_tmaxpos_realm, robust = TRUE, probs = c(0.025, 0.975)) %>%
  as.data.frame() %>% 
  rownames_to_column(var = "term") %>%
  filter(grepl(":tmaxpos", term)) %>%
  mutate(realm = gsub(":tmaxpos", "", term),
         realm = gsub("realm", "", realm),
         realm = factor(realm, levels = c("Terrestrial", "Freshwater", "Marine"))) %>%
  mutate(label = paste('beta == ', round(Estimate, 3),'~(', round(Q2.5, 3), "*','~", round(Q97.5, 3),')'))

## fitted values along thermal position
fitted_oc_tmaxpos <- fitted(brm_oc_tmaxpos_realm, newdata = newdata, re_formula = NA, ndraws = 1000, probs = c(0.025, 0.1, 0.9, 0.975)) %>%
  as_tibble()  %>% 
  bind_cols(newdata) %>% 
  mutate(realm = factor(realm, levels = c("Terrestrial", "Freshwater", "Marine")))

fitted_ac_tmaxpos <- fitted(brm_ac_tmaxpos_realm, newdata = newdata, re_formula = NA, ndraws = 1000, probs = c(0.025, 0.1, 0.9, 0.975)) %>%
  as_tibble()  %>% 
  bind_cols(newdata) %>% 
  mutate(realm = factor(realm, levels = c("Terrestrial", "Freshwater", "Marine")))


## draw figures
# colors for realm
realm_col <- c("Terrestrial" = "#D95F02", "Freshwater" = "#44bbff", "Marine" = "#7570B3")
realm_col <- c("All" = "#000000", "Terrestrial" = "#A65628", "Freshwater" = "#33A02C", "Marine" = "#1F78B4")

plot_oc <- ggplot(fitted_oc_tmaxpos) +
  facet_wrap(~ realm, nrow = 1) + 
  geom_hline(yintercept = 0, lty = 2, linewidth = 0.3) + 
  geom_ribbon(aes(x =  tmaxpos, ymin = Q2.5, ymax = Q97.5, fill = realm), alpha = 0.3) + 
  geom_line(aes(x = tmaxpos, y = Estimate, color = realm), lwd = 1) + 
  geom_text(data = beta_oc_tmaxpos, mapping = aes(x= 0.12, y = Inf, label = label), 
            hjust = 0, vjust = 1.5,  parse = T, size = 2.46) + 
  labs(x = "Baseline thermal position", y = "Occupancy change") + 
  scale_colour_manual(name = NULL, values = realm_col) + 
  scale_fill_manual(name = NULL, values = realm_col) + 
  theme_bw()  +
  theme(legend.position = "no", 
        text = element_text(size =9),
        strip.text = element_text(size = 9),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 9), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

plot_ac <- ggplot(fitted_ac_tmaxpos) +
  facet_wrap(~ realm, nrow = 1) + 
  geom_hline(yintercept = 0, lty = 2, linewidth = 0.3) + 
  geom_ribbon(aes(x =  tmaxpos, ymin = Q2.5, ymax = Q97.5, fill = realm), alpha = 0.3) + 
  geom_line(aes(x = tmaxpos, y = Estimate,  color = realm), lwd = 1) + 
  geom_text(data = beta_ac_tmaxpos, mapping = aes(x= 0.12, y = Inf, label = label), 
            hjust = 0, vjust = 1.5,  parse = T, size = 2.46) + 
  # scale_y_continuous(limits = c(-2.5,2.5)) +
  labs(x = "Baseline thermal position", y = "Abundance change") + 
  scale_colour_manual(name = NULL, values = realm_col) + 
  scale_fill_manual(name = NULL, values = realm_col) + 
  theme_bw()  +
  theme(legend.position = "no", 
        text = element_text(size =9),
        strip.text = element_text(size = 9),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 9), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

plot_oc_ac <- cowplot::plot_grid(plot_oc, plot_ac, nrow = 2, labels = c("a", "b"), label_size = 10)

ggsave(plot_oc_ac, file = "results/Supplement_Fig.S_effect_thermal_position_maximumTemperature.png", 
       width = 140, height = 110, units = 'mm', dpi = 600) 



## figure showing the interaction between baseline thermal position and temperature change
rm(list = ls())
gc()

load("results/data_input_to_models.RDATA") # observation data
brm_oc_tmaxpos_tmaxchange_realm <- readRDS("models/brm_oc_tmaxpos_tmaxchange_realm.rds")
brm_ac_tmaxpos_tmaxchange_realm <- readRDS("models/brm_ac_tmaxpos_tmaxchange_realm.rds")

## prepare new data to predict occupancy or abundance change
# sequence of thermal position
position_seq <- seq(0, 1, length = 101)

# temperature change 
tchange_ranges <- dat_meta %>%
  summarise(Mean = mean(change_tempmax/(duration_mean-1)),
            sd = sd(change_tempmax/(duration_mean-1)),
            low = round((Mean - 2*sd), 3),
            high = round((Mean + 2*sd), 3)) %>%
  dplyr::select(-sd)

tchange_seq <- seq(tchange_ranges$low, tchange_ranges$high, length = 101)
tchange_realm <- tibble(expand.grid(change_tempmax = tchange_seq, realm = c("Terrestrial", "Freshwater", "Marine")))

# new data for prediction
newdata <- tchange_realm %>% 
  left_join(expand.grid(change_tempmax = unique(tchange_realm$change_tempmax), tmaxpos = position_seq), relationship = "many-to-many") %>%
  left_join(oc_period %>% distinct(study, realm, n_samp, duration, extent_km2) %>% group_by(realm) %>% 
              summarise(n_samp = round(exp(mean(log(n_samp)))), # use the geometric mean 
                        duration = round(exp(mean(log(duration)))),
                        extent_km2 = round(exp(mean(log(extent_km2))), 1))) %>% 
  mutate(realm = factor(realm, levels = c("Terrestrial", "Freshwater", "Marine")))

## predicted occupancy and abundance change 
fitted_oc_tmaxpos <- fitted(brm_oc_tmaxpos_tmaxchange_realm, newdata = newdata, re_formula = NA, ndraws = 1000, probs = c(0.025, 0.975)) %>%
  as_tibble()  %>% 
  bind_cols(newdata)

fitted_ac_tmaxpos <- fitted(brm_ac_tmaxpos_tmaxchange_realm, newdata = newdata, re_formula = NA, ndraws = 1000, probs = c(0.025, 0.975)) %>%
  as_tibble()  %>% 
  bind_cols(newdata)

# text showing interaction effect
beta_oc_interaction_text <- fixef(brm_oc_tmaxpos_tmaxchange_realm) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "term") %>%
  filter(grepl("tmaxpos:change_tempmax", term)) %>%
  mutate(realm = gsub(":tmaxpos:change_tempmax", "", term),
         realm = gsub("realm", "", realm),
         realm = factor(realm, levels = c("Terrestrial", "Freshwater", "Marine"))) %>%
  mutate(label = paste('beta == ', round(Estimate, 2),'~(', round(Q2.5, 2), "*','~", round(Q97.5, 2),')'))

beta_ac_interaction_text <- fixef(brm_ac_tmaxpos_tmaxchange_realm) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "term") %>%
  filter(grepl("tmaxpos:change_tempmax", term)) %>%
  mutate(realm = gsub(":tmaxpos:change_tempmax", "", term),
         realm = gsub("realm", "", realm),
         realm = factor(realm, levels = c("Terrestrial", "Freshwater", "Marine"))) %>%
  mutate(label = paste('beta == ', round(Estimate, 2),'~(', round(Q2.5, 2), "*','~", round(Q97.5, 2),')'))


## draw figures
plot_oc <- ggplot(fitted_oc_tmaxpos) + 
  facet_wrap(~ realm, scales = "fixed") +
  geom_line(aes(x = tmaxpos, y = Estimate, color = change_tempmax, group = change_tempmax)) +
  geom_hline(yintercept = 0, lty =2, lwd = 0.4) + 
  geom_text(data = beta_oc_interaction_text, mapping = aes(x= 0.2, y = Inf, label = label), 
            hjust = 0, vjust = 1.5,  parse = T, size = 2.46) + 
  labs(tag = "a", x = "Baseline thermal position", y = "Occupancy change", 
       color = expression('Temperature change ('*degree*'C/year)')) + 
  scale_color_gradient2(breaks = c(-0.07, 0,0.07, 0.14),low = "blue", mid = "darkgray", high="red") +
  theme_bw() + 
  theme(legend.position = c(0.7, 1.3), 
        plot.margin = unit(c(0.6,0.2,0,0.2), "cm"),
        legend.direction = 'horizontal',
        legend.key.size = unit(4, 'mm'),
        legend.background = element_blank(),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 9),
        text = element_text(size = 10),
        strip.text = element_text(size = 9),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  guides(color = guide_colourbar(barwidth = 5, barheight = 0.5))

plot_ac <- ggplot(fitted_ac_tmaxpos) + 
  facet_wrap(~ realm, scales = "fixed") +
  geom_line(aes(x = tmaxpos, y = Estimate, color = change_tempmax, group = change_tempmax)) +
  geom_hline(yintercept = 0, lty =2, lwd = 0.4) + 
  geom_text(data = beta_ac_interaction_text, mapping = aes(x= 0.2, y = Inf, label = label), 
            hjust = 0, vjust = 1.5,  parse = T, size = 2.46) + 
  labs(tag = "b", x = "Baseline thermal position", y = "Abundance change") + 
  scale_color_gradient2(low = "blue", mid = "darkgray", high="red") +
  theme_bw() + 
  theme(legend.position = "no",
        plot.margin = unit(c(0,0.2,0.2,0.2), "cm"),
        text = element_text(size =10),
        strip.text = element_text(size = 9),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

cowplot::plot_grid(plot_oc, plot_ac,
                   nrow = 2, rel_heights = c(1.06, 1), align = "v")

ggsave(file = "results/supplement_Fig.S_interaction_maximumTemperature.png", 
       width = 150, height = 120, units = 'mm', dpi = 300)




######################
## sensitivity analyses based on monthly minimum temperature
## figure showing the relationship between occupancy/abundance change and thermal position 

rm(list = ls())
gc()

load("results/data_input_to_models.RDATA") # observation data
brm_oc_tminpos_realm <- readRDS("models/brm_oc_tminpos_realm.rds")
brm_ac_tminpos_realm <- readRDS("models/brm_ac_tminpos_realm.rds")


## new data for predicting relations across realms
position_seq <- seq(0, 1, length = 101) # sequence of thermal position

newdata <- oc_period %>% distinct(study, realm, n_samp, duration, extent_km2) %>% group_by(realm) %>% 
  summarise(n_samp = round(exp(mean(log(n_samp)))), # use the geometric mean 
            duration = round(exp(mean(log(duration)))),
            extent_km2 = round(exp(mean(log(extent_km2))), 1)) %>% 
  left_join(tibble(tminpos = rep(position_seq, 3), 
                   realm = rep(c("Terrestrial", "Freshwater", "Marine"), each = length(position_seq))))

## text showing fixed effect
beta_oc_tminpos <- fixef(brm_oc_tminpos_realm, robust = TRUE, probs = c(0.025, 0.975)) %>%
  as.data.frame() %>% 
  rownames_to_column(var = "term") %>%
  filter(grepl(":tminpos", term)) %>%
  mutate(realm = gsub(":tminpos", "", term),
         realm = gsub("realm", "", realm),
         realm = factor(realm, levels = c("Terrestrial", "Freshwater", "Marine"))) %>%
  mutate(label = paste('beta == ', round(Estimate, 3),'~(', round(Q2.5, 3), "*','~", round(Q97.5, 3),')'))

beta_ac_tminpos <- fixef(brm_ac_tminpos_realm, robust = TRUE, probs = c(0.025, 0.975)) %>%
  as.data.frame() %>% 
  rownames_to_column(var = "term") %>%
  filter(grepl(":tminpos", term)) %>%
  mutate(realm = gsub(":tminpos", "", term),
         realm = gsub("realm", "", realm),
         realm = factor(realm, levels = c("Terrestrial", "Freshwater", "Marine"))) %>%
  mutate(label = paste('beta == ', round(Estimate, 3),'~(', round(Q2.5, 3), "*','~", round(Q97.5, 3),')'))

## fitted values along thermal position
fitted_oc_tminpos <- fitted(brm_oc_tminpos_realm, newdata = newdata, re_formula = NA, ndraws = 1000, probs = c(0.025, 0.1, 0.9, 0.975)) %>%
  as_tibble()  %>% 
  bind_cols(newdata) %>% 
  mutate(realm = factor(realm, levels = c("Terrestrial", "Freshwater", "Marine")))

fitted_ac_tminpos <- fitted(brm_ac_tminpos_realm, newdata = newdata, re_formula = NA, ndraws = 1000, probs = c(0.025, 0.1, 0.9, 0.975)) %>%
  as_tibble()  %>% 
  bind_cols(newdata) %>% 
  mutate(realm = factor(realm, levels = c("Terrestrial", "Freshwater", "Marine")))


## draw figures
# colors for realm
realm_col <- c("Terrestrial" = "#D95F02", "Freshwater" = "#44bbff", "Marine" = "#7570B3")
realm_col <- c("All" = "#000000", "Terrestrial" = "#A65628", "Freshwater" = "#33A02C", "Marine" = "#1F78B4")

plot_oc <- ggplot(fitted_oc_tminpos) +
  facet_wrap(~ realm, nrow = 1) + 
  geom_hline(yintercept = 0, lty = 2, linewidth = 0.3) + 
  geom_ribbon(aes(x =  tminpos, ymin = Q2.5, ymax = Q97.5, fill = realm), alpha = 0.3) + 
  geom_line(aes(x = tminpos, y = Estimate, color = realm), lwd = 1) + 
  geom_text(data = beta_oc_tminpos, mapping = aes(x= 0.12, y = Inf, label = label), 
            hjust = 0, vjust = 1.5,  parse = T, size = 2.46) + 
  labs(x = "Baseline thermal position", y = "Occupancy change") + 
  scale_colour_manual(name = NULL, values = realm_col) + 
  scale_fill_manual(name = NULL, values = realm_col) + 
  theme_bw()  +
  theme(legend.position = "no", 
        text = element_text(size =9),
        strip.text = element_text(size = 9),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 9), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

plot_ac <- ggplot(fitted_ac_tminpos) +
  facet_wrap(~ realm, nrow = 1) + 
  geom_hline(yintercept = 0, lty = 2, linewidth = 0.3) + 
  geom_ribbon(aes(x =  tminpos, ymin = Q2.5, ymax = Q97.5, fill = realm), alpha = 0.3) + 
  geom_line(aes(x = tminpos, y = Estimate,  color = realm), lwd = 1) + 
  geom_text(data = beta_ac_tminpos, mapping = aes(x= 0.12, y = Inf, label = label), 
            hjust = 0, vjust = 1.5,  parse = T, size = 2.46) + 
  # scale_y_continuous(limits = c(-2.5,2.5)) +
  labs(x = "Baseline thermal position", y = "Abundance change") + 
  scale_colour_manual(name = NULL, values = realm_col) + 
  scale_fill_manual(name = NULL, values = realm_col) + 
  theme_bw()  +
  theme(legend.position = "no", 
        text = element_text(size =9),
        strip.text = element_text(size = 9),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 9), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

plot_oc_ac <- cowplot::plot_grid(plot_oc, plot_ac, nrow = 2, labels = c("a", "b"), label_size = 10)

ggsave(plot_oc_ac, file = "results/Supplement_Fig.S_effect_thermal_position_minimumTemperature.png", 
       width = 140, height = 110, units = 'mm', dpi = 600) 



## figure showing the interaction between baseline thermal position and temperature change
rm(list = ls())
gc()

load("results/data_input_to_models.RDATA") # observation data
brm_oc_tminpos_tminchange_realm <- readRDS("models/brm_oc_tminpos_tminchange_realm.rds")
brm_ac_tminpos_tminchange_realm <- readRDS("models/brm_ac_tminpos_tminchange_realm.rds")

## prepare new data to predict occupancy or abundance change
# sequence of thermal position
position_seq <- seq(0, 1, length = 101)

# temperature change 
tchange_ranges <- dat_meta %>%
  summarise(Mean = mean(change_tempmin/(duration_mean-1)),
            sd = sd(change_tempmin/(duration_mean-1)),
            low = round((Mean - 2*sd), 3),
            high = round((Mean + 2*sd), 3)) %>%
  dplyr::select(-sd)

tchange_seq <- seq(tchange_ranges$low, tchange_ranges$high, length = 101)
tchange_realm <- tibble(expand.grid(change_tempmin = tchange_seq, realm = c("Terrestrial", "Freshwater", "Marine")))

# new data for prediction
newdata <- tchange_realm %>% 
  left_join(expand.grid(change_tempmin = unique(tchange_realm$change_tempmin), tminpos = position_seq), relationship = "many-to-many") %>%
  left_join(oc_period %>% distinct(study, realm, n_samp, duration, extent_km2) %>% group_by(realm) %>% 
              summarise(n_samp = round(exp(mean(log(n_samp)))), # use the geometric mean 
                        duration = round(exp(mean(log(duration)))),
                        extent_km2 = round(exp(mean(log(extent_km2))), 1))) %>% 
  mutate(realm = factor(realm, levels = c("Terrestrial", "Freshwater", "Marine")))

## predicted occupancy and abundance change 
fitted_oc_tminpos <- fitted(brm_oc_tminpos_tminchange_realm, newdata = newdata, re_formula = NA, ndraws = 1000, probs = c(0.025, 0.975)) %>%
  as_tibble()  %>% 
  bind_cols(newdata)

fitted_ac_tminpos <- fitted(brm_ac_tminpos_tminchange_realm, newdata = newdata, re_formula = NA, ndraws = 1000, probs = c(0.025, 0.975)) %>%
  as_tibble()  %>% 
  bind_cols(newdata)

# text showing interaction effect
beta_oc_interaction_text <- fixef(brm_oc_tminpos_tminchange_realm) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "term") %>%
  filter(grepl("tminpos:change_tempmin", term)) %>%
  mutate(realm = gsub(":tminpos:change_tempmin", "", term),
         realm = gsub("realm", "", realm),
         realm = factor(realm, levels = c("Terrestrial", "Freshwater", "Marine"))) %>%
  mutate(label = paste('beta == ', round(Estimate, 2),'~(', round(Q2.5, 2), "*','~", round(Q97.5, 2),')'))

beta_ac_interaction_text <- fixef(brm_ac_tminpos_tminchange_realm) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "term") %>%
  filter(grepl("tminpos:change_tempmin", term)) %>%
  mutate(realm = gsub(":tminpos:change_tempmin", "", term),
         realm = gsub("realm", "", realm),
         realm = factor(realm, levels = c("Terrestrial", "Freshwater", "Marine"))) %>%
  mutate(label = paste('beta == ', round(Estimate, 2),'~(', round(Q2.5, 2), "*','~", round(Q97.5, 2),')'))


## draw figures
plot_oc <- ggplot(fitted_oc_tminpos) + 
  facet_wrap(~ realm, scales = "fixed") +
  geom_line(aes(x = tminpos, y = Estimate, color = change_tempmin, group = change_tempmin)) +
  geom_hline(yintercept = 0, lty =2, lwd = 0.4) + 
  geom_text(data = beta_oc_interaction_text, mapping = aes(x= 0.2, y = Inf, label = label), 
            hjust = 0, vjust = 1.5,  parse = T, size = 2.46) + 
  labs(tag = "a", x = "Baseline thermal position", y = "Occupancy change", 
       color = expression('Temperature change ('*degree*'C/year)')) + 
  scale_color_gradient2(low = "blue", mid = "darkgray", high="red") +
  theme_bw() + 
  theme(legend.position = c(0.7, 1.3), 
        plot.margin = unit(c(0.6,0.2,0,0.2), "cm"),
        legend.direction = 'horizontal',
        legend.key.size = unit(4, 'mm'),
        legend.background = element_blank(),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 9),
        text = element_text(size = 10),
        strip.text = element_text(size = 9),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  guides(color = guide_colourbar(barwidth = 4, barheight = 0.5))

plot_ac <- ggplot(fitted_ac_tminpos) + 
  facet_wrap(~ realm, scales = "fixed") +
  geom_line(aes(x = tminpos, y = Estimate, color = change_tempmin, group = change_tempmin)) +
  geom_hline(yintercept = 0, lty =2, lwd = 0.4) + 
  geom_text(data = beta_ac_interaction_text, mapping = aes(x= 0.2, y = Inf, label = label), 
            hjust = 0, vjust = 1.5,  parse = T, size = 2.46) + 
  labs(tag = "b", x = "Baseline thermal position", y = "Abundance change") + 
  scale_color_gradient2(low = "blue", mid = "darkgray", high="red") +
  theme_bw() + 
  theme(legend.position = "no",
        plot.margin = unit(c(0,0.2,0.2,0.2), "cm"),
        text = element_text(size =10),
        strip.text = element_text(size = 9),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

cowplot::plot_grid(plot_oc, plot_ac,
                   nrow = 2, rel_heights = c(1.06, 1), align = "v")

ggsave(file = "results/supplement_Fig.S_interaction_minimumTemperature.png", 
       width = 150, height = 120, units = 'mm', dpi = 300)



######################
## sensitivity analyses to analyze occupancy changes using only the data with abundance data
## figure showing the relationship between occupancychange and thermal position 

rm(list = ls())
gc()

load("results/data_input_to_models.RDATA") # observation data
brm_oc_tmeanpos_abundancedata <- readRDS("models/brm_oc_tmeanpos_abundancedata.rds")

## new data for predicting relations across realms
position_seq <- seq(0, 1, length = 101) # sequence of thermal position

newdata <- oc_period %>% distinct(study, realm, n_samp, duration, extent_km2) %>% group_by(realm) %>% 
  summarise(n_samp = round(exp(mean(log(n_samp)))), # use the geometric mean 
            duration = round(exp(mean(log(duration)))),
            extent_km2 = round(exp(mean(log(extent_km2))), 1)) %>% 
  left_join(tibble(tmeanpos = rep(position_seq, 3), 
                   realm = rep(c("Terrestrial", "Freshwater", "Marine"), each = length(position_seq))))

## text showing fixed effect
beta_oc_tmeanpos <- fixef(brm_oc_tmeanpos_abundancedata, robust = TRUE, probs = c(0.025, 0.975)) %>%
  as.data.frame() %>% 
  rownames_to_column(var = "term") %>%
  filter(grepl(":tmeanpos", term)) %>%
  mutate(realm = gsub(":tmeanpos", "", term),
         realm = gsub("realm", "", realm),
         realm = factor(realm, levels = c("Terrestrial", "Freshwater", "Marine"))) %>%
  mutate(label = paste('beta == ', round(Estimate, 3),'~(', round(Q2.5, 3), "*','~", round(Q97.5, 3),')'))

## fitted values along thermal position
fitted_oc_tmeanpos <- fitted(brm_oc_tmeanpos_abundancedata, newdata = newdata, re_formula = NA, ndraws = 1000, probs = c(0.025, 0.1, 0.9, 0.975)) %>%
  as_tibble()  %>% 
  bind_cols(newdata) %>% 
  mutate(realm = factor(realm, levels = c("Terrestrial", "Freshwater", "Marine")))

## draw figures
# colors for realm
realm_col <- c("Terrestrial" = "#D95F02", "Freshwater" = "#44bbff", "Marine" = "#7570B3")
realm_col <- c("All" = "#000000", "Terrestrial" = "#A65628", "Freshwater" = "#33A02C", "Marine" = "#1F78B4")

plot_oc <- ggplot(fitted_oc_tmeanpos) +
  facet_wrap(~ realm, nrow = 1) + 
  geom_hline(yintercept = 0, lty = 2, linewidth = 0.3) + 
  geom_ribbon(aes(x =  tmeanpos, ymin = Q2.5, ymax = Q97.5, fill = realm), alpha = 0.3) + 
  geom_line(aes(x = tmeanpos, y = Estimate, color = realm), lwd = 1) + 
  geom_text(data = beta_oc_tmeanpos, mapping = aes(x= 0.12, y = Inf, label = label), 
            hjust = 0, vjust = 1.5,  parse = T, size = 2.46) + 
  labs(x = "Baseline thermal position", y = "Occupancy change") + 
  scale_colour_manual(name = NULL, values = realm_col) + 
  scale_fill_manual(name = NULL, values = realm_col) + 
  theme_bw()  +
  theme(legend.position = "no", 
        text = element_text(size =9),
        strip.text = element_text(size = 9),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 9), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

ggsave(plot_oc, file = "results/Supplement_Fig.S_effect_thermal_position_SubsetAbundanceData.png", 
       width = 140, height = 50, units = 'mm', dpi = 600) 



## figure showing the interaction between baseline thermal position and temperature change
rm(list = ls())
gc()

load("results/data_input_to_models.RDATA") # observation data
brm_oc_tmeanpos_tmeanchange_abundancedata <- readRDS("models/brm_oc_tmeanpos_tmeanchange_abundancedata.rds")

## prepare new data to predict occupancy or abundance change
# sequence of thermal position
position_seq <- seq(0, 1, length = 101)

# temperature change 
tchange_ranges <- dat_meta %>%
  summarise(Mean = mean(change_tempmean/(duration_mean-1)),
            sd = sd(change_tempmean/(duration_mean-1)),
            low = round((Mean - 2*sd), 3),
            high = round((Mean + 2*sd), 3)) %>%
  dplyr::select(-sd)

tchange_seq <- seq(tchange_ranges$low, tchange_ranges$high, length = 101)
tchange_realm <- tibble(expand.grid(change_tempmean = tchange_seq, realm = c("Terrestrial", "Freshwater", "Marine")))

# new data for prediction
newdata <- tchange_realm %>% 
  left_join(expand.grid(change_tempmean = unique(tchange_realm$change_tempmean), tmeanpos = position_seq), relationship = "many-to-many") %>%
  left_join(oc_period %>% distinct(study, realm, n_samp, duration, extent_km2) %>% group_by(realm) %>% 
              summarise(n_samp = round(exp(mean(log(n_samp)))), # use the geometric mean 
                        duration = round(exp(mean(log(duration)))),
                        extent_km2 = round(exp(mean(log(extent_km2))), 1))) %>% 
  mutate(realm = factor(realm, levels = c("Terrestrial", "Freshwater", "Marine")))

## predicted occupancy and abundance change 
fitted_oc_tmeanpos <- fitted(brm_oc_tmeanpos_tmeanchange_abundancedata, newdata = newdata, re_formula = NA, ndraws = 1000, probs = c(0.025, 0.975)) %>%
  as_tibble()  %>% 
  bind_cols(newdata)

# text showing interaction effect
beta_oc_interaction_text <- fixef(brm_oc_tmeanpos_tmeanchange_abundancedata) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "term") %>%
  filter(grepl("tmeanpos:change_tempmean", term)) %>%
  mutate(realm = gsub(":tmeanpos:change_tempmean", "", term),
         realm = gsub("realm", "", realm),
         realm = factor(realm, levels = c("Terrestrial", "Freshwater", "Marine"))) %>%
  mutate(label = paste('beta == ', round(Estimate, 2),'~(', round(Q2.5, 2), "*','~", round(Q97.5, 2),')'))

## draw figures
plot_oc <- ggplot(fitted_oc_tmeanpos) + 
  facet_wrap(~ realm, scales = "fixed") +
  geom_line(aes(x = tmeanpos, y = Estimate, color = change_tempmean, group = change_tempmean)) +
  geom_hline(yintercept = 0, lty =2, lwd = 0.4) + 
  geom_text(data = beta_oc_interaction_text, mapping = aes(x= 0.1, y = Inf, label = label), 
            hjust = 0, vjust = 1.5,  parse = T, size = 2.46) + 
  labs(x = "Baseline thermal position", y = "Occupancy change", 
       color = expression('Temperature change ('*degree*'C/year)')) + 
  scale_color_gradient2(low = "blue", mid = "darkgray", high="red") +
  theme_bw() + 
  theme(legend.position = c(0.7, 1.3), 
        plot.margin = unit(c(0.9,0.2,0,0.2), "cm"),
        legend.direction = 'horizontal',
        legend.key.size = unit(4, 'mm'),
        legend.background = element_blank(),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 9),
        text = element_text(size = 10),
        strip.text = element_text(size = 9),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  guides(color = guide_colourbar(barwidth = 4, barheight = 0.5))

ggsave(plot_oc, file = "results/supplement_Fig.S_interaction_subsetAbundanceData.png", 
       width = 150, height = 60, units = 'mm', dpi = 300)



######################
## sensitivity analyses based on data with spatial extents smaller than 10000 km2
## figure showing the relationship between occupancy/abundance change and thermal position 

rm(list = ls())
gc()

load("results/data_input_to_models.RDATA") # observation data
brm_oc_tmeanpos_smregion <- readRDS("models/brm_oc_tmeanpos_smregion.rds")
brm_ac_tmeanpos_smregion <- readRDS("models/brm_ac_tmeanpos_smregion.rds")

## new data for predicting relations across realms
position_seq <- seq(0, 1, length = 101) # sequence of thermal position

newdata <- oc_period %>% 
  filter(extent_km2 < 10000) %>% 
  distinct(study, realm, n_samp, duration, extent_km2) %>% group_by(realm) %>% 
  summarise(n_samp = round(exp(mean(log(n_samp)))), # use the geometric mean 
            duration = round(exp(mean(log(duration)))),
            extent_km2 = round(exp(mean(log(extent_km2))), 1)) %>% 
  left_join(tibble(tmeanpos = rep(position_seq, 3), 
                   realm = rep(c("Terrestrial", "Freshwater", "Marine"), each = length(position_seq))))

## text showing fixed effect
beta_oc_tmeanpos <- fixef(brm_oc_tmeanpos_smregion, robust = TRUE, probs = c(0.025, 0.975)) %>%
  as.data.frame() %>% 
  rownames_to_column(var = "term") %>%
  filter(grepl(":tmeanpos", term)) %>%
  mutate(realm = gsub(":tmeanpos", "", term),
         realm = gsub("realm", "", realm),
         realm = factor(realm, levels = c("Terrestrial", "Freshwater", "Marine"))) %>%
  mutate(label = paste('beta == ', round(Estimate, 3),'~(', round(Q2.5, 3), "*','~", round(Q97.5, 3),')'))

beta_ac_tmeanpos <- fixef(brm_ac_tmeanpos_smregion, robust = TRUE, probs = c(0.025, 0.975)) %>%
  as.data.frame() %>% 
  rownames_to_column(var = "term") %>%
  filter(grepl(":tmeanpos", term)) %>%
  mutate(realm = gsub(":tmeanpos", "", term),
         realm = gsub("realm", "", realm),
         realm = factor(realm, levels = c("Terrestrial", "Freshwater", "Marine"))) %>%
  mutate(label = paste('beta == ', round(Estimate, 3),'~(', round(Q2.5, 3), "*','~", round(Q97.5, 3),')'))

## fitted values along thermal position
fitted_oc_tmeanpos <- fitted(brm_oc_tmeanpos_smregion, newdata = newdata, re_formula = NA, ndraws = 1000, probs = c(0.025, 0.1, 0.9, 0.975)) %>%
  as_tibble()  %>% 
  bind_cols(newdata) %>% 
  mutate(realm = factor(realm, levels = c("Terrestrial", "Freshwater", "Marine")))

fitted_ac_tmeanpos <- fitted(brm_ac_tmeanpos_smregion, newdata = newdata, re_formula = NA, ndraws = 1000, probs = c(0.025, 0.1, 0.9, 0.975)) %>%
  as_tibble()  %>% 
  bind_cols(newdata) %>% 
  mutate(realm = factor(realm, levels = c("Terrestrial", "Freshwater", "Marine")))


## draw figures
# colors for realm
realm_col <- c("Terrestrial" = "#D95F02", "Freshwater" = "#44bbff", "Marine" = "#7570B3")
realm_col <- c("All" = "#000000", "Terrestrial" = "#A65628", "Freshwater" = "#33A02C", "Marine" = "#1F78B4")

plot_oc <- ggplot(fitted_oc_tmeanpos) +
  facet_wrap(~ realm, nrow = 1) + 
  geom_hline(yintercept = 0, lty = 2, linewidth = 0.3) + 
  geom_ribbon(aes(x =  tmeanpos, ymin = Q2.5, ymax = Q97.5, fill = realm), alpha = 0.3) + 
  geom_line(aes(x = tmeanpos, y = Estimate, color = realm), lwd = 1) + 
  geom_text(data = beta_oc_tmeanpos, mapping = aes(x= 0.12, y = Inf, label = label), 
            hjust = 0, vjust = 1.5,  parse = T, size = 2.46) + 
  labs(x = "Baseline thermal position", y = "Occupancy change") + 
  scale_colour_manual(name = NULL, values = realm_col) + 
  scale_fill_manual(name = NULL, values = realm_col) + 
  theme_bw()  +
  theme(legend.position = "no", 
        text = element_text(size =9),
        strip.text = element_text(size = 9),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 9), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

plot_ac <- ggplot(fitted_ac_tmeanpos) +
  facet_wrap(~ realm, nrow = 1) + 
  geom_hline(yintercept = 0, lty = 2, linewidth = 0.3) + 
  geom_ribbon(aes(x =  tmeanpos, ymin = Q2.5, ymax = Q97.5, fill = realm), alpha = 0.3) + 
  geom_line(aes(x = tmeanpos, y = Estimate,  color = realm), lwd = 1) + 
  geom_text(data = beta_ac_tmeanpos, mapping = aes(x= 0.12, y = Inf, label = label), 
            hjust = 0, vjust = 1.5,  parse = T, size = 2.46) + 
  # scale_y_continuous(limits = c(-2.5,2.5)) +
  labs(x = "Baseline thermal position", y = "Abundance change") + 
  scale_colour_manual(name = NULL, values = realm_col) + 
  scale_fill_manual(name = NULL, values = realm_col) + 
  theme_bw()  +
  theme(legend.position = "no", 
        text = element_text(size =9),
        strip.text = element_text(size = 9),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 9), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

plot_oc_ac <- cowplot::plot_grid(plot_oc, plot_ac, nrow = 2, labels = c("a", "b"), label_size = 10)

ggsave(plot_oc_ac, file = "results/Supplement_Fig.S_effect_thermal_position_subsetSmallRegion.png", 
       width = 140, height = 110, units = 'mm', dpi = 600) 



## figure showing the interaction between baseline thermal position and temperature change
rm(list = ls())
gc()

load("results/data_input_to_models.RDATA") # observation data
brm_oc_tmeanpos_tmeanchange_smregion <- readRDS("models/brm_oc_tmeanpos_tmeanchange_smregion.rds")
brm_ac_tmeanpos_tmeanchange_smregion <- readRDS("models/brm_ac_tmeanpos_tmeanchange_smregion.rds")

## prepare new data to predict occupancy or abundance change
# sequence of thermal position
position_seq <- seq(0, 1, length = 101)

# temperature change 
tchange_ranges <- dat_meta %>%
  summarise(Mean = mean(change_tempmean/(duration_mean-1)),
            sd = sd(change_tempmean/(duration_mean-1)),
            low = round((Mean - 2*sd), 3),
            high = round((Mean + 2*sd), 3)) %>%
  dplyr::select(-sd)

tchange_seq <- seq(tchange_ranges$low, tchange_ranges$high, length = 101)
tchange_realm <- tibble(expand.grid(change_tempmean = tchange_seq, realm = c("Terrestrial", "Freshwater", "Marine")))

# new data for prediction
newdata <- tchange_realm %>% 
  left_join(expand.grid(change_tempmean = unique(tchange_realm$change_tempmean), tmeanpos = position_seq), relationship = "many-to-many") %>%
  left_join(oc_period %>% 
              filter(extent_km2 < 10000) %>% 
              distinct(study, realm, n_samp, duration, extent_km2) %>% group_by(realm) %>% 
              summarise(n_samp = round(exp(mean(log(n_samp)))), # use the geometric mean 
                        duration = round(exp(mean(log(duration)))),
                        extent_km2 = round(exp(mean(log(extent_km2))), 1))) %>% 
  mutate(realm = factor(realm, levels = c("Terrestrial", "Freshwater", "Marine")))

## predicted occupancy and abundance change 
fitted_oc_tmeanpos <- fitted(brm_oc_tmeanpos_tmeanchange_smregion, newdata = newdata, re_formula = NA, ndraws = 1000, probs = c(0.025, 0.975)) %>%
  as_tibble()  %>% 
  bind_cols(newdata)

fitted_ac_tmeanpos <- fitted(brm_ac_tmeanpos_tmeanchange_smregion, newdata = newdata, re_formula = NA, ndraws = 1000, probs = c(0.025, 0.975)) %>%
  as_tibble()  %>% 
  bind_cols(newdata)

# text showing interaction effect
beta_oc_interaction_text <- fixef(brm_oc_tmeanpos_tmeanchange_smregion) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "term") %>%
  filter(grepl("tmeanpos:change_tempmean", term)) %>%
  mutate(realm = gsub(":tmeanpos:change_tempmean", "", term),
         realm = gsub("realm", "", realm),
         realm = factor(realm, levels = c("Terrestrial", "Freshwater", "Marine"))) %>%
  mutate(label = paste('beta == ', round(Estimate, 2),'~(', round(Q2.5, 2), "*','~", round(Q97.5, 2),')'))

beta_ac_interaction_text <- fixef(brm_ac_tmeanpos_tmeanchange_smregion) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "term") %>%
  filter(grepl("tmeanpos:change_tempmean", term)) %>%
  mutate(realm = gsub(":tmeanpos:change_tempmean", "", term),
         realm = gsub("realm", "", realm),
         realm = factor(realm, levels = c("Terrestrial", "Freshwater", "Marine"))) %>%
  mutate(label = paste('beta == ', round(Estimate, 2),'~(', round(Q2.5, 2), "*','~", round(Q97.5, 2),')'))


## draw figures
plot_oc <- ggplot(fitted_oc_tmeanpos) + 
  facet_wrap(~ realm, scales = "fixed") +
  geom_line(aes(x = tmeanpos, y = Estimate, color = change_tempmean, group = change_tempmean)) +
  geom_hline(yintercept = 0, lty =2, lwd = 0.4) + 
  geom_text(data = beta_oc_interaction_text, mapping = aes(x= 0.1, y = Inf, label = label), 
            hjust = 0, vjust = 1.5,  parse = T, size = 2.46) + 
  labs(tag = "a", x = "Baseline thermal position", y = "Occupancy change", 
       color = expression('Temperature change ('*degree*'C/year)')) + 
  scale_color_gradient2(low = "blue", mid = "darkgray", high="red") +
  theme_bw() + 
  theme(legend.position = c(0.7, 1.3), 
        plot.margin = unit(c(0.6,0.2,0,0.2), "cm"),
        legend.direction = 'horizontal',
        legend.key.size = unit(4, 'mm'),
        legend.background = element_blank(),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 9),
        text = element_text(size = 10),
        strip.text = element_text(size = 9),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  guides(color = guide_colourbar(barwidth = 4, barheight = 0.5))

plot_ac <- ggplot(fitted_ac_tmeanpos) + 
  facet_wrap(~ realm, scales = "fixed") +
  geom_line(aes(x = tmeanpos, y = Estimate, color = change_tempmean, group = change_tempmean)) +
  geom_hline(yintercept = 0, lty =2, lwd = 0.4) + 
  geom_text(data = beta_ac_interaction_text, mapping = aes(x= 0.1, y = Inf, label = label), 
            hjust = 0, vjust = 1.5,  parse = T, size = 2.46) + 
  labs(tag = "b", x = "Baseline thermal position", y = "Abundance change") + 
  scale_color_gradient2(low = "blue", mid = "darkgray", high="red") +
  theme_bw() + 
  theme(legend.position = "no",
        plot.margin = unit(c(0,0.2,0.2,0.2), "cm"),
        text = element_text(size =10),
        strip.text = element_text(size = 9),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

plot_oc_ac <- cowplot::plot_grid(plot_oc, plot_ac,
                   nrow = 2, rel_heights = c(1.06, 1), align = "v")

ggsave(plot_oc_ac, file = "results/supplement_Fig.S_interaction_subsetSmallRegion.png", 
       width = 150, height = 120, units = 'mm', dpi = 300)



######################
## sensitivity analyses based on data sampled in the first and last year
## figure showing the relationship between occupancy/abundance change and thermal position 

rm(list = ls())
gc()

load("results/data_input_to_models.RDATA") # observation data
brm_oc_tmeanpos_2yr <- readRDS("models/brm_oc_tmeanpos_2yr.rds")
brm_ac_tmeanpos_2yr <- readRDS("models/brm_ac_tmeanpos_2yr.rds")

## new data for predicting relations across realms
position_seq <- seq(0, 1, length = 101) # sequence of thermal position

newdata <- oc_2yr %>% 
  distinct(study, realm, n_samp, duration, extent_km2) %>% group_by(realm) %>% 
  summarise(n_samp = round(exp(mean(log(n_samp)))), # use the geometric mean 
            duration = round(exp(mean(log(duration)))),
            extent_km2 = round(exp(mean(log(extent_km2))), 1)) %>% 
  left_join(tibble(tmeanpos = rep(position_seq, 3), 
                   realm = rep(c("Terrestrial", "Freshwater", "Marine"), each = length(position_seq))))

## text showing fixed effect
beta_oc_tmeanpos <- fixef(brm_oc_tmeanpos_2yr, robust = TRUE, probs = c(0.025, 0.975)) %>%
  as.data.frame() %>% 
  rownames_to_column(var = "term") %>%
  filter(grepl(":tmeanpos", term)) %>%
  mutate(realm = gsub(":tmeanpos", "", term),
         realm = gsub("realm", "", realm),
         realm = factor(realm, levels = c("Terrestrial", "Freshwater", "Marine"))) %>%
  mutate(label = paste('beta == ', round(Estimate, 3),'~(', round(Q2.5, 3), "*','~", round(Q97.5, 3),')'))

beta_ac_tmeanpos <- fixef(brm_ac_tmeanpos_2yr, robust = TRUE, probs = c(0.025, 0.975)) %>%
  as.data.frame() %>% 
  rownames_to_column(var = "term") %>%
  filter(grepl(":tmeanpos", term)) %>%
  mutate(realm = gsub(":tmeanpos", "", term),
         realm = gsub("realm", "", realm),
         realm = factor(realm, levels = c("Terrestrial", "Freshwater", "Marine"))) %>%
  mutate(label = paste('beta == ', round(Estimate, 3),'~(', round(Q2.5, 3), "*','~", round(Q97.5, 3),')'))

## fitted values along thermal position
fitted_oc_tmeanpos <- fitted(brm_oc_tmeanpos_2yr, newdata = newdata, re_formula = NA, ndraws = 1000, probs = c(0.025, 0.1, 0.9, 0.975)) %>%
  as_tibble()  %>% 
  bind_cols(newdata) %>% 
  mutate(realm = factor(realm, levels = c("Terrestrial", "Freshwater", "Marine")))

fitted_ac_tmeanpos <- fitted(brm_ac_tmeanpos_2yr, newdata = newdata, re_formula = NA, ndraws = 1000, probs = c(0.025, 0.1, 0.9, 0.975)) %>%
  as_tibble()  %>% 
  bind_cols(newdata) %>% 
  mutate(realm = factor(realm, levels = c("Terrestrial", "Freshwater", "Marine")))


## draw figures
# colors for realm
realm_col <- c("Terrestrial" = "#D95F02", "Freshwater" = "#44bbff", "Marine" = "#7570B3")
realm_col <- c("All" = "#000000", "Terrestrial" = "#A65628", "Freshwater" = "#33A02C", "Marine" = "#1F78B4")

plot_oc <- ggplot(fitted_oc_tmeanpos) +
  facet_wrap(~ realm, nrow = 1) + 
  geom_hline(yintercept = 0, lty = 2, linewidth = 0.3) + 
  geom_ribbon(aes(x =  tmeanpos, ymin = Q2.5, ymax = Q97.5, fill = realm), alpha = 0.3) + 
  geom_line(aes(x = tmeanpos, y = Estimate, color = realm), lwd = 1) + 
  geom_text(data = beta_oc_tmeanpos, mapping = aes(x= 0.12, y = Inf, label = label), 
            hjust = 0, vjust = 1.5,  parse = T, size = 2.46) + 
  labs(x = "Baseline thermal position", y = "Occupancy change") + 
  scale_colour_manual(name = NULL, values = realm_col) + 
  scale_fill_manual(name = NULL, values = realm_col) + 
  theme_bw()  +
  theme(legend.position = "no", 
        text = element_text(size =9),
        strip.text = element_text(size = 9),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 9), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

plot_ac <- ggplot(fitted_ac_tmeanpos) +
  facet_wrap(~ realm, nrow = 1) + 
  geom_hline(yintercept = 0, lty = 2, linewidth = 0.3) + 
  geom_ribbon(aes(x =  tmeanpos, ymin = Q2.5, ymax = Q97.5, fill = realm), alpha = 0.3) + 
  geom_line(aes(x = tmeanpos, y = Estimate,  color = realm), lwd = 1) + 
  geom_text(data = beta_ac_tmeanpos, mapping = aes(x= 0.12, y = Inf, label = label), 
            hjust = 0, vjust = 1.5,  parse = T, size = 2.46) + 
  # scale_y_continuous(limits = c(-2.5,2.5)) +
  labs(x = "Baseline thermal position", y = "Abundance change") + 
  scale_colour_manual(name = NULL, values = realm_col) + 
  scale_fill_manual(name = NULL, values = realm_col) + 
  theme_bw()  +
  theme(legend.position = "no", 
        text = element_text(size =9),
        strip.text = element_text(size = 9),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 9), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

plot_oc_ac <- cowplot::plot_grid(plot_oc, plot_ac, nrow = 2, labels = c("a", "b"), label_size = 10)

ggsave(plot_oc_ac, file = "results/Supplement_Fig.S_effect_thermal_position_first-lastYearData.png", 
       width = 140, height = 110, units = 'mm', dpi = 600) 



## figure showing the interaction between baseline thermal position and temperature change
rm(list = ls())
gc()

load("results/data_input_to_models.RDATA") # observation data
brm_oc_tmeanpos_tmeanchange_2yr <- readRDS("models/brm_oc_tmeanpos_tmeanchange_2yr.rds")
brm_ac_tmeanpos_tmeanchange_2yr <- readRDS("models/brm_ac_tmeanpos_tmeanchange_2yr.rds")

## prepare new data to predict occupancy or abundance change
# sequence of thermal position
position_seq <- seq(0, 1, length = 101)

# temperature change 
tchange_ranges <- dat_meta %>%
  summarise(Mean = mean(change_tempmean_2yr/(duration-1)),
            sd = sd(change_tempmean_2yr/(duration-1)),
            low = round((Mean - 2*sd), 3),
            high = round((Mean + 2*sd), 3)) %>%
  dplyr::select(-sd)

tchange_seq <- seq(tchange_ranges$low, tchange_ranges$high, length = 101)
tchange_realm <- tibble(expand.grid(change_tempmean = tchange_seq, realm = c("Terrestrial", "Freshwater", "Marine")))

# new data for prediction
newdata <- tchange_realm %>% 
  left_join(expand.grid(change_tempmean = unique(tchange_realm$change_tempmean), tmeanpos = position_seq), relationship = "many-to-many") %>%
  left_join(oc_2yr %>% 
              distinct(study, realm, n_samp, duration, extent_km2) %>% group_by(realm) %>% 
              summarise(n_samp = round(exp(mean(log(n_samp)))), # use the geometric mean 
                        duration = round(exp(mean(log(duration)))),
                        extent_km2 = round(exp(mean(log(extent_km2))), 1))) %>% 
  mutate(realm = factor(realm, levels = c("Terrestrial", "Freshwater", "Marine")))

## predicted occupancy and abundance change 
fitted_oc_tmeanpos <- fitted(brm_oc_tmeanpos_tmeanchange_2yr, newdata = newdata, re_formula = NA, ndraws = 1000, probs = c(0.025, 0.975)) %>%
  as_tibble()  %>% 
  bind_cols(newdata)

fitted_ac_tmeanpos <- fitted(brm_ac_tmeanpos_tmeanchange_2yr, newdata = newdata, re_formula = NA, ndraws = 1000, probs = c(0.025, 0.975)) %>%
  as_tibble()  %>% 
  bind_cols(newdata)

# text showing interaction effect
beta_oc_interaction_text <- fixef(brm_oc_tmeanpos_tmeanchange_2yr) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "term") %>%
  filter(grepl("tmeanpos:change_tempmean", term)) %>%
  mutate(realm = gsub(":tmeanpos:change_tempmean", "", term),
         realm = gsub("realm", "", realm),
         realm = factor(realm, levels = c("Terrestrial", "Freshwater", "Marine"))) %>%
  mutate(label = paste('beta == ', round(Estimate, 2),'~(', round(Q2.5, 2), "*','~", round(Q97.5, 2),')'))

beta_ac_interaction_text <- fixef(brm_ac_tmeanpos_tmeanchange_2yr) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "term") %>%
  filter(grepl("tmeanpos:change_tempmean", term)) %>%
  mutate(realm = gsub(":tmeanpos:change_tempmean", "", term),
         realm = gsub("realm", "", realm),
         realm = factor(realm, levels = c("Terrestrial", "Freshwater", "Marine"))) %>%
  mutate(label = paste('beta == ', round(Estimate, 2),'~(', round(Q2.5, 2), "*','~", round(Q97.5, 2),')'))


## draw figures
plot_oc <- ggplot(fitted_oc_tmeanpos) + 
  facet_wrap(~ realm, scales = "fixed") +
  geom_line(aes(x = tmeanpos, y = Estimate, color = change_tempmean, group = change_tempmean)) +
  geom_hline(yintercept = 0, lty =2, lwd = 0.4) + 
  geom_text(data = beta_oc_interaction_text, mapping = aes(x= 0.1, y = Inf, label = label), 
            hjust = 0, vjust = 1.5,  parse = T, size = 2.46) + 
  labs(tag = "a", x = "Baseline thermal position", y = "Occupancy change", 
       color = expression('Temperature change ('*degree*'C/year)')) + 
  scale_color_gradient2(low = "blue", mid = "darkgray", high="red") +
  theme_bw() + 
  theme(legend.position = c(0.7, 1.3), 
        plot.margin = unit(c(0.6,0.2,0,0.2), "cm"),
        legend.direction = 'horizontal',
        legend.key.size = unit(4, 'mm'),
        legend.background = element_blank(),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 9),
        text = element_text(size = 10),
        strip.text = element_text(size = 9),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  guides(color = guide_colourbar(barwidth = 4, barheight = 0.5))

plot_ac <- ggplot(fitted_ac_tmeanpos) + 
  facet_wrap(~ realm, scales = "fixed") +
  geom_line(aes(x = tmeanpos, y = Estimate, color = change_tempmean, group = change_tempmean)) +
  geom_hline(yintercept = 0, lty =2, lwd = 0.4) + 
  geom_text(data = beta_ac_interaction_text, mapping = aes(x= 0.1, y = Inf, label = label), 
            hjust = 0, vjust = 1.5,  parse = T, size = 2.46) + 
  labs(tag = "b", x = "Baseline thermal position", y = "Abundance change") + 
  scale_color_gradient2(low = "blue", mid = "darkgray", high="red") +
  theme_bw() + 
  theme(legend.position = "no",
        plot.margin = unit(c(0,0.2,0.2,0.2), "cm"),
        text = element_text(size =10),
        strip.text = element_text(size = 9),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

plot_oc_ac <- cowplot::plot_grid(plot_oc, plot_ac,
                                 nrow = 2, rel_heights = c(1.06, 1), align = "v")

ggsave(plot_oc_ac, file = "results/supplement_Fig.S_interaction_first-lastYearData.png", 
       width = 150, height = 120, units = 'mm', dpi = 300)



######################
## sensitivity analyses based on thermal  position estimated as 1% and 99% quantile of observed temperatures
## figure showing the relationship between occupancy/abundance change and thermal position 

rm(list = ls())
gc()

load("results/data_input_to_models.RDATA") # observation data
brm_oc_tmeanpos98_realm <- readRDS("models/brm_oc_tmeanpos98_realm.rds")
brm_ac_tmeanpos98_realm <- readRDS("models/brm_ac_tmeanpos98_realm.rds")

## new data for predicting relations across realms
position_seq <- seq(0, 1, length = 101) # sequence of thermal position

newdata <- oc_period %>% 
  distinct(study, realm, n_samp, duration, extent_km2) %>% group_by(realm) %>% 
  summarise(n_samp = round(exp(mean(log(n_samp)))), # use the geometric mean 
            duration = round(exp(mean(log(duration)))),
            extent_km2 = round(exp(mean(log(extent_km2))), 1)) %>% 
  left_join(tibble(tmeanpos98 = rep(position_seq, 3), 
                   realm = rep(c("Terrestrial", "Freshwater", "Marine"), each = length(position_seq))))

## text showing fixed effect
beta_oc_tmeanpos <- fixef(brm_oc_tmeanpos98_realm, robust = TRUE, probs = c(0.025, 0.975)) %>%
  as.data.frame() %>% 
  rownames_to_column(var = "term") %>%
  filter(grepl(":tmeanpos", term)) %>%
  mutate(realm = gsub(":tmeanpos98", "", term),
         realm = gsub("realm", "", realm),
         realm = factor(realm, levels = c("Terrestrial", "Freshwater", "Marine"))) %>%
  mutate(label = paste('beta == ', round(Estimate, 3),'~(', round(Q2.5, 3), "*','~", round(Q97.5, 3),')'))

beta_ac_tmeanpos <- fixef(brm_ac_tmeanpos98_realm, robust = TRUE, probs = c(0.025, 0.975)) %>%
  as.data.frame() %>% 
  rownames_to_column(var = "term") %>%
  filter(grepl(":tmeanpos", term)) %>%
  mutate(realm = gsub(":tmeanpos98", "", term),
         realm = gsub("realm", "", realm),
         realm = factor(realm, levels = c("Terrestrial", "Freshwater", "Marine"))) %>%
  mutate(label = paste('beta == ', round(Estimate, 3),'~(', round(Q2.5, 3), "*','~", round(Q97.5, 3),')'))

## fitted values along thermal position
fitted_oc_tmeanpos <- fitted(brm_oc_tmeanpos98_realm, newdata = newdata, re_formula = NA, ndraws = 1000, probs = c(0.025, 0.1, 0.9, 0.975)) %>%
  as_tibble()  %>% 
  bind_cols(newdata) %>% 
  mutate(realm = factor(realm, levels = c("Terrestrial", "Freshwater", "Marine")))

fitted_ac_tmeanpos <- fitted(brm_ac_tmeanpos98_realm, newdata = newdata, re_formula = NA, ndraws = 1000, probs = c(0.025, 0.1, 0.9, 0.975)) %>%
  as_tibble()  %>% 
  bind_cols(newdata) %>% 
  mutate(realm = factor(realm, levels = c("Terrestrial", "Freshwater", "Marine")))


## draw figures
# colors for realm
realm_col <- c("Terrestrial" = "#D95F02", "Freshwater" = "#44bbff", "Marine" = "#7570B3")
realm_col <- c("All" = "#000000", "Terrestrial" = "#A65628", "Freshwater" = "#33A02C", "Marine" = "#1F78B4")

plot_oc <- ggplot(fitted_oc_tmeanpos) +
  facet_wrap(~ realm, nrow = 1) + 
  geom_hline(yintercept = 0, lty = 2, linewidth = 0.3) + 
  geom_ribbon(aes(x =  tmeanpos98, ymin = Q2.5, ymax = Q97.5, fill = realm), alpha = 0.3) + 
  geom_line(aes(x = tmeanpos98, y = Estimate, color = realm), lwd = 1) + 
  geom_text(data = beta_oc_tmeanpos, mapping = aes(x= 0.12, y = Inf, label = label), 
            hjust = 0, vjust = 1.5,  parse = T, size = 2.46) + 
  labs(x = "Baseline thermal position", y = "Occupancy change") + 
  scale_colour_manual(name = NULL, values = realm_col) + 
  scale_fill_manual(name = NULL, values = realm_col) + 
  theme_bw()  +
  theme(legend.position = "no", 
        text = element_text(size =9),
        strip.text = element_text(size = 9),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 9), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

plot_ac <- ggplot(fitted_ac_tmeanpos) +
  facet_wrap(~ realm, nrow = 1) + 
  geom_hline(yintercept = 0, lty = 2, linewidth = 0.3) + 
  geom_ribbon(aes(x =  tmeanpos98, ymin = Q2.5, ymax = Q97.5, fill = realm), alpha = 0.3) + 
  geom_line(aes(x = tmeanpos98, y = Estimate,  color = realm), lwd = 1) + 
  geom_text(data = beta_ac_tmeanpos, mapping = aes(x= 0.12, y = Inf, label = label), 
            hjust = 0, vjust = 1.5,  parse = T, size = 2.46) + 
  # scale_y_continuous(limits = c(-2.5,2.5)) +
  labs(x = "Baseline thermal position", y = "Abundance change") + 
  scale_colour_manual(name = NULL, values = realm_col) + 
  scale_fill_manual(name = NULL, values = realm_col) + 
  theme_bw()  +
  theme(legend.position = "no", 
        text = element_text(size =9),
        strip.text = element_text(size = 9),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 9), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

plot_oc_ac <- cowplot::plot_grid(plot_oc, plot_ac, nrow = 2, labels = c("a", "b"), label_size = 10)

ggsave(plot_oc_ac, file = "results/Supplement_Fig.S_effect_thermal_position_2.5-97.5quantiles.png", 
       width = 140, height = 110, units = 'mm', dpi = 600) 



## figure showing the interaction between baseline thermal position and temperature change
rm(list = ls())
gc()

load("results/data_input_to_models.RDATA") # observation data
brm_oc_tmeanpos98_tmeanchange_realm <- readRDS("models/brm_oc_tmeanpos98_tmeanchange_realm.rds")
brm_ac_tmeanpos98_tmeanchange_realm <- readRDS("models/brm_ac_tmeanpos98_tmeanchange_realm.rds")

## prepare new data to predict occupancy or abundance change
# sequence of thermal position
position_seq <- seq(0, 1, length = 101)

# temperature change 
tchange_ranges <- dat_meta %>%
  summarise(Mean = mean(change_tempmean/(duration_mean-1)),
            sd = sd(change_tempmean/(duration_mean-1)),
            low = round((Mean - 2*sd), 3),
            high = round((Mean + 2*sd), 3)) %>%
  dplyr::select(-sd)

tchange_seq <- seq(tchange_ranges$low, tchange_ranges$high, length = 101)
tchange_realm <- tibble(expand.grid(change_tempmean = tchange_seq, realm = c("Terrestrial", "Freshwater", "Marine")))

# new data for prediction
newdata <- tchange_realm %>% 
  left_join(expand.grid(change_tempmean = unique(tchange_realm$change_tempmean), tmeanpos98 = position_seq), relationship = "many-to-many") %>%
  left_join(oc_period %>% 
              distinct(study, realm, n_samp, duration, extent_km2) %>% group_by(realm) %>% 
              summarise(n_samp = round(exp(mean(log(n_samp)))), # use the geometric mean 
                        duration = round(exp(mean(log(duration)))),
                        extent_km2 = round(exp(mean(log(extent_km2))), 1))) %>% 
  mutate(realm = factor(realm, levels = c("Terrestrial", "Freshwater", "Marine")))

## predicted occupancy and abundance change 
fitted_oc_tmeanpos <- fitted(brm_oc_tmeanpos98_tmeanchange_realm, newdata = newdata, re_formula = NA, ndraws = 1000, probs = c(0.025, 0.975)) %>%
  as_tibble()  %>% 
  bind_cols(newdata)

fitted_ac_tmeanpos <- fitted(brm_ac_tmeanpos98_tmeanchange_realm, newdata = newdata, re_formula = NA, ndraws = 1000, probs = c(0.025, 0.975)) %>%
  as_tibble()  %>% 
  bind_cols(newdata)

# text showing interaction effect
beta_oc_interaction_text <- fixef(brm_oc_tmeanpos98_tmeanchange_realm) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "term") %>%
  filter(grepl("tmeanpos98:change_tempmean", term)) %>%
  mutate(realm = gsub(":tmeanpos98:change_tempmean", "", term),
         realm = gsub("realm", "", realm),
         realm = factor(realm, levels = c("Terrestrial", "Freshwater", "Marine"))) %>%
  mutate(label = paste('beta == ', round(Estimate, 2),'~(', round(Q2.5, 2), "*','~", round(Q97.5, 2),')'))

beta_ac_interaction_text <- fixef(brm_ac_tmeanpos98_tmeanchange_realm) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "term") %>%
  filter(grepl("tmeanpos98:change_tempmean", term)) %>%
  mutate(realm = gsub(":tmeanpos98:change_tempmean", "", term),
         realm = gsub("realm", "", realm),
         realm = factor(realm, levels = c("Terrestrial", "Freshwater", "Marine"))) %>%
  mutate(label = paste('beta == ', round(Estimate, 2),'~(', round(Q2.5, 2), "*','~", round(Q97.5, 2),')'))


## draw figures
plot_oc <- ggplot(fitted_oc_tmeanpos) + 
  facet_wrap(~ realm, scales = "fixed") +
  geom_line(aes(x = tmeanpos98, y = Estimate, color = change_tempmean, group = change_tempmean)) +
  geom_hline(yintercept = 0, lty =2, lwd = 0.4) + 
  geom_text(data = beta_oc_interaction_text, mapping = aes(x= 0.1, y = Inf, label = label), 
            hjust = 0, vjust = 1.5,  parse = T, size = 2.46) + 
  labs(tag = "a", x = "Baseline thermal position", y = "Occupancy change", 
       color = expression('Temperature change ('*degree*'C/year)')) + 
  scale_color_gradient2(low = "blue", mid = "darkgray", high="red") +
  theme_bw() + 
  theme(legend.position = c(0.7, 1.3), 
        plot.margin = unit(c(0.6,0.2,0,0.2), "cm"),
        legend.direction = 'horizontal',
        legend.key.size = unit(4, 'mm'),
        legend.background = element_blank(),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 9),
        text = element_text(size = 10),
        strip.text = element_text(size = 9),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  guides(color = guide_colourbar(barwidth = 4, barheight = 0.5))

plot_ac <- ggplot(fitted_ac_tmeanpos) + 
  facet_wrap(~ realm, scales = "fixed") +
  geom_line(aes(x = tmeanpos98, y = Estimate, color = change_tempmean, group = change_tempmean)) +
  geom_hline(yintercept = 0, lty =2, lwd = 0.4) + 
  geom_text(data = beta_ac_interaction_text, mapping = aes(x= 0.1, y = Inf, label = label), 
            hjust = 0, vjust = 1.5,  parse = T, size = 2.46) + 
  labs(tag = "b", x = "Baseline thermal position", y = "Abundance change") + 
  scale_color_gradient2(low = "blue", mid = "darkgray", high="red") +
  theme_bw() + 
  theme(legend.position = "no",
        plot.margin = unit(c(0,0.2,0.2,0.2), "cm"),
        text = element_text(size =10),
        strip.text = element_text(size = 9),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

plot_oc_ac <- cowplot::plot_grid(plot_oc, plot_ac,
                                 nrow = 2, rel_heights = c(1.06, 1), align = "v")

ggsave(plot_oc_ac, file = "results/supplement_Fig.S_interaction_2.5-97.5quantiles.png", 
       width = 150, height = 120, units = 'mm', dpi = 300)




######################
## sensitivity analyses based on the temperature at the sampling years
## figure showing the relationship between occupancy/abundance change and thermal position 

rm(list = ls())
gc()

load("results/data_input_to_models.RDATA") # observation data
brm_oc_tmeanpos_sampyr <- readRDS("models/brm_oc_tmeanpos_sampyr.rds")
brm_ac_tmeanpos_sampyr <- readRDS("models/brm_ac_tmeanpos_sampyr.rds")


## new data for predicting relations across realms
position_seq <- seq(0, 1, length = 101) # sequence of thermal position

newdata <- oc_period %>% distinct(study, realm, n_samp, duration, extent_km2) %>% group_by(realm) %>% 
  summarise(n_samp = round(exp(mean(log(n_samp)))), # use the geometric mean 
            duration = round(exp(mean(log(duration)))),
            extent_km2 = round(exp(mean(log(extent_km2))), 1)) %>% 
  left_join(tibble(tmeanpos_sampyr = rep(position_seq, 3), 
                   realm = rep(c("Terrestrial", "Freshwater", "Marine"), each = length(position_seq))))

## text showing fixed effect
beta_oc_tmeanpos <- fixef(brm_oc_tmeanpos_sampyr, robust = TRUE, probs = c(0.025, 0.975)) %>%
  as.data.frame() %>% 
  rownames_to_column(var = "term") %>%
  filter(grepl(":tmeanpos_sampyr", term)) %>%
  mutate(realm = gsub(":tmeanpos_sampyr", "", term),
         realm = gsub("realm", "", realm),
         realm = factor(realm, levels = c("Terrestrial", "Freshwater", "Marine"))) %>%
  mutate(label = paste('beta == ', round(Estimate, 3),'~(', round(Q2.5, 3), "*','~", round(Q97.5, 3),')'))

beta_ac_tmeanpos <- fixef(brm_ac_tmeanpos_sampyr, robust = TRUE, probs = c(0.025, 0.975)) %>%
  as.data.frame() %>% 
  rownames_to_column(var = "term") %>%
  filter(grepl(":tmeanpos_sampyr", term)) %>%
  mutate(realm = gsub(":tmeanpos_sampyr", "", term),
         realm = gsub("realm", "", realm),
         realm = factor(realm, levels = c("Terrestrial", "Freshwater", "Marine"))) %>%
  mutate(label = paste('beta == ', round(Estimate, 3),'~(', round(Q2.5, 3), "*','~", round(Q97.5, 3),')'))

## fitted values along thermal position
fitted_oc_tmeanpos <- fitted(brm_oc_tmeanpos_sampyr, newdata = newdata, re_formula = NA, ndraws = 1000, probs = c(0.025, 0.1, 0.9, 0.975)) %>%
  as_tibble()  %>% 
  bind_cols(newdata) %>% 
  mutate(realm = factor(realm, levels = c("Terrestrial", "Freshwater", "Marine")))

fitted_ac_tmeanpos <- fitted(brm_ac_tmeanpos_sampyr, newdata = newdata, re_formula = NA, ndraws = 1000, probs = c(0.025, 0.1, 0.9, 0.975)) %>%
  as_tibble()  %>% 
  bind_cols(newdata) %>% 
  mutate(realm = factor(realm, levels = c("Terrestrial", "Freshwater", "Marine")))


## draw figures
# colors for realm
realm_col <- c("Terrestrial" = "#D95F02", "Freshwater" = "#44bbff", "Marine" = "#7570B3")
realm_col <- c("All" = "#000000", "Terrestrial" = "#A65628", "Freshwater" = "#33A02C", "Marine" = "#1F78B4")

plot_oc <- ggplot(fitted_oc_tmeanpos) +
  facet_wrap(~ realm, nrow = 1) + 
  geom_hline(yintercept = 0, lty = 2, linewidth = 0.3) + 
  geom_ribbon(aes(x =  tmeanpos_sampyr, ymin = Q2.5, ymax = Q97.5, fill = realm), alpha = 0.3) + 
  geom_line(aes(x = tmeanpos_sampyr, y = Estimate, color = realm), lwd = 1) + 
  geom_text(data = beta_oc_tmeanpos, mapping = aes(x= 0.12, y = Inf, label = label), 
            hjust = 0, vjust = 1.5,  parse = T, size = 2.46) + 
  labs(x = "Baseline thermal position", y = "Occupancy change") + 
  scale_colour_manual(name = NULL, values = realm_col) + 
  scale_fill_manual(name = NULL, values = realm_col) + 
  theme_bw()  +
  theme(legend.position = "no", 
        text = element_text(size =9),
        strip.text = element_text(size = 9),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 9), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

plot_ac <- ggplot(fitted_ac_tmeanpos) +
  facet_wrap(~ realm, nrow = 1) + 
  geom_hline(yintercept = 0, lty = 2, linewidth = 0.3) + 
  geom_ribbon(aes(x =  tmeanpos_sampyr, ymin = Q2.5, ymax = Q97.5, fill = realm), alpha = 0.3) + 
  geom_line(aes(x = tmeanpos_sampyr, y = Estimate,  color = realm), lwd = 1) + 
  geom_text(data = beta_ac_tmeanpos, mapping = aes(x= 0.12, y = Inf, label = label), 
            hjust = 0, vjust = 1.5,  parse = T, size = 2.46) + 
  # scale_y_continuous(limits = c(-2.5,2.5)) +
  labs(x = "Baseline thermal position", y = "Abundance change") + 
  scale_colour_manual(name = NULL, values = realm_col) + 
  scale_fill_manual(name = NULL, values = realm_col) + 
  theme_bw()  +
  theme(legend.position = "no", 
        text = element_text(size =9),
        strip.text = element_text(size = 9),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 9), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

plot_oc_ac <- cowplot::plot_grid(plot_oc, plot_ac, nrow = 2, labels = c("a", "b"), label_size = 10)

ggsave(plot_oc_ac, file = "results/Supplement_Fig.S_effect_thermal_position_sampyrTemperature.png", 
       width = 140, height = 110, units = 'mm', dpi = 600) 



## figure showing the interaction between baseline thermal position and temperature change
rm(list = ls())
gc()

load("results/data_input_to_models.RDATA") # observation data
brm_oc_tmeanpos_tmeanchange_sampyr <- readRDS("models/brm_oc_tmeanpos_tmeanchange_sampyr.rds")
brm_ac_tmeanpos_tmeanchange_sampyr <- readRDS("models/brm_ac_tmeanpos_tmeanchange_sampyr.rds")

## prepare new data to predict occupancy or abundance change
# sequence of thermal position
position_seq <- seq(0, 1, length = 101)

# temperature change 
tchange_ranges <- dat_meta %>%
  summarise(Mean = mean(change_tempmean_sampyr/(duration_mean-1)),
            sd = sd(change_tempmean_sampyr/(duration_mean-1)),
            low = round((Mean - 2*sd), 3),
            high = round((Mean + 2*sd), 3)) %>%
  dplyr::select(-sd)

tchange_seq <- seq(tchange_ranges$low, tchange_ranges$high, length = 101)
tchange_realm <- tibble(expand.grid(change_tempmean_sampyr = tchange_seq, realm = c("Terrestrial", "Freshwater", "Marine")))

# new data for prediction
newdata <- tchange_realm %>% 
  left_join(expand.grid(change_tempmean_sampyr = unique(tchange_realm$change_tempmean_sampyr), tmeanpos_sampyr = position_seq), relationship = "many-to-many") %>%
  left_join(oc_period %>% distinct(study, realm, n_samp, duration, extent_km2) %>% group_by(realm) %>% 
              summarise(n_samp = round(exp(mean(log(n_samp)))), # use the geometric mean 
                        duration = round(exp(mean(log(duration)))),
                        extent_km2 = round(exp(mean(log(extent_km2))), 1))) %>% 
  mutate(realm = factor(realm, levels = c("Terrestrial", "Freshwater", "Marine")))

## predicted occupancy and abundance change 
fitted_oc_tmeanpos <- fitted(brm_oc_tmeanpos_tmeanchange_sampyr, newdata = newdata, re_formula = NA, ndraws = 1000, probs = c(0.025, 0.975)) %>%
  as_tibble()  %>% 
  bind_cols(newdata)

fitted_ac_tmeanpos <- fitted(brm_ac_tmeanpos_tmeanchange_sampyr, newdata = newdata, re_formula = NA, ndraws = 1000, probs = c(0.025, 0.975)) %>%
  as_tibble()  %>% 
  bind_cols(newdata)

# text showing interaction effect
beta_oc_interaction_text <- fixef(brm_oc_tmeanpos_tmeanchange_sampyr) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "term") %>%
  filter(grepl("tmeanpos_sampyr:change_tempmean_sampyr", term)) %>%
  mutate(realm = gsub(":tmeanpos_sampyr:change_tempmean_sampyr", "", term),
         realm = gsub("realm", "", realm),
         realm = factor(realm, levels = c("Terrestrial", "Freshwater", "Marine"))) %>%
  mutate(label = paste('beta == ', round(Estimate, 2),'~(', round(Q2.5, 2), "*','~", round(Q97.5, 2),')'))

beta_ac_interaction_text <- fixef(brm_ac_tmeanpos_tmeanchange_sampyr) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "term") %>%
  filter(grepl("tmeanpos_sampyr:change_tempmean_sampyr", term)) %>%
  mutate(realm = gsub(":tmeanpos_sampyr:change_tempmean_sampyr", "", term),
         realm = gsub("realm", "", realm),
         realm = factor(realm, levels = c("Terrestrial", "Freshwater", "Marine"))) %>%
  mutate(label = paste('beta == ', round(Estimate, 2),'~(', round(Q2.5, 2), "*','~", round(Q97.5, 2),')'))


## draw figures
plot_oc <- ggplot(fitted_oc_tmeanpos) + 
  facet_wrap(~ realm, scales = "fixed") +
  geom_line(aes(x = tmeanpos_sampyr, y = Estimate, color = change_tempmean_sampyr, group = change_tempmean_sampyr)) +
  geom_hline(yintercept = 0, lty =2, lwd = 0.4) + 
  geom_text(data = beta_oc_interaction_text, mapping = aes(x= 0.1, y = Inf, label = label), 
            hjust = 0, vjust = 1.5,  parse = T, size = 2.46) + 
  labs(tag = "a", x = "Baseline thermal position", y = "Occupancy change", 
       color = expression('Temperature change ('*degree*'C/year)')) + 
  scale_color_gradient2(low = "blue", mid = "darkgray", high="red") +
  theme_bw() + 
  theme(legend.position = c(0.7, 1.3), 
        plot.margin = unit(c(0.6,0.2,0,0.2), "cm"),
        legend.direction = 'horizontal',
        legend.key.size = unit(4, 'mm'),
        legend.background = element_blank(),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 9),
        text = element_text(size = 10),
        strip.text = element_text(size = 9),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  guides(color = guide_colourbar(barwidth = 4, barheight = 0.5))

plot_ac <- ggplot(fitted_ac_tmeanpos) + 
  facet_wrap(~ realm, scales = "fixed") +
  geom_line(aes(x = tmeanpos_sampyr, y = Estimate, color = change_tempmean_sampyr, group = change_tempmean_sampyr)) +
  geom_hline(yintercept = 0, lty =2, lwd = 0.4) + 
  geom_text(data = beta_ac_interaction_text, mapping = aes(x= 0.1, y = Inf, label = label), 
            hjust = 0, vjust = 1.5,  parse = T, size = 2.46) + 
  labs(tag = "b", x = "Baseline thermal position", y = "Abundance change") + 
  scale_color_gradient2(low = "blue", mid = "darkgray", high="red") +
  theme_bw() + 
  theme(legend.position = "no",
        plot.margin = unit(c(0,0.2,0.2,0.2), "cm"),
        text = element_text(size =10),
        strip.text = element_text(size = 9),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

plot_oc_ac <- cowplot::plot_grid(plot_oc, plot_ac,
                   nrow = 2, rel_heights = c(1.06, 1), align = "v")

ggsave(plot_oc_ac, file = "results/supplement_Fig.S_interaction_sampyrTemperature.png", 
       width = 150, height = 120, units = 'mm', dpi = 300)



######################
##  fixed effects of baseline thermal position and interactions with temperature change across realms
rm(list = ls())
gc()

brm_oc_tmeanpos_tmeanchange_realm <- readRDS("models/brm_oc_tmeanpos_tmeanchange_realm.rds")
brm_ac_tmeanpos_tmeanchange_realm <- readRDS("models/brm_ac_tmeanpos_tmeanchange_realm.rds")

# fixed effects on occupancy change
fixef_oc_tmeanpos_tmeanchange_realm <- fixef(brm_oc_tmeanpos_tmeanchange_realm, probs = c(0.025, 0.975, 0.05, 0.95, 0.1, 0.9)) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "term") %>%
  filter(!grepl("sigma", term)) %>%
  mutate(term =  gsub(pattern = "realm", "", .[, "term"])) %>%
  mutate(term = gsub(pattern = "tmeanpos:change_tempmean", "Interaction", .[, "term"])) %>%
  separate(term, c( "realm", "predictor"), sep = ":") %>% 
  mutate(predictor = ifelse(is.na(predictor), "Intercept", predictor)) %>%
  mutate(Estimate = ifelse(predictor %in% c("Intercept", "tmeanpos"), 10*Estimate, Estimate),
         Q2.5 = ifelse(predictor %in% c("Intercept", "tmeanpos"), 10*Q2.5, Q2.5),
         Q97.5 = ifelse(predictor %in% c("Intercept", "tmeanpos"), 10*Q97.5, Q97.5),
         Q5 = ifelse(predictor %in% c("Intercept", "tmeanpos"), 10*Q5, Q5),
         Q95 = ifelse(predictor %in% c("Intercept", "tmeanpos"), 10*Q95, Q95),
         Q10 = ifelse(predictor %in% c("Intercept", "tmeanpos"), 10*Q10, Q10),
         Q90 = ifelse(predictor %in% c("Intercept", "tmeanpos"), 10*Q90, Q90)) %>%        
  mutate(predictor = factor(predictor, levels = c("Intercept", "tmeanpos", "change_tempmean", "Interaction"), 
                            labels = c("10*Intercept", "10*Thermal\nposition", "Temperature\nchange", "Interaction")),
         realm = factor(realm, levels = c("Terrestrial", "Freshwater", "Marine")))

# fixed effects on abundance change
fixef_ac_tmeanpos_tmeanchange_realm <- fixef(brm_ac_tmeanpos_tmeanchange_realm, probs = c(0.025, 0.975, 0.05, 0.95, 0.1, 0.9)) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "term") %>%
  filter(!grepl("sigma", term)) %>%
  mutate(term =  gsub(pattern = "realm", "", .[, "term"])) %>%
  mutate(term = gsub(pattern = "tmeanpos:change_tempmean", "Interaction", .[, "term"])) %>%
  separate(term, c( "realm", "predictor"), sep = ":") %>% 
  mutate(predictor = ifelse(is.na(predictor), "Intercept", predictor)) %>%
  mutate(Estimate = ifelse(predictor %in% c("Intercept", "tmeanpos"), 10*Estimate, Estimate),
         Q2.5 = ifelse(predictor %in% c("Intercept", "tmeanpos"), 10*Q2.5, Q2.5),
         Q97.5 = ifelse(predictor %in% c("Intercept", "tmeanpos"), 10*Q97.5, Q97.5),
         Q5 = ifelse(predictor %in% c("Intercept", "tmeanpos"), 10*Q5, Q5),
         Q95 = ifelse(predictor %in% c("Intercept", "tmeanpos"), 10*Q95, Q95),
         Q10 = ifelse(predictor %in% c("Intercept", "tmeanpos"), 10*Q10, Q10),
         Q90 = ifelse(predictor %in% c("Intercept", "tmeanpos"), 10*Q90, Q90)) %>%       
  mutate(predictor = factor(predictor, levels = c("Intercept", "tmeanpos", "change_tempmean", "Interaction"), 
                            labels = c("10*Intercept", "10*Thermal\nposition", "Temperature\nchange", "Interaction")),
         realm = factor(realm, levels = c("Terrestrial", "Freshwater", "Marine")))


# generate figures
library(ggstance)

# colors for realm
realm_col <- c("Terrestrial" = "#D95F02", "Freshwater" = "#44bbff", "Marine" = "#7570B3")

plot_fixef_oc <- ggplot(data = fixef_oc_tmeanpos_tmeanchange_realm) + 
  geom_linerangeh(aes(xmin = Q2.5, xmax = Q97.5, y = fct_rev(predictor), group = fct_rev(realm), color = realm),
                  position = position_dodge(width = 0.8, preserve = "total"), size =0.3, alpha = 1) + 
  geom_linerangeh(aes(xmin = Q5, xmax = Q95, y = fct_rev(predictor), group = fct_rev(realm), color = realm),
                  position = position_dodge(width = 0.8, preserve = "total"), size =0.6, alpha = 1) + 
  geom_linerangeh(aes(xmin = Q10, xmax = Q90, y = fct_rev(predictor), group = fct_rev(realm), color = realm),
                  position = position_dodge(width = 0.8, preserve = "total"), size =1, alpha = 1) + 
  geom_point(aes(x = Estimate, y = fct_rev(predictor), group = fct_rev(realm), colour = realm, shape = realm), 
             position = position_dodge(width = 0.8, preserve = "total"), size = 2, show.legend = T) + 
  geom_vline(xintercept = 0, lty = 2) + 
  labs(x = "Effect size on occupancy change", y = NULL, color = NULL, shape = NULL) + 
  scale_color_manual(values = realm_col) + 
  scale_shape_manual(values = c("Terrestrial" = 15, "Freshwater" = 16, "Marine" = 17)) + 
  #scale_shape_manual(values = c("Terrestrial" = 22, "Freshwater" = 21, "Marine" = 24)) + 
  theme_classic() +  
  theme(legend.position = c(0.25, 0.7),
        legend.background = element_blank(),
        legend.text = element_text(size = 7),
        axis.title.x = element_text(size = 9),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 9))

plot_fixef_ac <- ggplot(data = fixef_ac_tmeanpos_tmeanchange_realm) + 
  geom_linerangeh(aes(xmin = Q2.5, xmax = Q97.5, y = fct_rev(predictor), group = fct_rev(realm), color = realm),
                  position = position_dodge(width = 0.8, preserve = "total"), size =0.3, alpha = 1) + 
  geom_linerangeh(aes(xmin = Q5, xmax = Q95, y = fct_rev(predictor), group = fct_rev(realm), color = realm),
                  position = position_dodge(width = 0.8, preserve = "total"), size =0.6, alpha = 1) + 
  geom_linerangeh(aes(xmin = Q10, xmax = Q90, y = fct_rev(predictor), group = fct_rev(realm), color = realm),
                  position = position_dodge(width = 0.8, preserve = "total"), size =1, alpha = 1) + 
  geom_point(aes(x = Estimate, y = fct_rev(predictor), group = fct_rev(realm), color = realm, shape = realm), 
             position = position_dodge(width = 0.8, preserve = "total"), size = 2, show.legend = T) + 
  geom_vline(xintercept = 0, lty = 2) + 
  labs(x = "Effect size on abundance change", y = NULL, color = NULL, shape = NULL) + 
  scale_color_manual(values = realm_col) + 
  scale_shape_manual(values = c("Terrestrial" = 15, "Freshwater" = 16, "Marine" = 17)) + 
  theme_classic() +  
  theme(legend.position = "n",
        legend.background = element_blank(),
        legend.text = element_text(size = 7),
        axis.title.x = element_text(size = 9),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 9))

plot_fixef_oc_ac <- cowplot::plot_grid(plot_fixef_oc, plot_fixef_ac, nrow = 1, ncol= 2,
                                       labels = c("a", "b"), label_size = 10)

ggsave(plot_fixef_oc_ac, file = "results/supplement_Fig.S_fixef_thermal_position_realm.png", 
       width = 160, height = 70, units = 'mm', dpi = 300)



######################
## histogram of regional annual temperature changes across realms

rm(list = ls())
gc()

load("results/data_input_to_models.RDATA") # observation data

dat_meta <- dat_meta %>%
  mutate(realm = factor(realm, levels = c("Terrestrial", "Freshwater", "Marine")))
         
# colors for realm
realm_col <- c("Terrestrial" = "#D95F02", "Freshwater" = "#44bbff", "Marine" = "#7570B3")
realm_col <- c("All" = "#000000", "Terrestrial" = "#A65628", "Freshwater" = "#33A02C", "Marine" = "#1F78B4")

ggplot(dat_meta %>% mutate(change_tempmean = change_tempmean/(duration_mean -1))) +
  geom_histogram(aes(change_tempmean, fill = realm), alpha = 1) + 
  geom_vline(xintercept = 0, lty = 2) + 
  geom_vline(data = dat_meta  %>% 
               mutate(change_tempmean = change_tempmean/(duration_mean -1)) %>% 
               group_by(realm) %>% 
               summarise(change_tempmean = mean(change_tempmean)),
             aes(xintercept = change_tempmean, color = realm), lty = 2, size = 1) + 
  scale_fill_manual(name = NULL, values = realm_col) + 
  scale_colour_manual(name = NULL, values = realm_col) + 
  guides(color = FALSE) + 
  labs(x= expression('Temperature change ('*degree*'C/year)'), y = "Number of regions" )+ 
  #labs(x= bquote('Temperature change ('*degree*'C)'), ylab = "Number of regions" )+ 
  theme_classic() +  
  theme(legend.position = c(0.25, 0.7),
        legend.background = element_blank(),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8),
        legend.text = element_text(size = 10))

ggsave(file = "results/supplement_Fig.S_Temperature_change_periods_realms.png", 
       width = 90, height = 90, units = 'mm', dpi = 300)




################################
## visualize how occupancy/abundance changes vary with extent, duration, and n-sammples

rm(list = ls())
gc()

load("results/data_input_to_models.RDATA") # observation data

dat_meta <- dat_meta %>%
  mutate(realm = factor(realm, levels = c("Terrestrial", "Freshwater", "Marine")))

# colors for realm
realm_col <- c("Terrestrial" = "#D95F02", "Freshwater" = "#44bbff", "Marine" = "#7570B3")
realm_col <- c("All" = "#000000", "Terrestrial" = "#A65628", "Freshwater" = "#33A02C", "Marine" = "#1F78B4")

# plot
p_oc_extent <- ggplot(oc_period) + geom_bin2d(aes(x = extent_km2, y = occup_change)) + scale_x_log10() +
  labs( x = bquote('Spatial extent ('*km^2*')'), y = "Occupancy change", fill = "# populations") + 
  theme_classic() + theme(axis.title = element_text(size = 10))

p_oc_duration <- ggplot(oc_period) + geom_bin2d(aes(x = duration, y = occup_change)) + scale_x_log10() +
  labs( x = "Duration of sampling (years)", y = "Occupancy change", fill = "# populations") + 
  theme_classic() + theme(axis.title = element_text(size = 10))

p_oc_nsamp <- ggplot(oc_period) + geom_bin2d(aes(x = n_samp, y = occup_change)) + scale_x_log10() +
  labs( x = "Number of samples", y = "Occupancy change", fill = "# populations") + 
  theme_classic() + theme(axis.title = element_text(size = 10))

p_ac_extent <- ggplot(ac_period) + geom_bin2d(aes(x = extent_km2, y = N_change_log)) + scale_x_log10() +
  labs( x = bquote('Spatial extent ('*km^2*')'), y = "Abundance change", fill = "# populations") + 
  theme_classic() + theme(axis.title = element_text(size = 10))

p_ac_duration <- ggplot(ac_period) + geom_bin2d(aes(x = duration, y = N_change_log)) + scale_x_log10() +
  labs( x = "Duration of sampling (years)", y = "Abundance change", fill = "# populations") + 
  theme_classic() + theme(axis.title = element_text(size = 10))

p_ac_nsamp <- ggplot(ac_period) + geom_bin2d(aes(x = n_samp, y = N_change_log)) + scale_x_log10() +
  labs( x = "Number of samples", y = "Abundance change", fill = "# populations") + 
  theme_classic() + theme(axis.title = element_text(size = 10))

left <- cowplot::plot_grid(p_oc_extent, p_oc_duration, p_oc_nsamp, nrow = 3, labels = c("a", "c", "e"), label_size = 10)
right <- cowplot::plot_grid(p_ac_extent, p_ac_duration, p_ac_nsamp, nrow = 3, labels = c("b", "d", "f"), label_size = 10)
plot_oc_ac <- cowplot::plot_grid(left, right, ncol =2)

ggsave(plot_oc_ac, file = "results/supplement_Fig.S_OccupancyAbundanceChange_vs_extent_duration_nsamp.png", 
       width = 210, height = 210, units = 'mm', dpi = 300)




###################################
## Figure showing frequency distribution of study characteristics
# cent_lat, sprich, extent_km2, nsample_used, duration, start_year

rm(list = ls())
gc()

load("results/data_input_to_models.RDATA") # observation data

dat_meta <- dat_meta %>%
  mutate(realm = factor(realm, levels = c("Terrestrial", "Freshwater", "Marine")))

# colors for realm
realm_col <- c("Terrestrial" = "#D95F02", "Freshwater" = "#44bbff", "Marine" = "#7570B3")
realm_col <- c("All" = "#000000", "Terrestrial" = "#A65628", "Freshwater" = "#33A02C", "Marine" = "#1F78B4")

plot_meta_latitude <- ggplot(data = dat_meta) +
  geom_histogram(aes(cent_lat, fill = realm)) + 
  #scale_y_continuous(expand = expansion(mult = c(0.03, 0.05), add = c(0, 0))) + 
  labs(x = expression("Latitude ("*degree*")"), y = NULL, fill = NULL, tag = "a") + 
  scale_fill_manual(name = NULL, values = realm_col) + 
  theme_classic() +
  theme(legend.position = "no", 
        axis.title = element_text(size = 9, face = "plain"),
        axis.text = element_text(size = 8, face = "plain"),
        plot.tag = element_text(size = 9, face = "bold"))

plot_meta_richness <- ggplot(data = dat_meta) +
  geom_histogram(aes(sprich_thermal, fill = realm)) + 
  scale_x_log10() +
  labs(x = "Regional species richness", y = NULL, tag = "b") + 
  scale_fill_manual(name = NULL, values = realm_col) + 
  theme_classic() +
  theme(legend.position = "no", 
        axis.title = element_text(size = 9, face = "plain"),
        axis.text = element_text(size = 8, face = "plain"),
        plot.tag = element_text(size = 9, face = "bold"))

plot_meta_extent <- ggplot(data = dat_meta) +
  geom_histogram(aes(extent_km2, fill = realm)) + 
  scale_x_log10(breaks = c(10^0, 10^3, 10^6)) + 
  labs(x = bquote('Spatial extent ('*km^2*')'), y = NULL, tag = "c") + 
  scale_fill_manual(name = NULL, values = realm_col) + 
  theme_classic() +
  theme(legend.position = "no", 
        axis.title = element_text(size = 9, face = "plain"),
        axis.text = element_text(size = 8, face = "plain"),
        plot.tag = element_text(size = 9, face = "bold"))

plot_meta_nsite <- ggplot(data = dat_meta) +
  geom_histogram(aes(n_site, fill = realm)) + 
  scale_x_log10() +
  labs(x = "Number of sites", y = NULL, tag = "d") + 
  scale_fill_manual(name = NULL, values = realm_col) + 
  theme_classic() +
  theme(legend.position = "no", 
        axis.title = element_text(size = 9, face = "plain"),
        axis.text = element_text(size = 8, face = "plain"),
        plot.tag = element_text(size = 9, face = "bold"))

plot_meta_duration <- ggplot(data = dat_meta) +
  geom_histogram(aes(duration, fill = realm)) + 
  scale_x_log10() +
  labs(x = "Duration of sampling (years)", y = NULL, tag = "e") + 
  scale_fill_manual(name = NULL, values = realm_col) + 
  theme_classic() +
  theme(legend.position = "no", 
        axis.title = element_text(size = 9, face = "plain"),
        axis.text = element_text(size = 8, face = "plain"),
        plot.tag = element_text(size = 9, face = "bold"))

plot_meta_start <- ggplot(data = dat_meta) +
  geom_histogram(aes(start_year, fill = realm)) + 
  labs(x = "First year of sampling", y = NULL, tag = "f") + 
  scale_fill_manual(name = NULL, values = realm_col) + 
  theme_classic() +
  theme(legend.position = c(0.28, 0.85),
        legend.text = element_text(size = 8),
        legend.key.size = unit(4, "mm"),
        axis.title = element_text(size = 9, face = "plain"),
        axis.text = element_text(size = 8, face = "plain"),
        plot.tag = element_text(size = 9, face = "bold"))

# save the figure
left <- cowplot::plot_grid(NULL)
right <- cowplot::plot_grid(plot_meta_latitude, plot_meta_richness, plot_meta_extent,
                            plot_meta_nsite, plot_meta_duration, plot_meta_start,
                            nrow = 2, align = "hv") 

plot_study_characteristics <- cowplot::plot_grid(left, right, nrow = 1, rel_widths = c(0.025, 1)) + 
  cowplot::draw_label("Number of regions", x = 0.01, angle = 90, size = 10) +
  theme(plot.background = element_rect(fill = "white", colour = NA))

ggsave(plot_study_characteristics, file = "results/supplement_Fig.S_histogram_study_characteristics.png", width = 180, height = 125, units = 'mm')




######################
## histogram of regional annual temperature changes and study characteristics
# cent_lat, sprich, extent_km2, nsample_used, duration, start_year

rm(list = ls())
gc()

load("results/data_input_to_models.RDATA") # observation data

dat_meta <- dat_meta %>%
  mutate(realm = factor(realm, levels = c("Terrestrial", "Freshwater", "Marine")))

# colors for realm
realm_col <- c("Terrestrial" = "#D95F02", "Freshwater" = "#44bbff", "Marine" = "#7570B3")
realm_col <- c("All" = "#000000", "Terrestrial" = "#A65628", "Freshwater" = "#33A02C", "Marine" = "#1F78B4")

plot_meta_latitude <- ggplot(data = dat_meta) +
  geom_histogram(aes(cent_lat, fill = realm), alpha = 0.8) + 
  #scale_y_continuous(expand = expansion(mult = c(0.03, 0.05), add = c(0, 0))) + 
  labs(x = expression("Latitude ("*degree*")"), y = NULL, fill = NULL, tag = "a") + 
  scale_fill_manual(name = NULL, values = realm_col) + 
  theme_classic() +
  theme(legend.position = "no", 
        axis.title = element_text(size = 9, face = "plain"),
        axis.text = element_text(size = 8, face = "plain"),
        plot.tag = element_text(size = 10, face = "bold"))

plot_meta_richness <- ggplot(data = dat_meta) +
  geom_histogram(aes(sprich_thermal, fill = realm), alpha = 0.8) + 
  scale_x_log10() +
  labs(x = "Regional species richness", y = NULL, tag = "b") + 
  scale_fill_manual(name = NULL, values = realm_col) + 
  theme_classic() +
  theme(legend.position = "no", 
        axis.title = element_text(size = 9, face = "plain"),
        axis.text = element_text(size = 8, face = "plain"),
        plot.tag = element_text(size = 10, face = "bold"))

plot_meta_extent <- ggplot(data = dat_meta) +
  geom_histogram(aes(extent_km2, fill = realm), alpha = 0.8) + 
  scale_x_log10(breaks = c(10^0, 10^3, 10^6)) + 
  labs(x = bquote('Spatial extent ('*km^2*')'), y = NULL, tag = "c") + 
  scale_fill_manual(name = NULL, values = realm_col) + 
  theme_classic() +
  theme(legend.position = "no", 
        axis.title = element_text(size = 9, face = "plain"),
        axis.text = element_text(size = 8, face = "plain"),
        plot.tag = element_text(size = 10, face = "bold"))

plot_meta_nsite <- ggplot(data = dat_meta) +
  geom_histogram(aes(n_site, fill = realm), alpha = 0.8) + 
  scale_x_log10() +
  labs(x = "Number of sites", y = NULL, tag = "d") + 
  scale_fill_manual(name = NULL, values = realm_col) + 
  theme_classic() +
  theme(legend.position = "no", 
        axis.title = element_text(size = 9, face = "plain"),
        axis.text = element_text(size = 8, face = "plain"),
        plot.tag = element_text(size = 10, face = "bold"))

plot_meta_duration <- ggplot(data = dat_meta) +
  geom_histogram(aes(duration, fill = realm), alpha = 0.8) + 
  scale_x_log10() +
  labs(x = "Duration of sampling (years)", y = NULL, tag = "e") + 
  scale_fill_manual(name = NULL, values = realm_col) + 
  theme_classic() +
  theme(legend.position = "no", 
        axis.title = element_text(size = 9, face = "plain"),
        axis.text = element_text(size = 8, face = "plain"),
        plot.tag = element_text(size = 10, face = "bold"))

plot_meta_start <- ggplot(data = dat_meta) +
  geom_histogram(aes(start_year, fill = realm), alpha = 0.8) + 
  labs(x = "First year of sampling", y = NULL, tag = "f") + 
  scale_fill_manual(name = NULL, values = realm_col) + 
  theme_classic() +
  theme(legend.position = "no",
        legend.text = element_text(size = 8),
        legend.key.size = unit(4, "mm"),
        axis.title = element_text(size = 9, face = "plain"),
        axis.text = element_text(size = 8, face = "plain"),
        plot.tag = element_text(size = 10, face = "bold"))

# regional temperature
plot_meta_tempchange <- ggplot(dat_meta %>% mutate(change_tempmean = change_tempmean/(duration_mean -1))) +
  geom_histogram(aes(change_tempmean, fill = realm), alpha = 0.8) + 
  geom_vline(xintercept = 0, lty = 2) + 
  geom_vline(data = dat_meta  %>% 
               mutate(change_tempmean = change_tempmean/(duration_mean -1)) %>% 
               group_by(realm) %>% 
               summarise(change_tempmean = mean(change_tempmean)),
             aes(xintercept = change_tempmean, color = realm), lty = 2, size = 1) + 
  scale_fill_manual(name = NULL, values = realm_col) + 
  scale_colour_manual(name = NULL, values = realm_col) + 
  guides(color = FALSE) + 
  labs(x= expression('Temperature change ('*degree*'C/year)'),  y = NULL, tag = "g")+ 
  #labs(x= bquote('Temperature change ('*degree*'C)'), ylab = "Number of regions" )+ 
  theme_classic() +  
  theme(legend.position = "no",
        axis.title = element_text(size = 9),
        axis.text = element_text(size = 8),
        plot.tag = element_text(size = 10, face = "bold"))

# extract the legend realm
plot_legend_realm <- ggplot(data = dat_meta) +
  geom_histogram(aes(start_year, fill = realm), alpha = 0.8) + 
  labs(x = "First year of sampling", y = NULL) + 
  scale_fill_manual(name = NULL, values = realm_col) + 
  theme_classic() +
  theme(legend.text = element_text(size = 10),
        legend.key.size = unit(4, "mm"),
        axis.title = element_text(size = 9, face = "plain"),
        axis.text = element_text(size = 8, face = "plain"),
        plot.tag = element_text(size = 10, face = "bold"))

legend_realm <- get_legend(plot_legend_realm)


# save the figure
left <- cowplot::plot_grid(NULL)
right <- cowplot::plot_grid(plot_meta_latitude, plot_meta_richness, plot_meta_extent,
                            plot_meta_nsite, plot_meta_duration, plot_meta_start, 
                            plot_meta_tempchange,legend_realm,
                            nrow = 2) 

plot_study_characteristics <- cowplot::plot_grid(left, right, nrow = 1, rel_widths = c(0.02, 1)) + 
  cowplot::draw_label("Number of regions", x = 0.01, angle = 90, size = 10) +
  theme(plot.background = element_rect(fill = "white", colour = NA))

ggsave(plot_study_characteristics, 
       file = "results/supplement_Fig.S_histogram_study_characteristics_tempchange.png", 
       width = 200, height = 125, units = 'mm')



###########
##  histogram of study characteristics + density plots of thermal position and occupancy and abundance changes

density_tmeanpos <- ggplot(oc_period, aes(tmeanpos)) +
  geom_density(aes(y = after_stat(density),fill = realm), alpha = 0.8, position = "identity") + 
  geom_vline(data = oc_period %>% group_by(realm) %>% summarise(tmeanpos = mean(tmeanpos)),
             aes(xintercept = tmeanpos, color = realm), lty = 2, lwd = 1) + 
  scale_fill_manual(name = NULL, values = realm_col) + 
  scale_colour_manual(name = NULL, values = realm_col) + 
  guides(color = FALSE) + 
  labs(x = "Baseline thermal position", tag = "h") + 
  theme_classic() +  
  theme(legend.position = "no",
        legend.text = element_text(size = 8),
        axis.title = element_text(size = 9, face = "plain"),
        axis.text = element_text(size = 8, face = "plain"),
        plot.tag = element_text(size = 10, face = "bold"))

density_ocupancy_change <- ggplot(oc_period, aes(occup_change_logit)) +
  geom_density(aes(y = after_stat(density),fill = realm), alpha = 0.8, position = "identity") + 
  geom_vline(data = oc_period %>% group_by(realm) %>% summarise(occup_change_logit = mean(occup_change_logit)),
             aes(xintercept = occup_change_logit, color = realm), lty = 2, lwd = 1) + 
  xlim(-0.5, 0.5) +
  scale_fill_manual(name = NULL, values = realm_col) + 
  scale_colour_manual(name = NULL, values = realm_col) + 
  guides(color = FALSE) + 
  labs(x = "Occupancy change", tag = "i") + 
  theme_classic() +  
  theme(legend.position = "no",
        legend.text = element_text(size = 8),
        axis.title = element_text(size = 9, face = "plain"),
        axis.text = element_text(size = 8, face = "plain"),
        plot.tag = element_text(size = 10, face = "bold"))

density_abundance_change <- ggplot(ac_period, aes(N_change_log)) +
  geom_density(aes(y = after_stat(density),fill = realm), alpha = 0.8, position = "identity") + 
  geom_vline(data = ac_period %>% group_by(realm) %>% summarise(N_change_log = mean(N_change_log)),
             aes(xintercept = N_change_log, color = realm), lty = 2, lwd = 1) + 
  xlim(-0.5, 0.5) +
  scale_fill_manual(name = NULL, values = realm_col) + 
  scale_colour_manual(name = NULL, values = realm_col) + 
  guides(color = FALSE) + 
  labs(x = "Abundance change", tag = "j") + 
  theme_classic() +  
  theme(legend.position = "no",
        legend.text = element_text(size = 8),
        axis.title = element_text(size = 9, face = "plain"),
        axis.text = element_text(size = 8, face = "plain"),
        plot.tag = element_text(size = 10, face = "bold"))


bottom <- cowplot::plot_grid(density_tmeanpos, density_ocupancy_change, density_abundance_change, nrow = 1) 

plot_study_population_characteristics <- cowplot::plot_grid(plot_study_characteristics, bottom, nrow = 2, rel_heights = c(2, 1))

ggsave(plot_study_population_characteristics,
       file = "results/supplement_Fig.S_histogram_study_population_characteristics.png", 
       width = 200, height = 188, units = 'mm')



###################################
## Figure showing the distribution of species' thermal range and compare them between realms

rm(list = ls())
gc()

load("data/Species_thermal_limits.RDATA")
load("data/GBIF/Species_distribution_summary.RDATA")
load("data/Combined_assemblages.RDATA")

# colors for realm
realm_col <- c("Terrestrial" = "#D95F02", "Freshwater" = "#44bbff", "Marine" = "#7570B3")
realm_col <- c("All" = "#000000", "Terrestrial" = "#A65628", "Freshwater" = "#33A02C", "Marine" = "#1F78B4")

splimit <- dat %>% 
  distinct(species, specieskey, realm) %>% 
  left_join(splimit) %>% 
  mutate(realm = factor(realm, levels = c("terrestrial", "freshwater", "marine"),
                        labels = c("Terrestrial", "Freshwater", "Marine")))

splimit <- splimit %>% 
  left_join(spsuma %>% select(specieskey, n_clean, n_noyear, aoo10))


ggplot(splimit, aes(tempmean_range)) +
  geom_density(aes(y = after_stat(density),fill = realm), alpha = 0.8, position = "identity") + 
  # mean thermal range: terrestrial - 10.5, freshwater - 10.6, marine - 9.7
  geom_vline(data = splimit %>% group_by(realm) %>% summarise(tempmean_range = mean(tempmean_range)),
             aes(xintercept = tempmean_range, color = realm), lty = 2, lwd = 1) + 
  scale_fill_manual(name = NULL, values = realm_col) + 
  scale_colour_manual(name = NULL, values = realm_col) + 
  guides(color = FALSE) + 
  labs(x = expression('Oberved species thermal range ('*degree*'C)')) + 
  theme_classic() +  
  theme(legend.position = c(0.7, 0.7),
        legend.text = element_text(size = 10),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 9))

ggsave(file = "results/supplement_Fig.S_Thermal_range_realms.png", 
       width = 120, height = 90, units = 'mm', dpi = 300)


ggplot(splimit, aes(tempmean_range)) +
  geom_density(aes(y = after_stat(count),fill = realm), alpha = 0.8, position = "identity") + 
  # mean thermal range: terrestrial - 10.5, freshwater - 10.6, marine - 9.7
  geom_vline(data = splimit %>% group_by(realm) %>% summarise(tempmean_range = mean(tempmean_range)),
             aes(xintercept = tempmean_range, color = realm), lty = 2, lwd = 1) + 
  scale_fill_manual(name = NULL, values = realm_col) + 
  scale_colour_manual(name = NULL, values = realm_col) + 
  guides(color = FALSE) + 
  labs(x = expression('Oberved species thermal range ('*degree*'C)'),
       y = "Number of species") + 
  theme_classic() +  
  theme(legend.position = c(0.7, 0.7),
        legend.text = element_text(size = 10),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 9))

ggsave(file = "results/supplement_Fig.S_Thermal_range_realms_nspecies.png", 
       width = 120, height = 90, units = 'mm', dpi = 300)


ggplot(splimit) +
  geom_point(aes(x = nocc_climate, y = tempmean_range), alpha = 0.5) +
  scale_x_log10()

ggplot(splimit) + 
  geom_bin2d(aes(x = aoo10, y = tempmean_high, fill = log10(after_stat(count)))) +
  scale_x_log10() + 
  scale_fill_continuous(type = "viridis") + 
  labs(fill = bquote(Log[10]*'(number of populations)')) + 
  theme_bw() +
  theme(legend.position= "top",
        legend.key.size = unit(3, "mm"),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 8),
        legend.margin = margin(0,0,0,0),
        legend.box.spacing = unit(1, units = 'mm'),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())



###################################
## Figure showing the number of occurrences at the cold and warm distributional edges (10% of gid-cells with lowest and highest temperatures) 

rm(list = ls())
gc()

load("data/species_occurrences_at_thermal_edges.RDATA")
load("data/Combined_assemblages.RDATA")

nocc_thermal_egdes <- dat %>% 
  distinct(species, specieskey, realm) %>% 
  left_join(nocc_thermal_egdes) %>% 
  mutate(realm = factor(realm, levels = c("terrestrial", "freshwater", "marine"),
                        labels = c("Terrestrial", "Freshwater", "Marine")))

x <- dat %>% 
  distinct(species, specieskey, realm) %>% 
  left_join(nocc_thermal_egdes) %>% 
  mutate(nocc_avg_cold = 2*nocc_cold/nocc_10km,
         nocc_avg_warm = 2*nocc_warm/nocc_10km,)

x %>% 
  filter(nocc_10km > 20) %>% 
  group_by(realm) %>% 
  summarise(nocc_mean_cold = mean(nocc_avg_cold),
         nocc_mean_warm = mean(nocc_avg_warm),
         nocc_median_cold = median(nocc_avg_cold),
         nocc_median_warm = median(nocc_avg_warm))


ggplot(nocc_thermal_egdes) + 
  geom_bin2d(aes(x = nocc_cold, y = nocc_warm, fill = after_stat(count))) +
  geom_abline(aes(intercept = 0, slope = 1), col = "red", lty =2, lwd = 1) +
  scale_x_log10() + 
  scale_y_log10() + 
  scale_fill_continuous(type = "viridis") + 
  labs(x = "# occurrences at cold thermal range",
       y = "# occurrences at warm thermal range",
       fill = 'Number of species') + 
  theme_bw() +
  theme(legend.position= c(0.3, 0.75),
        legend.background = element_blank(),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave(file = "results/supplement_Fig.S_occurrences_at_thermal_edges.png", 
       width = 100, height = 100, units = 'mm', dpi = 300)


ggplot(nocc_thermal_egdes) + 
  facet_wrap(~ realm) +
  geom_bin2d(aes(x = nocc_cold, y = nocc_warm, fill = after_stat(count))) +
  geom_abline(aes(intercept = 0, slope = 1), col = "red", lty =2, lwd = 1) +
  scale_x_log10() + 
  scale_y_log10() + 
  scale_fill_continuous(type = "viridis") + 
  labs(x = "# occurrences at cold thermal range",
       y = "# occurrences at warm thermal range",
       fill = 'Number of species') + 
  theme_bw() +
  theme(legend.position = c(0.8, 1.2), 
        plot.margin = unit(c(1,0.2,0.2,0.2), "cm"),
        legend.direction = 'horizontal',
        legend.key.size = unit(4, 'mm'),
        legend.background = element_blank(),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 9),
        text = element_text(size = 10),
        strip.text = element_text(size = 9),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  guides(color = guide_colourbar(barwidth = 4, barheight = 0.5))

ggsave(file = "results/supplement_Fig.S_occurrences_at_thermal_edges_realms.png", 
       width = 180, height = 80, units = 'mm', dpi = 300)



####################################
## save summary of models as tables   
rm(list = ls())
gc()

load("results/data_input_to_models.RDATA") # observation data
brm_oc_tmeanpos <- readRDS("models/brm_oc_tmeanpos.rds")
brm_ac_tmeanpos <- readRDS("models/brm_ac_tmeanpos.rds")
brm_oc_tmeanpos_realm <- readRDS("models/brm_oc_tmeanpos_realm.rds")
brm_ac_tmeanpos_realm <- readRDS("models/brm_ac_tmeanpos_realm.rds")
brm_oc_tmeanpos_tmeanchange_realm <- readRDS("models/brm_oc_tmeanpos_tmeanchange_realm.rds")
brm_ac_tmeanpos_tmeanchange_realm <- readRDS("models/brm_ac_tmeanpos_tmeanchange_realm.rds")

brm_slope_oc_tmeanchange_realm <- readRDS("models/brm_slope_oc_tmeanchange_realm.rds")
brm_slope_ac_tmeanchange_realm <- readRDS("models/brm_slope_ac_tmeanchange_realm.rds")
brm_slope_oc_covariates_realm <- readRDS("models/brm_slope_oc_covariates_realm.rds")
brm_slope_ac_covariates_realm <- readRDS("models/brm_slope_ac_covariates_realm.rds")


summary(brm_oc_tmeanpos_realm)[["ngrps"]]
summary(brm_oc_tmeanpos_realm)[["nobs"]]
oc_period %>% distinct(realm, study) %>% group_by(realm) %>% summarise(n())
oc_period %>% group_by(realm) %>% summarise(n())

# model summaries
## overall model
brmsummary_oc_tmeanpos <- summary(brm_oc_tmeanpos)[["fixed"]] %>%
  as.data.frame() %>%
  rownames_to_column(var = "Parameter") %>%
  mutate(Parameter = gsub("realm", "", Parameter)) %>%
  mutate(across(Estimate:'Rhat', \(x) round(x, 3))) %>%
  mutate(across(Bulk_ESS:Tail_ESS, \(x) round(x, 0)))

brmsummary_ac_tmeanpos <- summary(brm_ac_tmeanpos)[["fixed"]] %>%
  as.data.frame() %>%
  rownames_to_column(var = "Parameter") %>%
  mutate(Parameter = gsub("realm", "", Parameter)) %>%
  mutate(across(Estimate:'Rhat', \(x) round(x, 3))) %>%
  mutate(across(Bulk_ESS:Tail_ESS, \(x) round(x, 0)))

# models include differences among realms
brmsummary_oc_tmeanpos_realm <- summary(brm_oc_tmeanpos_realm)[["fixed"]] %>%
  as.data.frame() %>%
  rownames_to_column(var = "Parameter") %>%
  mutate(Parameter = gsub("realm", "", Parameter)) %>%
  mutate(across(Estimate:'Rhat', \(x) round(x, 3))) %>%
  mutate(across(Bulk_ESS:Tail_ESS, \(x) round(x, 0)))

brmsummary_ac_tmeanpos_realm <- summary(brm_ac_tmeanpos_realm)[["fixed"]] %>%
  as.data.frame() %>%
  rownames_to_column(var = "Parameter") %>%
  mutate(Parameter = gsub("realm", "", Parameter)) %>%
  mutate(across(Estimate:'Rhat', \(x) round(x, 3))) %>%
  mutate(across(Bulk_ESS:Tail_ESS, \(x) round(x, 0)))

# models include an interaction between thermal position and temperature change
brmsummary_oc_tmeanpos_tmeanchange_realm <- summary(brm_oc_tmeanpos_tmeanchange_realm)[["fixed"]] %>%
  as.data.frame() %>%
  rownames_to_column(var = "Parameter") %>%
  mutate(Parameter = gsub("realm", "", Parameter)) %>%
  mutate(across(Estimate:'Rhat', \(x) round(x, 3))) %>%
  mutate(across(Bulk_ESS:Tail_ESS, \(x) round(x, 0)))

brmsummary_ac_tmeanpos_tmeanchange_realm <- summary(brm_ac_tmeanpos_tmeanchange_realm)[["fixed"]] %>%
  as.data.frame() %>%
  rownames_to_column(var = "Parameter") %>%
  mutate(Parameter = gsub("realm", "", Parameter)) %>%
  mutate(across(Estimate:'Rhat', \(x) round(x, 3))) %>%
  mutate(across(Bulk_ESS:Tail_ESS, \(x) round(x, 0)))


# models analyzing variation in study-level slopes
brmsummary_slope_oc_tmeanchange_realm <- summary(brm_slope_oc_tmeanchange_realm)[["fixed"]] %>%
  as.data.frame() %>%
  rownames_to_column(var = "Parameter") %>%
  mutate(Parameter = gsub("realm", "", Parameter)) %>%
  mutate(across(Estimate:'Rhat', \(x) round(x, 3))) %>%
  mutate(across(Bulk_ESS:Tail_ESS, \(x) round(x, 0)))

brmsummary_slope_ac_tmeanchange_realm <- summary(brm_slope_ac_tmeanchange_realm)[["fixed"]] %>%
  as.data.frame() %>%
  rownames_to_column(var = "Parameter") %>%
  mutate(Parameter = gsub("realm", "", Parameter)) %>%
  mutate(across(Estimate:'Rhat', \(x) round(x, 3))) %>%
  mutate(across(Bulk_ESS:Tail_ESS, \(x) round(x, 0)))

brmsummary_slope_oc_covariates_realm <- summary(brm_slope_oc_covariates_realm)[["fixed"]] %>%
  as.data.frame() %>%
  rownames_to_column(var = "Parameter") %>%
  mutate(Parameter = gsub("realm", "", Parameter)) %>%
  mutate(across(Estimate:'Rhat', \(x) round(x, 4))) %>%
  mutate(across(Bulk_ESS:Tail_ESS, \(x) round(x, 0)))

brmsummary_slope_ac_covariates_realm <- summary(brm_slope_ac_covariates_realm)[["fixed"]] %>%
  as.data.frame() %>%
  rownames_to_column(var = "Parameter") %>%
  mutate(Parameter = gsub("realm", "", Parameter)) %>%
  mutate(across(Estimate:'Rhat', \(x) round(x, 4))) %>%
  mutate(across(Bulk_ESS:Tail_ESS, \(x) round(x, 0)))

# save tables
dir.create("results/tables", showWarnings = FALSE)
write.csv(brmsummary_oc_tmeanpos, "Results/tables/model_summary_brm_oc_tmeanpos.csv")
write.csv(brmsummary_ac_tmeanpos, "Results/tables/model_summary_brm_ac_tmeanpos.csv")
write.csv(brmsummary_oc_tmeanpos_realm, "Results/tables/model_summary_brm_oc_tmeanpos_realm.csv")
write.csv(brmsummary_ac_tmeanpos_realm, "Results/tables/model_summary_brm_ac_tmeanpos_realm.csv")
write.csv(brmsummary_oc_tmeanpos_tmeanchange_realm, "Results/tables/model_summary_brm_oc_tmeanpos_tmeanchange_realm.csv")
write.csv(brmsummary_ac_tmeanpos_tmeanchange_realm, "Results/tables/model_summary_brm_ac_tmeanpos_tmeanchange_realm.csv")
write.csv(brmsummary_slope_oc_tmeanchange_realm, "Results/tables/model_summary_brm_slope_oc_tmeanchange_realm.csv")
write.csv(brmsummary_slope_ac_tmeanchange_realm, "Results/tables/model_summary_brm_slope_ac_tmeanchange_realm.csv")
write.csv(brmsummary_slope_oc_covariates_realm, "Results/tables/model_summary_brm_slope_oc_covariates_realm.csv")
write.csv(brmsummary_slope_ac_covariates_realm, "Results/tables/model_summary_brm_slope_ac_covariates_realm.csv")
