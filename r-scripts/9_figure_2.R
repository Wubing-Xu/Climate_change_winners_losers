
rm(list = ls())
gc()

# load packages
needed_libs <- c("tidyverse","ggplot2", "brms", "tidybayes", "abind", "rnaturalearth", "sf", "cowplot", "ggridges")

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

# input the observation data and models
load("results/data_input_to_models.RDATA")


#######################
## figure 2: maps show regional-level temperature change and slopes of occupancy/abundance change aganist baseline thermal position

brm_oc_tmeanpos <- readRDS("models/brm_oc_tmeanpos.rds")
brm_ac_tmeanpos <- readRDS("models/brm_ac_tmeanpos.rds")
brm_oc_tmeanpos_realm <- readRDS("models/brm_oc_tmeanpos_realm.rds")
brm_ac_tmeanpos_realm <- readRDS("models/brm_ac_tmeanpos_realm.rds")

## study-level temperature changes
dat_tempchange <- oc_period %>% 
  select(c(study:n_samp, duration_mean, taxon_new:change_tempmean, realm_taxa, realm_climate)) %>% distinct() %>% 
  left_join(dat_meta %>% dplyr::select(study, cent_lat, cent_long)) %>%
  relocate(study, database, studyID, study_name)


## predicted occupancy and abundance changes

# sequence of thermal position
position_seq <- seq(0, 1, length = 101)

# new data for predicting the overall relations
newdata_global <- tibble(tmeanpos = position_seq, 
                         oc_period %>% distinct(study, realm, n_samp, duration, extent_km2) %>% 
                           summarise(n_samp = round(exp(mean(log(n_samp)))), # use the geometric mean 
                                     duration = round(exp(mean(log(duration)))),
                                     extent_km2 = round(exp(mean(log(extent_km2))), 1)))

# new data for predicting relations across realms
newdata_realm <- oc_period %>% distinct(study, realm, n_samp, duration, extent_km2) %>% group_by(realm) %>% 
  summarise(n_samp = round(exp(mean(log(n_samp)))), # use the geometric mean 
            duration = round(exp(mean(log(duration)))),
            extent_km2 = round(exp(mean(log(extent_km2))), 1)) %>% 
  left_join(tibble(tmeanpos = rep(position_seq, 3), 
                   realm = rep(c("Terrestrial", "Freshwater", "Marine"), each = length(position_seq))))


## text showing fixed effect
## for occupancy change 
beta_oc_tmeanpos_global <- fixef(brm_oc_tmeanpos) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "term") %>%
  filter(grepl("tmeanpos", term)) %>%
  mutate(label = paste0("beta == ", round(Estimate, 3),'~(', round(Q2.5, 3), "*','~", round(Q97.5, 3),')'))

beta_oc_tmeanpos_realm <- fixef(brm_oc_tmeanpos_realm, robust = TRUE, probs = c(0.025, 0.975)) %>%
  as.data.frame() %>% 
  rownames_to_column(var = "term") %>%
  filter(grepl(":tmeanpos", term)) %>%
  mutate(realm = gsub(":tmeanpos", "", term),
         realm = gsub("realm", "", realm),
         realm = factor(realm, levels = c("Terrestrial", "Freshwater", "Marine"))) %>%
  mutate(label = paste0("beta == ", round(Estimate, 3),'~(', round(Q2.5, 3), "*','~", round(Q97.5, 3),')'))

beta_oc_tmeanpos <- data.frame(realm = "All", beta_oc_tmeanpos_global) %>% 
  bind_rows(beta_oc_tmeanpos_realm) %>% 
  mutate(realm = factor(realm, levels = c("All", "Terrestrial", "Freshwater", "Marine")))


## for abundance change
beta_ac_tmeanpos_global <- fixef(brm_ac_tmeanpos) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "term") %>%
  filter(grepl("tmeanpos", term)) %>%
  mutate(label = paste0('italic(beta) == ', round(Estimate, 3),'~(', round(Q2.5, 3), "*','~", round(Q97.5, 3),')'))

beta_ac_tmeanpos_realm <- fixef(brm_ac_tmeanpos_realm, robust = TRUE, probs = c(0.025, 0.975)) %>%
  as.data.frame() %>% 
  rownames_to_column(var = "term") %>%
  filter(grepl(":tmeanpos", term)) %>%
  mutate(realm = gsub(":tmeanpos", "", term),
         realm = gsub("realm", "", realm),
         realm = factor(realm, levels = c("Terrestrial", "Freshwater", "Marine"))) %>%
  mutate(label = paste0('italic(beta) == ', round(Estimate, 3),'~(', round(Q2.5, 3), "*','~", round(Q97.5, 3),')'))

beta_ac_tmeanpos <- data.frame(realm = "All", beta_ac_tmeanpos_global) %>% 
  bind_rows(beta_ac_tmeanpos_realm) %>% 
  mutate(realm = factor(realm, levels = c("All", "Terrestrial", "Freshwater", "Marine")))


## number of studies used for each realm
nstudy_oc <- oc_period %>% 
  distinct(study, realm, n_site) %>% 
  group_by(realm) %>%
  summarise(nstudy = n(),
            nsite = sum(n_site)) %>% 
  bind_rows(summarise(.,
              across(where(is.numeric), sum),
              realm = "All")) %>%
mutate(label = paste0("italic(n) == ", nstudy),
       realm = factor(realm, levels = c("All", "Terrestrial", "Freshwater", "Marine")))

nstudy_ac <- ac_period %>% 
  distinct(study, realm, n_site) %>% 
  group_by(realm) %>%
  summarise(nstudy = n(),
            nsite = sum(n_site)) %>% 
  bind_rows(summarise(.,
                      across(where(is.numeric), sum),
                      realm = "All")) %>%
  mutate(label = paste0("italic(n) == ", nstudy),
         realm = factor(realm, levels = c("All", "Terrestrial", "Freshwater", "Marine")))


## fitted values along thermal position
## for occupancy change
fitted_oc_tmeanpos_global <- fitted(brm_oc_tmeanpos, newdata = newdata_global, re_formula = NA, ndraws = 1000, probs = c(0.025, 0.1, 0.9, 0.975)) %>%
  as_tibble()  %>% 
  bind_cols(newdata_global)

fitted_oc_tmeanpos_realm <- fitted(brm_oc_tmeanpos_realm, newdata = newdata_realm, re_formula = NA, ndraws = 1000, probs = c(0.025, 0.1, 0.9, 0.975)) %>%
  as_tibble()  %>% 
  bind_cols(newdata_realm)

fitted_oc_tmeanpos <- tibble(realm = "All", fitted_oc_tmeanpos_global) %>% 
  bind_rows(fitted_oc_tmeanpos_realm) %>% 
  mutate(realm = factor(realm, levels = c("All", "Terrestrial", "Freshwater", "Marine")))


## for abundance change
fitted_ac_tmeanpos_global <- fitted(brm_ac_tmeanpos, newdata = newdata_global, re_formula = NA, ndraws = 1000, probs = c(0.025, 0.1, 0.9, 0.975)) %>%
  as_tibble()  %>% 
  bind_cols(newdata_global)

fitted_ac_tmeanpos_realm <- fitted(brm_ac_tmeanpos_realm, newdata = newdata_realm, re_formula = NA, ndraws = 1000, probs = c(0.025, 0.1, 0.9, 0.975)) %>%
  as_tibble()  %>% 
  bind_cols(newdata_realm)

fitted_ac_tmeanpos <- tibble(realm = "All", fitted_ac_tmeanpos_global) %>% 
  bind_rows(fitted_ac_tmeanpos_realm) %>% 
  mutate(realm = factor(realm, levels = c("All", "Terrestrial", "Freshwater", "Marine")))



## draw figure

# colors for realm
realm_col <- c("All" = "#000000", "Terrestrial" = "#D95F02", "Freshwater" = "#44bbff", "Marine" = "#7570B3")
realm_col <- c("All" = "#000000", "Terrestrial" = "#E78AC3", "Freshwater" = "#FFD92F", "Marine" = "#66C2A5")
realm_col <- c("All" = "#000000", "Terrestrial" = "#A65628", "Freshwater" = "#33A02C", "Marine" = "#1F78B4")

# world map 
worldmap <- ne_countries(scale = 'medium', type = 'map_units',  returnclass = 'sf')
worldmap <- worldmap %>% filter(worldmap$admin != "Antarctica")

# generate a spatial point file containing temperature changes
dat_tempchange_points <- dat_tempchange %>% 
  dplyr::select(cent_long, cent_lat, realm, change_tempmean) %>% 
  mutate(change_tempmean = ifelse(change_tempmean <= -0.1, -0.1, ifelse(change_tempmean >= 0.1, 0.1, change_tempmean))) %>%
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
  geom_point(data = dat_tempchange, aes(x = cent_long, y = cent_lat, shape = realm), fill = "gray50") + 
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


###### map of temperature change
map_tempchange_global <- ggplot() + 
  geom_sf(data = worldmap, fill = "gray95", colour = "gray50", size = 0.2) + 
  geom_sf(data = namerica, fill = NA, colour = "black") + 
  geom_sf(data = europe, fill = NA, colour = "black") + 
  geom_sf(data = dat_tempchange_points, aes(fill = change_tempmean, shape = realm), color = "black",
          size = 1.5, alpha = 0.8) + 
  coord_sf(crs = st_crs('+proj=moll'), ylim = c(-8172663, 9020048), expand = FALSE) + 
  scale_fill_gradient2(low = "blue", mid = "darkgray", high = "red") + #breaks = c(-2, -1, 0, 1, 2), 
  scale_shape_manual(values = c("Terrestrial" = 22, "Freshwater" = 21, "Marine" = 24)) + 
#  annotate("text", x =1.125*-11046300, y= 0.80*6386580, label = "North\nAmerica", size = 2.46) +
#  annotate("text", x =1.75*-1469518, y= 0.90*6876759, label = "Europe", size = 2.46) +
  theme_bw() +
  theme(plot.margin = unit(c(0,-0.5,0,-1), "cm"),
        legend.position= c(0.55, 0.10), 
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
  guides(shape = "none") + 
  guides(fill = guide_colourbar(title = expression('Temperature change ('*degree*'C/year)'), title.position = "top", barwidth = 5))

histgram_tempchange <- ggplot(data = dat_tempchange) +
  geom_histogram(aes(change_tempmean), fill = "gray60") +
  geom_vline(xintercept = 0, lty = 2) + 
  geom_vline(xintercept = mean(dat_tempchange$change_tempmean), linetype =1) + 
  labs(x= expression(Delta*'T ('*degree*'C/year)'), y = "Number of regions" )+ 
  theme_classic() +  
  theme(axis.title = element_text(size = 7),
        axis.text = element_text(size = 6),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank())

map_tempchange_global <- ggdraw(map_tempchange_global) +
  draw_plot(histgram_tempchange, x = 0.01, y = 0.01, width = 0.3, height = 0.5)

# map of north america
map_tempchange_namerica <- ggplot() + 
  geom_sf(data = worldmap, fill = "gray95", colour = "gray50", size = 0.2) + 
  geom_sf(data = dat_tempchange_points, aes(fill = change_tempmean, shape = realm), color = "black",
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
map_tempchange_europe <- ggplot() + 
  geom_sf(data = worldmap, fill = "gray95", colour = "gray50", size = 0.2) + 
  geom_sf(data = dat_tempchange_points, aes(fill = change_tempmean, shape = realm), color = "black",
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

right <- cowplot::plot_grid(map_tempchange_namerica, map_tempchange_europe, ncol =1, nrow = 2, rel_heights = c(1, 1.135))
map_tempchange <- cowplot::plot_grid(map_tempchange_global, right,  ncol =2, nrow = 1, rel_widths = c(2.52, 1),
                                 hjust = 0, vjust = 1) +
  theme(plot.background = element_rect(fill = "white", colour = NA))

# add the legend of realm
map_tempchange_legend <- cowplot::plot_grid(legend_realm, map_tempchange, ncol =1, nrow = 2, rel_heights = c(0.1, 1)) +
  theme(plot.background = element_rect(fill = "white", colour = NA))

ggsave(map_tempchange_legend, file = "results/Fig.2a.png", 
       width = 180, height = 72.8, units = 'mm', dpi = 600) # 2.4725:1


## occupancy changes
plot_oc <- ggplot(fitted_oc_tmeanpos) +
  facet_wrap(~ realm, nrow = 1) + 
  geom_hline(yintercept = 0, lty = 2, linewidth = 0.3) + 
  geom_ribbon(aes(x =  tmeanpos, ymin = Q2.5, ymax = Q97.5, fill = realm), alpha = 0.3) + 
  geom_line(aes(x = tmeanpos, y = Estimate, color = realm), lwd = 1) + 
  geom_text(data = beta_oc_tmeanpos, mapping = aes(x= 0.12, y = Inf, label = label), 
            hjust = 0, vjust = 1.5,  parse = T, size = 2.46) + 
  geom_text(data = nstudy_oc, mapping = aes(x= -Inf, y = -Inf, label = label), 
            hjust = -0.2, vjust = -0.5,  parse = T, size = 2.46) + 
  # scale_y_continuous(limits = c(-0.28,0.4)) +
  labs(x = "Baseline thermal position", y = "Occupancy change") + 
  scale_colour_manual(name = NULL, values = realm_col) + 
  scale_fill_manual(name = NULL, values = realm_col) + 
  theme_bw()  +
  theme(legend.position = "none", 
        text = element_text(size =9),
        strip.text = element_text(size = 9),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 9), 
        # panel.background = element_rect(fill = 'transparent', colour = NA),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())


## abundance changes
plot_ac <- ggplot(fitted_ac_tmeanpos) +
  facet_wrap(~ realm, nrow = 1) + 
  geom_hline(yintercept = 0, lty = 2, linewidth = 0.3) + 
  geom_ribbon(aes(x =  tmeanpos, ymin = Q2.5, ymax = Q97.5, fill = realm), alpha = 0.3) + 
  geom_line(aes(x = tmeanpos, y = Estimate,  color = realm), lwd = 1) + 
  geom_text(data = beta_ac_tmeanpos, mapping = aes(x= 0.12, y = Inf, label = label), 
            hjust = 0, vjust = 1.5,  parse = TRUE, size = 2.46) + 
  geom_text(data = nstudy_ac, mapping = aes(x= -Inf, y = -Inf, label = label), 
            hjust = -0.2, vjust = -0.5,  parse = T, size = 2.46) + 
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


# combine panels
plot_tempchange_oc_ac <- cowplot::plot_grid(map_tempchange_legend, plot_oc, plot_ac, 
                                            labels = c("a", "b", "c"), label_size = 10,
                                            nrow = 3, rel_heights = c(0.728, 0.5, 0.5))

ggsave(plot_tempchange_oc_ac, file = "results/Fig.2.png", 
       width = 180, height = 172, units = 'mm', dpi = 600) # 1:0.96

ggsave(plot_tempchange_oc_ac, file = "results/Fig.2.pdf", 
       width = 180, height = 170, units = 'mm') # 1:0.9233379




########################
## explortory analyses

quantile(oc_period$occup_change_logit, probs = c(0.1,0.9))
ggplot(oc_period %>% filter(occup_change_logit >-0.218 & occup_change_logit < 0.2498)) + 
  geom_bin2d(aes(x = tmeanpos, y = occup_change_logit, fill = log10(..count..))) 

ggplot(oc_period %>% filter(occup_change_logit >-0.1 & occup_change_logit < 0.1)) + 
  # facet_wrap(~realm, nrow = 3) + 
  geom_bin2d(aes(x = tmeanpos, y = occup_change_logit, fill = log10(after_stat(count)))) +
  scale_fill_continuous(type = "viridis") + 
  labs(x = "Baseline thermal position",
       y = "Occupancy change",
       fill = bquote(Log[10]*'(number of populations)')) + 
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


ggplot(ac_period %>% filter(N_change_log >-0.1 & N_change_log < 0.1)) + 
  geom_bin2d(aes(x = tmeanpos, y = N_change_log, fill = log10(after_stat(count)))) +
  scale_fill_continuous(type = "viridis") + 
  labs(x = "Baseline thermal position",
       y = "Abundance change",
       fill = bquote(Log[10]*'(number of populations)')) + 
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

