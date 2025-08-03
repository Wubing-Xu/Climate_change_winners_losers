
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


######################
# Fig. 4: predicted occupancy/abundance change across thermal positions and temperature changes

brm_oc_tmeanpos_tmeanchange_realm <- readRDS("models/brm_oc_tmeanpos_tmeanchange_realm.rds")
brm_ac_tmeanpos_tmeanchange_realm <- readRDS("models/brm_ac_tmeanpos_tmeanchange_realm.rds")

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
fitted_oc_tmeanpos <- fitted(brm_oc_tmeanpos_tmeanchange_realm, newdata = newdata, re_formula = NA, ndraws = 1000, probs = c(0.025, 0.975)) %>%
  as_tibble()  %>% 
  bind_cols(newdata)

fitted_ac_tmeanpos <- fitted(brm_ac_tmeanpos_tmeanchange_realm, newdata = newdata, re_formula = NA, ndraws = 1000, probs = c(0.025, 0.975)) %>%
  as_tibble()  %>% 
  bind_cols(newdata)

# text showing interaction effect
beta_oc_interaction_text <- fixef(brm_oc_tmeanpos_tmeanchange_realm) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "term") %>%
  filter(grepl("tmeanpos:change_tempmean", term)) %>%
  mutate(realm = gsub(":tmeanpos:change_tempmean", "", term),
         realm = gsub("realm", "", realm),
         realm = factor(realm, levels = c("Terrestrial", "Freshwater", "Marine"))) %>%
  mutate(label = paste('beta == ', round(Estimate, 2),'~(', round(Q2.5, 2), "*','~", round(Q97.5, 2),')'))

beta_ac_interaction_text <- fixef(brm_ac_tmeanpos_tmeanchange_realm) %>% 
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
  scale_y_continuous(limits = c(-0.085, 0.09), breaks = c( -0.08, 0, 0.08)) +
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
  scale_y_continuous(limits = c(-0.055, 0.065), breaks = c( -0.05, 0, 0.05)) +
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

ggsave(file = "results/Fig.4.png", 
       width = 150, height = 120, units = 'mm', dpi = 300)

ggsave(file = "results/Fig.4.pdf", 
       width = 150, height = 120, units = 'mm')
