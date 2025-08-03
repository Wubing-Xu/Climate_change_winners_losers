
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


##############################
## figure 3: 

# the regional-level slopes of baseline thermal position was extracted in the r-script for fig. 2
load("results/coef_oc_ac_tmeanpos.RDATA")

## fit models firstly: include uncertainty of the study-level estimate in the model 
brm_slope_oc_tmeanchange_realm <- brm(bf(estimate_slope|se(se_slope) ~ 0 + realm + change_tempmean:realm + (1|study)),
                                      data = coef_oc_tmeanpos,
                                      control = list(adapt_delta = 0.9, max_treedepth = 10),
                                      cores = 4, chains = 4, iter = 8000, thin = 2,
                                      file = "models/brm_slope_oc_tmeanchange_realm")
brm_slope_oc_tmeanchange_realm
fixef(brm_slope_oc_tmeanchange_realm, prob = c(0.025, 0.975, 0.1, 0.9))

brm_slope_ac_tmeanchange_realm <- brm(bf(estimate_slope|se(se_slope) ~ 0 + realm + change_tempmean:realm + (1|study)),
                                      data = coef_ac_tmeanpos,
                                      control = list(adapt_delta = 0.9, max_treedepth = 10),
                                      cores = 4, chains = 4, iter = 8000, thin = 2,
                                      file = "models/brm_slope_ac_tmeanchange_realm")
brm_slope_ac_tmeanchange_realm
fixef(brm_slope_ac_tmeanchange_realm, prob = c(0.025, 0.975, 0.1, 0.9))



## get the fitted values along temperature change
fitted_slope_oc_tmeanchange_realm <- cbind(brm_slope_oc_tmeanchange_realm$data,
                                           fitted(brm_slope_oc_tmeanchange_realm, re_formula = NA)) %>% 
  as_tibble() %>% 
  mutate(realm = factor(realm, levels = c("Terrestrial", "Freshwater", "Marine")))

fitted_slope_ac_tmeanchange_realm <- cbind(brm_slope_ac_tmeanchange_realm$data,
                                           fitted(brm_slope_ac_tmeanchange_realm, re_formula = NA)) %>% 
  as_tibble() %>% 
  mutate(realm = factor(realm, levels = c("Terrestrial", "Freshwater", "Marine")))

## get the slopes of regressions
beta_slope_oc_tmeanchange_text <- fixef(brm_slope_oc_tmeanchange_realm) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "term") %>%
  filter(grepl("change_tempmean", term)) %>%
  mutate(realm = gsub(":change_tempmean", "", term),
         realm = gsub("realm", "", realm),
         realm = factor(realm, levels = c("Terrestrial", "Freshwater", "Marine"))) %>%
  mutate(label = paste('beta == ', round(Estimate, 3),'~(', round(Q2.5, 3), "*','~", round(Q97.5, 3),')'))

beta_slope_ac_tmeanchange_text <- fixef(brm_slope_ac_tmeanchange_realm) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "term") %>%
  filter(grepl("change_tempmean", term)) %>%
  mutate(realm = gsub(":change_tempmean", "", term),
         realm = gsub("realm", "", realm),
         realm = factor(realm, levels = c("Terrestrial", "Freshwater", "Marine"))) %>%
  mutate(label = paste('beta == ', round(Estimate, 3),'~(', round(Q2.5, 3), "*','~", round(Q97.5, 3),')'))


## draw figures
# colors for realm
realm_col <- c("Terrestrial" = "#D95F02", "Freshwater" = "#44bbff", "Marine" = "#7570B3")
realm_col <- c("All" = "#000000", "Terrestrial" = "#A65628", "Freshwater" = "#33A02C", "Marine" = "#1F78B4")


plot_coef_oc_tmeanpos_tmeanchange <- ggplot(fitted_slope_oc_tmeanchange_realm) + 
  facet_wrap(~ realm, scales = "free_x") + 
  geom_point(aes(x = change_tempmean, y = estimate_slope, fill = realm, size = 1/se_slope), shape = 21, color = "black", alpha = 0.7) +
  geom_ribbon(aes(x = change_tempmean, ymin = Q2.5, ymax = Q97.5, fill = realm), alpha = 0.3) + 
  geom_line(aes(x = change_tempmean, y = Estimate, color = realm), lwd = 1) + 
  geom_vline(xintercept = 0, linetype =2, colour = "gray", lwd = 0.4) + 
  geom_hline(yintercept = 0, linetype =2, colour = "gray", lwd = 0.4) + 
  geom_text(data = beta_slope_oc_tmeanchange_text, mapping = aes(x= 0.161, y = Inf, label = label), 
            hjust = 1, vjust = 1.3,  parse = T, size = 2.46) + 
  scale_y_continuous(limits = c(-0.2, 0.2), breaks = c( -0.15, 0, 0.15)) + # not show a value below -0.2 in the marine
  labs(y = "Occupancy change ~\n f(baseline thermal position)", x = expression('Temperature change ('*degree*'C/year)')) + 
  scale_colour_manual(name = NULL, values = realm_col) + 
  scale_fill_manual(name = NULL, values = realm_col) + 
  scale_size(range = c(0.3,2.5)) +
  theme_bw() + 
  theme(legend.position = "n", 
        #strip.text.x  = element_blank(),
        #line = element_line(size = 0.4),
        #strip.background = element_blank(),
        strip.text = element_text(size = 9),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 9),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

plot_coef_ac_tmeanpos_tmeanchange <- ggplot(fitted_slope_ac_tmeanchange_realm) + 
  facet_wrap(~ realm, scales = "free_x") + 
  geom_point(aes(x = change_tempmean, y = estimate_slope, fill = realm, size = 1/se_slope), shape = 21, color = "black", alpha = 0.7) +
  geom_ribbon(aes(x = change_tempmean, ymin = Q2.5, ymax = Q97.5, fill = realm), alpha = 0.3) + 
  geom_line(aes(x = change_tempmean, y = Estimate, color = realm), lwd = 1) + 
  geom_vline(xintercept = 0, linetype =2, colour = "gray", lwd = 0.4) + 
  geom_hline(yintercept = 0, linetype =2, colour = "gray", lwd = 0.4) + 
  geom_text(data = beta_slope_ac_tmeanchange_text, mapping = aes(x= 0.161, y = Inf, label = label), 
            hjust = 1, vjust = 1.3,  parse = T, size = 2.46) + 
  scale_y_continuous(limits = c(-0.12, 0.12), breaks = c( -0.08, 0, 0.08)) + 
  labs(y = "Abundance change ~\n f(baseline thermal position)", x = expression('Temperature change ('*degree*'C/year)')) + 
  scale_colour_manual(name = NULL, values = realm_col) + 
  scale_fill_manual(name = NULL, values = realm_col) + 
  scale_size(range = c(0.3, 2.5)) +
  theme_bw() + 
  theme(legend.position = "n", 
        #strip.text.x  = element_blank(),
        #line = element_line(size = 0.4),
        #strip.background = element_blank(),
        strip.text = element_text(size = 9),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 9),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())


cowplot::plot_grid(plot_coef_oc_tmeanpos_tmeanchange, plot_coef_ac_tmeanpos_tmeanchange,
                   nrow = 2, align = "hv", labels = c("a","b"), label_size = 10)

ggsave(file = "results/Fig.3.png", width = 160, height = 120, units = 'mm', dpi = 300)
ggsave(file = "results/Fig.3.pdf", width = 160, height = 120, units = 'mm')
