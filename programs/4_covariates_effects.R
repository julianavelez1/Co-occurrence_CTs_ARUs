# 4_covariates_effects.R: this script plots point estimates and credible intervals        


rm(list=ls()) 

library(posterior)
library(tidybayes)
library(spOccupancy)
library(tidyverse)
library(here)

# Read best models --------------------------------------------------------

best_mods <- readRDS(here("results", "best_mods", "fits_best_models.RDS"))

summary(best_mods$`julian_dry_rainy_spatial/dry_rainy_au_web`)


# Species names -----------------------------------------------------------

sp.names <- c("Cattle","Brocket deer", "Collared peccary", 
              "Lowland tapir","Poaching",
              "White-lipped peccary", "White-tailed deer")

# Occupancy estimates ---------------------------------------------------------------

occ_ct <- best_mods$`julian_dry_rainy_spatial/dry_rainy_ct_web`$psi.samples
mean_occ_ct <- apply(occ_ct, 1:2, mean)

mean_occ_ct <-  data.frame(mean_occ_ct)
names(mean_occ_ct) <- sp.names
mean_occ_ct <- mean_occ_ct %>% mutate("Method" = "CTs")

occ_au <- best_mods$`julian_dry_rainy_spatial/dry_rainy_au_web`$psi.samples
mean_occ_au <- apply(occ_au, 1:2, mean)

mean_occ_au <-  data.frame(mean_occ_au) 
names(mean_occ_au) <- sp.names
mean_occ_au <- mean_occ_au %>% mutate("Method" = "ARUs")

mean_occ_both <- rbind(mean_occ_au, mean_occ_ct)

# Get occupancy values 

ci_occ <- mean_occ_both %>% 
  gather(key = "Species", value = "Occupancy", -Method) %>% 
  group_by(Species, Method) %>% median_qi(Occupancy, .width = 0.89)



mean_occ_p <- mean_occ_both %>% 
  gather(key = "Species", value = "Occupancy", -Method) %>% 
  mutate(Species = factor(Species, levels = c("Cattle", "Poaching",
                                              "Brocket deer",
                                              "Collared peccary", 
                                              "Lowland tapir",
                                              "White-lipped peccary", 
                                              "White-tailed deer"
  ))) %>% 
  filter(!(Species %in% c("Brocket deer",
                          "Collared peccary", 
                          "Lowland tapir",
                          "White-lipped peccary", 
                          "White-tailed deer") &
             Method %in% "ARUs")) %>%
  ggplot(aes(Species, Occupancy, fill = Method, shape = Method, colour = Species)) +
  stat_pointinterval(position = "dodge",
                     point_interval = "median_qi",
                     .width = c(0.5, 0.89), point_size = 3.5) +
  scale_shape_manual(values = c(8,16)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size = rel(1.3)),
        axis.title = element_text(size = rel(1.5)),
        legend.position = "none"
        ) +
  coord_flip() +
  labs(y = expression("Occupancy probability ("*psi*")")) +
  scale_color_viridis_d(guide = "none")


# Occupancy betas posterior distributions ---------------------------------------------


beta_samples <- best_mods %>% 
  map(~.x %>% 
        pluck("beta.samples"))


## Distance to water -------------------------------------------------------

occ_posterior_water <- beta_samples %>% 
  map(~.x %>% 
        as_draws_df() %>% 
        select(contains(("water"))) %>% 
        rename_all(~sp.names) %>% 
        gather(species, value) %>%
        rename(Species = species, 
               Estimates = value)
  ) %>% 
  map_at(~.x %>%
           mutate(Method = "ARUs"), .at = "julian_dry_rainy_spatial/dry_rainy_au_web") %>% 
  map_at(~.x %>% 
           mutate(Method = "CTs"), .at = "julian_dry_rainy_spatial/dry_rainy_ct_web") %>% 
  bind_rows() %>% 
  mutate_if(is.numeric, ~round(., 1))

ci_water <- occ_posterior_water %>% 
  group_by(Species, Method) %>% 
  median_qi(Estimates, .width = 0.89)


disturbance <- c("Cattle", "Poaching")

ybreaks_occ <- seq(-4, 4, 2)
ylims_occ <- c(-4.6, 4.1)

# Plot

p_water_wild <- occ_posterior_water %>% 
  filter(!Species %in% disturbance) %>%
  filter(Method %in% "CTs") %>%
  ggplot(aes(Species, Estimates, fill = Method, colour = Species)) +
  stat_pointinterval(position = "dodge",
                     point_interval = "median_qi",
                     .width = c(0.5, 0.89),
                     point_size = 3.5,
                     shape = 16)  +
  geom_hline(aes(yintercept = 0), linetype = 2) +
  scale_y_continuous(breaks = ybreaks_occ) +
  coord_flip(ylim = ylims_occ) +# flip axes
  theme(axis.title.y=element_blank(),
        axis.text = element_text(size=16),
        axis.title = element_text(size=18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 18)) +
  ylab('Distance to water') + 
  ggtitle("Occupancy probability") +
  scale_color_viridis_d()

p_water_dist <- occ_posterior_water %>% 
  filter(Species %in% disturbance) %>%
  ggplot(aes(Species, Estimates, 
             fill = Method, shape = Method, colour = Species)) +
  stat_pointinterval(position = "dodge",
                     point_interval = "median_qi",
                     .width = c(0.5, 0.89),
                     point_size = 3.5)  +
  scale_shape_manual(values = c(8,16)) +
  geom_hline(aes(yintercept = 0), linetype = 2)  +
  scale_y_continuous(breaks = ybreaks_occ) +
  coord_flip(ylim = ylims_occ)  +# flip axes
  theme(axis.title.y=element_blank(),
        axis.text = element_text(size=16),
        axis.title = element_text(size=18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 18)) +
  ylab('Distance to water') + 
  ggtitle("Occupancy probability") +
  scale_color_viridis_d(begin = 0, end = 0.5, guide = "none")


## Distance to forest edge -------------------------------------------------


occ_posterior_edge <- beta_samples %>% 
  map(~.x %>% 
        as_draws_df() %>% 
        select(contains(("edge"))) %>% 
        rename_all(~sp.names) %>% 
        gather(species, value) %>%
        rename(Species = species, 
               Estimates = value)
  ) %>% 
  map_at(~.x %>%
           mutate(Method = "ARUs"), .at = "julian_dry_rainy_spatial/dry_rainy_au_web") %>% 
  map_at(~.x %>% 
           mutate(Method = "CTs"), .at = "julian_dry_rainy_spatial/dry_rainy_ct_web") %>% 
  bind_rows() %>% 
  mutate_if(is.numeric, ~round(., 1))

ci_edge <- occ_posterior_edge %>% 
  group_by(Species, Method) %>% median_qi(Estimates, .width = 0.89)

# Plot

p_edge_wild <- occ_posterior_edge %>% 
  filter(!Species %in% disturbance) %>%
  filter(Method %in% "CTs") %>%
  ggplot(aes(Species, Estimates, fill = Method, colour = Species)) +
  stat_pointinterval(position = "dodge",
                     point_interval = "median_qi",
                     .width = c(0.5, 0.89),
                     point_size = 3.5,
                     shape = 16)  +
  geom_hline(aes(yintercept = 0), linetype = 2) +
  scale_y_continuous(breaks = ybreaks_occ) +
  coord_flip(ylim = ylims_occ) +# flip axes
  theme(axis.title.y=element_blank(),
        axis.text = element_text(size=16),
        axis.title = element_text(size=18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 18)) +
  ylab('Distance to forest edge') + 
  ggtitle("Occupancy probability") +
  scale_color_viridis_d()

p_edge_dist <- occ_posterior_edge %>% 
  filter(Species %in% disturbance) %>%
  ggplot(aes(Species, Estimates, fill = Method, shape = Method, colour = Species)) +
  stat_pointinterval(position = "dodge",
                     point_interval = "median_qi",
                     .width = c(0.5, 0.89),
                     point_size = 3.5)  +
  scale_shape_manual(values = c(8, 16)) +
  geom_hline(aes(yintercept = 0), linetype = 2)  +
  scale_y_continuous(breaks = ybreaks_occ) +
  coord_flip(ylim = ylims_occ)  + # flip axes
  theme(axis.title.y=element_blank(),
        axis.text = element_text(size=16),
        axis.title = element_text(size=18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 18)) +
  ylab('Distance to forest edge') + 
  ggtitle("Occupancy probability") +
  scale_color_viridis_d(begin = 0, end = 0.5, guide = "none")

## Basal area -------------------------------------------------


occ_posterior_ba <- beta_samples %>% 
  map(~.x %>% 
        as_draws_df() %>% 
        select(contains(("basal_area"))) %>% 
        rename_all(~sp.names) %>% 
        gather(species, value) %>%
        rename(Species = species, 
               Estimates = value)
  ) %>% 
  map_at(~.x %>%
           mutate(Method = "ARUs"), .at = "julian_dry_rainy_spatial/dry_rainy_au_web") %>% 
  map_at(~.x %>% 
           mutate(Method = "CTs"), .at = "julian_dry_rainy_spatial/dry_rainy_ct_web") %>% 
  bind_rows() %>% 
  mutate_if(is.numeric, ~round(., 1))

ci_ba_occ <- occ_posterior_ba %>% 
  group_by(Species, Method) %>% median_qi(Estimates, .width = 0.89)

# Plot

p_ba_wild_occ <- occ_posterior_ba %>% 
  filter(!Species %in% disturbance) %>%
  filter(Method %in% "CTs") %>%
  ggplot(aes(Species, Estimates, fill = Method, colour = Species)) +
  stat_pointinterval(position = "dodge",
                     point_interval = "median_qi",
                     .width = c(0.5, 0.89),
                     point_size = 3.5,
                     shape = 16)  +
  geom_hline(aes(yintercept = 0), linetype = 2) +
  scale_y_continuous(breaks = ybreaks_occ) +
  coord_flip(ylim = ylims_occ) +# flip axes
  theme(axis.title.y=element_blank(),
        axis.text = element_text(size=16),
        axis.title = element_text(size=18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 18)) +
  ylab('Basal area') + 
  ggtitle("Occupancy probability") +
  scale_color_viridis_d()

p_ba_dist_occ <- occ_posterior_ba %>% 
  filter(Species %in% disturbance) %>%
  ggplot(aes(Species, Estimates, fill = Method, shape = Method, colour = Species)) +
  stat_pointinterval(position = "dodge",
                     point_interval = "median_qi",
                     .width = c(0.5, 0.89),
                     point_size = 3.5)  +
  scale_shape_manual(values = c(8, 16)) +
  geom_hline(aes(yintercept = 0), linetype = 2)  +
  scale_y_continuous(breaks = ybreaks_occ) +
  coord_flip(ylim = ylims_occ)  + # flip axes
  theme(axis.title.y=element_blank(),
        axis.text = element_text(size=16),
        axis.title = element_text(size=18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 18)) +
  ylab('Basal area') + 
  ggtitle("Occupancy probability") +
  scale_color_viridis_d(begin = 0, end = 0.5, guide = "none")


# Detection alphas posterior distribution ----------------------------------------

alpha_samples <- best_mods %>% 
  map(~.x %>% 
        pluck("alpha.samples"))


## Shrub cover -------------------------------------------------------------

det_posterior_shrubs <- alpha_samples %>% 
  map(~.x %>% 
        as_draws_df() %>% 
        select(contains(("shrubs"))) %>% 
        rename_all(~sp.names) %>% 
        gather(species, value) %>%
        rename(Species = species, 
               Estimates = value)
  ) %>% 
  map_at(~.x %>%
           mutate(Method = "ARUs"), .at = "julian_dry_rainy_spatial/dry_rainy_au_web") %>% 
  map_at(~.x %>% 
           mutate(Method = "CTs"), .at = "julian_dry_rainy_spatial/dry_rainy_ct_web") %>% 
  bind_rows() %>% 
  mutate_if(is.numeric, ~round(., 1))

ci_shrubs <- det_posterior_shrubs %>% 
  group_by(Species, Method) %>% median_qi(Estimates, .width = 0.89)

# Plot

ybreaks_det <- seq(-1, 1, 0.5)
ylims_det <- c(-1, 1)

p_shrubs_wild <- det_posterior_shrubs %>% 
  filter(!Species %in% disturbance) %>%
  filter(Method %in% "CTs") %>% 
  ggplot(aes(Species, Estimates, 
             fill = Method, shape = Method, colour = Species)) +
  stat_pointinterval(position = "dodge",
                     point_interval = "median_qi",
                     .width = c(0.5, 0.89),
                     point_size = 3.5,
                     shape = 16)  +
  geom_hline(aes(yintercept = 0), linetype = 2) +
  scale_y_continuous(breaks = ybreaks_det) +
  coord_flip(ylim = ylims_det) +# flip axes
  theme(axis.title.y=element_blank(),
        axis.text = element_text(size=16),
        axis.title = element_text(size=18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 18)) +
  ylab('Shrubs aereal coverage') + 
  ggtitle("Detection probability") +
  scale_color_viridis_d()


p_shrubs_dist <- det_posterior_shrubs %>% 
  filter(Species %in% disturbance) %>%
  ggplot(aes(Species, Estimates, fill = Method, shape = Method, colour = Species)) +
  stat_pointinterval(position = "dodge",
                     point_interval = "median_qi",
                     .width = c(0.5, 0.89),
                     point_size = 3.5)  +
  scale_shape_manual(values = c(8,16)) +
  geom_hline(aes(yintercept = 0), linetype = 2)  +
  scale_y_continuous(breaks = ybreaks_det) +
  coord_flip(ylim = ylims_det)  +# flip axes
  theme(axis.title.y=element_blank(),
        axis.text = element_text(size=16),
        axis.title = element_text(size=18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 18)) +
  ylab('Shrubs aereal coverage') + 
  ggtitle("Detection probability") +
  scale_color_viridis_d(begin = 0, end = 0.5, guide = "none")


## Basal area --------------------------------------------------------------


det_posterior_ba <- alpha_samples %>% 
  map(~.x %>% 
        as_draws_df() %>% 
        select(contains(("basal_area"))) %>% 
        rename_all(~sp.names) %>% 
        gather(species, value) %>%
        rename(Species = species, 
               Estimates = value)) %>% 
  map_at(~.x %>%
           mutate(Method = "ARUs"), .at = "julian_dry_rainy_spatial/dry_rainy_au_web") %>% 
  map_at(~.x %>% 
           mutate(Method = "CTs"), .at = "julian_dry_rainy_spatial/dry_rainy_ct_web") %>% 
  bind_rows() %>% 
  mutate_if(is.numeric, ~round(., 1))

ci_ba_det <- det_posterior_ba %>% 
  group_by(Species, Method) %>% median_qi(Estimates, .width = 0.89)

# Plot


p_ba_wild_det <- det_posterior_ba %>% 
  filter(!Species %in% disturbance) %>%
  filter(Method %in% "CTs") %>% 
  ggplot(aes(Species, Estimates, 
             fill = Method, shape = Method, colour = Species)) +
  stat_pointinterval(position = "dodge",
                     point_interval = "median_qi",
                     .width = c(0.5, 0.89),
                     point_size = 3.5,
                     shape = 16)  +
  geom_hline(aes(yintercept = 0), linetype = 2) +
  scale_y_continuous(breaks = ybreaks_det) +
  coord_flip(ylim = ylims_det) +# flip axes
  theme(axis.title.y=element_blank(),
        axis.text = element_text(size=16),
        axis.title = element_text(size=18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 18)) +
  ylab('Basal area') + 
  ggtitle("Detection probability") +
  scale_color_viridis_d()

p_ba_dist_det <- det_posterior_ba %>% 
  filter(Species %in% disturbance) %>% 
  ggplot(aes(Species, Estimates, fill = Method, shape = Method, colour = Species)) +
  stat_pointinterval(position = "dodge",
                     point_interval = "median_qi",
                     .width = c(0.5, 0.89),
                     point_size = 3.5)  +
  scale_shape_manual(values = c(8,16)) +
  geom_hline(aes(yintercept = 0), linetype = 2)  +
  scale_y_continuous(breaks = ybreaks_det) +
  coord_flip(ylim = ylims_det)  +# flip axes
  theme(axis.title.y=element_blank(),
        axis.text = element_text(size=16),
        axis.title = element_text(size=18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 18)) +
  ylab('Basal area') + 
  ggtitle("Detection probability") +
  scale_color_viridis_d(begin = 0, end = 0.5, guide = "none")


# Tabulate estimates ------------------------------------------------------

# Function to get credible intervals
get_CI <- function(mcmc_obj) {
  quants <- summary(mcmc_obj, quantiles = c(0.055, 0.5, 0.945)) # 89% CI
  return(quants$quantiles)
}

# Betas table - example for CTs-ARUs data set
occ_betas <- get_CI(best_mods[[1]]$beta.samples) %>% # species-specific betas
  as.data.frame() %>% 
  rownames_to_column() %>% 
  separate(rowname, sep = "-", into = c("Parameter", "Species")) %>%
  gather(Quantile, Value, -Parameter, -Species) %>% 
  mutate(Value = signif(Value, 2)) %>% 
  spread(Quantile, Value) %>% 
  mutate(content = paste0(`50%`, ' (', `5.5%`,',', `94.5%`,')')) %>% 
  dplyr::select(Species, Parameter, content) %>% 
  spread(Species, content)


# Save output -------------------------------------------------------------

# Estimates
saveRDS(ci_occ, "results/estimates/ci_occ.RDS")
saveRDS(ci_water, "results/estimates/ci_water.RDS")
saveRDS(ci_edge, "results/estimates/ci_edge.RDS")
saveRDS(ci_ba_occ, "results/estimates/ci_ba_occ.RDS")

saveRDS(ci_shrubs, "results/estimates/ci_shrubs.RDS")
saveRDS(ci_ba_det, "results/estimates/ci_ba_det.RDS")


# Plots
saveRDS(mean_occ_p, "figures/mean_occ_p.RDS")

saveRDS(p_water_wild, "figures/occ_water_wild.RDS")
saveRDS(p_edge_wild, "figures/occ_edge_wild.RDS")
saveRDS(p_ba_wild_occ, "figures/occ_ba_wild.RDS")

saveRDS(p_water_dist, "figures/occ_water_disturb.RDS")
saveRDS(p_edge_dist, "figures/occ_edge_disturb.RDS")
saveRDS(p_ba_dist_occ, "figures/occ_ba_disturb.RDS")


saveRDS(p_shrubs_wild, "figures/det_shrubs_wild.RDS")
saveRDS(p_ba_wild_det, "figures/det_ba_wild.RDS")

saveRDS(p_shrubs_dist, "figures/det_shrubs_dist.RDS")
saveRDS(p_ba_dist_det, "figures/det_ba_dist.RDS")
