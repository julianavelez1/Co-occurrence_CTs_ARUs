# 5_predictions.R: this script gets model predictions (supplementary information).

rm(list = ls())

library(spOccupancy)
library(here)
library(tidyverse)

set.seed(7337)

load(here('data', 'data_bundles.rdata'))

rm(list=setdiff(ls(), c("data_dry_rainy_au", "data_dry_rainy_ct")))

covs.list <- data_dry_rainy_au

occ_covs_means <- apply(covs.list$occ.covs, 2, mean)
occ_covs_sd <- apply(covs.list$occ.covs, 2, sd)

det_covs_means <- sapply(covs.list$det.covs, mean)
det_covs_sd <- sapply(covs.list$det.covs, sd)

sp.codes <- c("Cattle",
              "Brocket Deer",
              "Collared peccary", 
              "Lowland tapir","Poaching",
              "White-lipped peccary", "White-tailed deer")


# Function for predictons -------------------------------------------------


get_preds <- function(path_to_model = "/Volumes/NO NAME/dry_rainy_au_we-4-chain-2022-10-19.R", covariate = c('water', 'edge', 'basal_area')) {
  
  match.arg(covariate)
  tmp <- load(path_to_model)
  model_name <- 'out.2.sfMsPGOcc'
  
  assign(model_name, get(tmp)) 
  
  # need coordinates for each prediction location, we're using the centroid of the 52 sites
  coords.0 <- apply(out.2.sfMsPGOcc$coords, 2, mean)
  
  # c_range <- range(scale(covs.list$occ.covs$cover))
  w_range <- range(scale(covs.list$occ.covs$dist_water_500))
  e_range <- range(scale(covs.list$occ.covs$dist_edge_500))
  b_range <- range(scale(covs.list$occ.covs$basal_area))
  
  # cover <- seq(c_range[1], c_range[2], length.out = 100)
  dist_water_500 <- seq(w_range[1], w_range[2], length.out = 100)
  dist_edge_500 <- seq(e_range[1], e_range[2], length.out = 100)
  basal_area <- seq(b_range[1], b_range[2], length.out = 100)

  
  if (covariate == 'water') {
    X.0 <- cbind(1, dist_water_500, dist_edge_500 = 0, basal_area = 0)
    covariate_unscaled <- ((dist_water_500 *occ_covs_sd['dist_water_500']) + occ_covs_means['dist_water_500']) / 1000
  }
  
  if (covariate == 'edge') {
    X.0 <- cbind(1, dist_water_500 = 0, dist_edge_500, basal_area = 0)
    covariate_unscaled <- ((dist_edge_500 *occ_covs_sd['dist_edge_500']) + occ_covs_means['dist_edge_500']) / 1000
  }
  
  if (covariate == 'basal_area') {
    X.0 <- cbind(1, dist_water_500 = 0, dist_edge_500 = 0, basal_area)
    covariate_unscaled <- ((basal_area *occ_covs_sd['basal_area']) + occ_covs_means['basal_area']) / 1000
  }
  
  coords.0 <- cbind(rep(coords.0[1] ,nrow(X.0)), rep(coords.0[2],nrow(X.0)))
  
  test <- predict(out.2.sfMsPGOcc, X.0, coords.0, ignore.RE = T)
  
  colnames(test$psi.0.samples) <- sp.codes
  
  test_mean <- t(apply(test$psi.0.samples, 2:3, mean))
  test_q0.025 <- t(apply(test$psi.0.samples, 2:3, quantile, 0.025))
  test_q0.975 <- t(apply(test$psi.0.samples, 2:3, quantile, 0.975))

  pred_summary <- do.call(cbind, list(data.frame(covariate = covariate_unscaled, test_mean) %>% gather(Species, psi_mean, -covariate), data.frame(test_q0.025) %>% gather(Species, psi_lwr) %>% dplyr::select(-Species), data.frame(test_q0.975) %>% gather(Species, psi_upr) %>% dplyr::select(-Species)))
  
return(pred_summary)
  
}


# CTs-ARUs ----------------------------------------------------------------

## Distance to water

full_preds_dry_rainy_au <- get_preds(here("results", "spOccupancy_fits",
                                          "julian_dry_rainy_spatial",
                                          "dry_rainy_au_web-4-chain-2023-07-24.R"),
                                     covariate = 'water') %>% 
  mutate(Species = str_replace_all(Species, "\\.", " "))

water_au <- full_preds_dry_rainy_au %>% 
  ggplot(aes(covariate, psi_mean, group = Species)) + 
  geom_ribbon(aes(ymin = psi_lwr, ymax = psi_upr, fill = Species), 
              alpha = 0.2) + 
  geom_line(aes(colour = Species)) + 
  labs(x = 'Distance to water (km)', y = expression(psi)) + 
  facet_wrap(~Species) + scale_colour_viridis_d() + 
  scale_fill_viridis_d() + 
  theme(panel.grid.major = element_blank(),
                                 panel.grid.minor = element_blank(),
                                 panel.background = element_blank(),
        legend.position = "none")

## Distance to forest edge

full_preds_dry_rainy_au <- get_preds(here("results", "spOccupancy_fits",
                                          "julian_dry_rainy_spatial",
                                          "dry_rainy_au_web-4-chain-2023-07-24.R"),
                                     covariate = 'edge') %>% mutate(Species = str_replace_all(Species, "\\.", " "))

edge_au <- full_preds_dry_rainy_au %>% 
  ggplot(aes(covariate, psi_mean, group = Species)) + 
  geom_ribbon(aes(ymin = psi_lwr, ymax = psi_upr, fill = Species), 
              alpha = 0.2) + 
  geom_line(aes(colour = Species)) + 
  labs(x = 'Distance to edge (km)', y = expression(psi)) + 
  facet_wrap(~Species) + scale_colour_viridis_d() + 
  scale_fill_viridis_d() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "none")

full_preds_dry_rainy_au <- get_preds(here("results", "spOccupancy_fits",
                                          "julian_dry_rainy_spatial",
                                          "dry_rainy_au_web-4-chain-2023-07-24.R"),
                                     covariate = 'basal_area') %>% mutate(Species = str_replace_all(Species, "\\.", " "))

## Basal area

basal_area_au <- full_preds_dry_rainy_au %>% 
  ggplot(aes(covariate, psi_mean, group = Species)) + 
  geom_ribbon(aes(ymin = psi_lwr, ymax = psi_upr, fill = Species), 
              alpha = 0.2) + 
  geom_line(aes(colour = Species)) + 
  labs(x = 'Basal area of woody vegetation', y = expression(psi)) + 
  facet_wrap(~Species) + scale_colour_viridis_d() + 
  scale_fill_viridis_d() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "none")


# CTs ---------------------------------------------------------------------

## Distance to water

full_preds_dry_rainy_ct <- get_preds(here("results", "spOccupancy_fits",
                                          "julian_dry_rainy_spatial", 
                                          "dry_rainy_ct_web-4-chain-2023-07-25.R"), 
                                          covariate = 'water') %>% 
  mutate(Species = str_replace_all(Species, "\\.", " "))

water_ct <- full_preds_dry_rainy_ct %>% 
  ggplot(aes(covariate, psi_mean, group = Species)) + 
  geom_ribbon(aes(ymin = psi_lwr, ymax = psi_upr, fill = Species), 
              alpha = 0.2) + 
  geom_line(aes(colour = Species)) + 
  labs(x = 'Distance to water (km)', y = expression(psi)) + 
  facet_wrap(~Species) + scale_colour_viridis_d() + 
  scale_fill_viridis_d() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "none")

## Distance to forest edge

full_preds_dry_rainy_ct <- get_preds(here("results", "spOccupancy_fits",
                                          "julian_dry_rainy_spatial", 
                                          "dry_rainy_ct_web-4-chain-2023-07-25.R"),
                                     covariate = 'edge') %>% mutate(Species = str_replace_all(Species, "\\.", " "))

edge_ct <- full_preds_dry_rainy_ct %>% 
  ggplot(aes(covariate, psi_mean, group = Species)) + 
  geom_ribbon(aes(ymin = psi_lwr, ymax = psi_upr, fill = Species), 
              alpha = 0.2) + 
  geom_line(aes(colour = Species)) + 
  labs(x = 'Distance to edge (km)', y = expression(psi)) + 
  facet_wrap(~Species) + scale_colour_viridis_d() + 
  scale_fill_viridis_d() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "none")

## Basal area

full_preds_dry_rainy_ct <- get_preds(here("results", "spOccupancy_fits",
                                          "julian_dry_rainy_spatial", 
                                          "dry_rainy_ct_web-4-chain-2023-07-25.R"),
                                     covariate = 'basal_area') %>% 
  mutate(Species = str_replace_all(Species, "\\.", " "))

basal_area_ct <- full_preds_dry_rainy_ct %>% 
  ggplot(aes(covariate, psi_mean, group = Species)) + 
  geom_ribbon(aes(ymin = psi_lwr, ymax = psi_upr, fill = Species), 
              alpha = 0.2) + 
  geom_line(aes(colour = Species)) + 
  labs(x = 'Basal area of woody vegetation', y = expression(psi)) + 
  facet_wrap(~Species) + scale_colour_viridis_d() + 
  scale_fill_viridis_d() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "none")

# Save data ---------------------------------------------------------------

# CTs-ARUs
saveRDS(water_au, "figures/marginal_water_au.RDS")
saveRDS(edge_au, "figures/marginal_edge_au.RDS")
saveRDS(basal_area_au, "figures/marginal_basal_area_au.RDS")


# CTs
saveRDS(water_ct, "figures/marginal_water_ct.RDS")
saveRDS(edge_ct, "figures/marginal_edge_ct.RDS")
saveRDS(basal_area_ct, "figures/marginal_basal_area_ct.RDS")


