# 3_waic_convergence.R: this script obtains the Widely Applicable Information 
#                       Criterion and selects the best-performing model.
#                       It also performs convergence assessment by inspecting
#                       traceplots, Rhat values and the effective sample size (ESS).

rm(list = ls())

library(spOccupancy)
library(tidyverse)
library(here)


# Read model fits ---------------------------------------------------------

path_to_mods <- 'results/spOccupancy_fits/'
fits <- list.files(path_to_mods, pattern = 'dry_rainy*', recursive = TRUE)

read_fits <- function(path = fits){
  tmp <- load(paste0(path_to_mods, path))
  objname <- strsplit(path, '-4-chain')[[1]][1]
  assign(objname, get(tmp))
  out <- list(get(objname))
  names(out) <- objname
  return(out)
}

fits <- lapply(fits, read_fits)
fits <- unlist(fits, recursive = F)

# Split fits names --------------------------------------------------------

model_name_components <- strsplit(names(fits), '_')
dataset <- sapply(model_name_components, function(x) paste(x[1:(length(x) - 1)], 
                                                           collapse = '_'))
model <- sapply(model_name_components, function(x) x[length(x)])


# Apply WAICOcc to model fits and select best model ------------------------------------------------

waic_df <- fits %>% 
  map_df(waicOcc) %>% 
  mutate(dataset_model = names(fits),
         dataset = dataset,
         model = model,
         dataset = str_remove(dataset, "julian_dry_rainy_"))


# Get best model WAIC by data set and type of model (i.e., spatial or nonspatial)

remove <- c("nonspatial/", "spatial/")

best_mods <- waic_df %>% 
  mutate(dataset = str_remove(dataset, paste(remove, collapse = "|"))) %>% 
  group_by(dataset) %>% 
  slice_min(WAIC, n=3) # filter the best 3 models

best_mods <- fits[which(names(fits) %in% best_mods$dataset_model)]

best_mods <- names(best_mods) %>% 
  str_detect('web') %>% 
  keep(best_mods, .)


# Inspect traceplots, r-hat values and effective sample sizes of the posterior samples.

# CTs

best_mod_ct <- names(best_mods) %>% 
  str_detect('ct') %>% 
  keep(best_mods, .)

summary(best_mod_ct[[1]])

par(mar = c(2, 2, 2, 2))
plot(best_mod_ct[[1]]$beta.comm.samples, density = FALSE)
plot(best_mod_ct[[1]]$beta.samples, density = FALSE)
plot(best_mod_ct[[1]]$alpha.samples, density = FALSE) 
plot(best_mod_ct[[1]]$lambda.samples, density = FALSE) 

# CTs-ARUs

best_mod_au <- names(best_mods) %>% 
  str_detect('au') %>% 
  keep(best_mods, .)

summary(best_mod_au[[1]])

plot(best_mod_au[[1]]$beta.comm.samples, density = FALSE)
plot(best_mod_au[[1]]$beta.samples, density = FALSE)
plot(best_mod_au[[1]]$alpha.samples, density = FALSE) 
plot(best_mod_au[[1]]$lambda.samples, density = FALSE) 

# Save WAIC and best models -----------------------------------------------

#saveRDS(waic_df, file = here("results", "waic", "waic_df.RDS"))
#saveRDS(best_mods, file = here("results", "best_mods", "fits_best_models.RDS"))

