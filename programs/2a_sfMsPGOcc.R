# 2a_sfMsPGOcc.R: this script runs a spatially-explicit joint species distribution model
#                   with imperfect detection using sfMsPGOcc from the spOccupancy 
#                   package.This model accounts for (1) imperfect detection,
#                   (2) spatial autocorrelation, and (3) species correlations.
#                   

rm(list = ls())
library(spOccupancy)
library(coda)
library(sf)
library(here)

# Read in the data --------------------------------------------------------
load(here("data", "data_bundles.rdata"))

objs <-  ls(pattern = 'data_')

data_set_names <- gsub('data_', '', objs)

# create list of individual data bundles
data_bundles <- lapply(objs, get)

names(data_bundles) <- data_set_names

data_bundles <- data_bundles[c('dry_rainy_au', 'dry_rainy_ct')]

# Wrap model prep and execution into a single function ----------------------------------------------------------

run_models <- function(data_list, data_set_name, outdir = 'results/spOccupancy_fits/julian_dry_rainy_spatial', n.batch = 6000, batch.length = 25, n.thin = 50, n.chains = 4) {
  
  # create directory for outputs if needed
  suppressWarnings(dir.create(outdir, recursive = T))
  
  # center and scale occupancy and detection covariates
  data_list$occ.covs <- apply(data_list$occ.covs, 2, scale)
  
  data_list$det.covs <- lapply(data_list$det.covs, scale)
  
  # Priors for spatial range parameters
  distances <- dist(data_list$coords)
  mean.dist <- mean(distances)
  min.dist <- min(distances)
  max.dist <- max(distances)
  
  # Priors ------------------------------
    prior.list <- list(beta.comm.normal = list(mean = 0, var = 2.72),
                       alpha.comm.normal = list(mean = 0, var = 2.72),
                       tau.sq.beta.ig = list(a = 0.1, b = 0.1),
                       tau.sq.alpha.ig = list(a = 0.1, b = 0.1), 
                       phi.unif = list(a = 3 / max.dist, b = 3 / min.dist))
    
   # Specify the models -----------------------------------------------------------
    n.neighbors <- 15
    n.factors <- 2
    cov.model <- "exponential"
    n.samples <- n.batch * batch.length
    n.burn <- 0.6*(n.batch*batch.length)
    
#### model with main effects of cover, distance to water, distance to forest edge and basal area

out_cweb <- sfMsPGOcc(occ.formula = ~ cover + dist_water_500 + dist_edge_500 + basal_area,
                         det.formula = ~ shrubs + basal_area,
                         data = data_list, priors = prior.list, n.neighbors = n.neighbors,
                         cov.model = cov.model,
                         NNGP = TRUE, n.factors = n.factors, n.batch = n.batch,
                         batch.length = batch.length, n.burn = n.burn,
                         accept.rate = 0.43, n.thin = n.thin,
                         n.omp.threads = 7, n.chains = n.chains, k.fold = 4,
                         k.fold.threads = 4, k.fold.seed = 73)
    
save(out_cweb, file = paste(outdir, "/", data_set_name, "_cweb-", n.chains, "-chain-", Sys.Date(), ".R", sep = ''))    
    
### model with main effects of cover, distance to water, and distance to forest edge
out_cwe <- sfMsPGOcc(occ.formula = ~ cover + dist_water_500 + dist_edge_500,
                     det.formula = ~ shrubs + basal_area,
                     data = data_list, priors = prior.list, n.neighbors = n.neighbors,
                     cov.model = cov.model,
                     NNGP = TRUE, n.factors = n.factors, n.batch = n.batch,
                     batch.length = batch.length, n.burn = n.burn,
                     accept.rate = 0.43, n.thin = n.thin,
                     n.omp.threads = 7, n.chains = n.chains, k.fold = 4,
                     k.fold.threads = 4, k.fold.seed = 73)

save(out_cwe, file = paste(outdir, "/", data_set_name, "_cwe-", n.chains, "-chain-", Sys.Date(), ".R", sep = ''))

#### model with main effects of cover, distance to water, and basal area
out_cwb <- sfMsPGOcc(occ.formula = ~ cover + dist_water_500 + basal_area,
                     det.formula = ~ shrubs + basal_area,
                     data = data_list, priors = prior.list, n.neighbors = n.neighbors,
                     cov.model = cov.model,
                     NNGP = TRUE, n.factors = n.factors, n.batch = n.batch,
                     batch.length = batch.length, n.burn = n.burn,
                     accept.rate = 0.43, n.thin = n.thin,
                     n.omp.threads = 7, n.chains = n.chains, k.fold = 4,
                     k.fold.threads = 4, k.fold.seed = 73)

save(out_cwb, file = paste(outdir, "/", data_set_name, "_cwb-", n.chains, "-chain-", Sys.Date(), ".R", sep = ''))

#### model with main effects of distance to water, distance to forest edge and basal area
out_web <- sfMsPGOcc(occ.formula = ~ dist_water_500 + dist_edge_500 + basal_area,
                     det.formula = ~ shrubs + basal_area,
                     data = data_list, priors = prior.list, n.neighbors = n.neighbors,
                     cov.model = cov.model,
                     NNGP = TRUE, n.factors = n.factors, n.batch = n.batch,
                     batch.length = batch.length, n.burn = n.burn,
                     accept.rate = 0.43, n.thin = n.thin,
                     n.omp.threads = 7, n.chains = n.chains, k.fold = 4,
                     k.fold.threads = 4, k.fold.seed = 73)

save(out_web, file = paste(outdir, "/", data_set_name, "_web-", n.chains, "-chain-", Sys.Date(), ".R", sep = ''))

#### model with main effects of cover, distance to forest edge and basal area
out_ceb <- sfMsPGOcc(occ.formula = ~ cover + dist_edge_500 + basal_area,
                     det.formula = ~ shrubs + basal_area,
                     data = data_list, priors = prior.list, n.neighbors = n.neighbors,
                     cov.model = cov.model,
                     NNGP = TRUE, n.factors = n.factors, n.batch = n.batch,
                     batch.length = batch.length, n.burn = n.burn,
                     accept.rate = 0.43, n.thin = n.thin,
                     n.omp.threads = 7, n.chains = n.chains, k.fold = 4,
                     k.fold.threads = 4, k.fold.seed = 73)

save(out_ceb, file = paste(outdir, "/", data_set_name, "_ceb-", n.chains, "-chain-", Sys.Date(), ".R", sep = ''))
   
### model with main effects of cover and distance to water
out_cw <- sfMsPGOcc(occ.formula = ~ cover + dist_water_500,
                    det.formula = ~ shrubs + basal_area,
                    data = data_list, priors = prior.list, n.neighbors = n.neighbors,
                    cov.model = cov.model,
                    NNGP = TRUE, n.factors = n.factors, n.batch = n.batch,
                    batch.length = batch.length, n.burn = n.burn,
                    accept.rate = 0.43, n.thin = n.thin,
                    n.omp.threads = 7, n.chains = n.chains, k.fold = 4,
                    k.fold.threads = 4, k.fold.seed = 73)

save(out_cw, file = paste(outdir, "/", data_set_name, "_cw-", n.chains, "-chain-", Sys.Date(), ".R", sep = ''))
    
### model with main effects of cover and distance to forest edge
out_ce <- sfMsPGOcc(occ.formula = ~ cover + dist_edge_500,
                    det.formula = ~ shrubs + basal_area,
                    data = data_list, priors = prior.list, n.neighbors = n.neighbors,
                    cov.model = cov.model,
                    NNGP = TRUE, n.factors = n.factors, n.batch = n.batch,
                    batch.length = batch.length, n.burn = n.burn,
                    accept.rate = 0.43, n.thin = n.thin,
                    n.omp.threads = 7, n.chains = n.chains, k.fold = 4,
                    k.fold.threads = 4, k.fold.seed = 73)

save(out_ce, file = paste(outdir, "/", data_set_name, "_ce-", n.chains, "-chain-", Sys.Date(), ".R", sep = ''))

#### model with main effects of cover and basal area
out_cb <- sfMsPGOcc(occ.formula = ~ cover + basal_area,
                    det.formula = ~ shrubs + basal_area,
                    data = data_list, priors = prior.list, n.neighbors = n.neighbors,
                    cov.model = cov.model,
                    NNGP = TRUE, n.factors = n.factors, n.batch = n.batch,
                    batch.length = batch.length, n.burn = n.burn,
                    accept.rate = 0.43, n.thin = n.thin,
                    n.omp.threads = 7, n.chains = n.chains, k.fold = 4,
                    k.fold.threads = 4, k.fold.seed = 73)

save(out_cb, file = paste(outdir, "/", data_set_name, "_cb-", n.chains, "-chain-", Sys.Date(), ".R", sep = ''))

### model with main effects of distance to water and distance to forest edge
out_we <- sfMsPGOcc(occ.formula = ~ dist_water_500 + dist_edge_500,
                      det.formula = ~ shrubs + basal_area,
                      data = data_list, priors = prior.list, n.neighbors = n.neighbors,
                      cov.model = cov.model, NNGP = TRUE,
                      n.factors = n.factors, n.batch = n.batch, batch.length = batch.length,
                      n.burn = n.burn, accept.rate = 0.43, n.thin = n.thin, n.omp.threads = 7,
                      n.chains = n.chains, k.fold = 4, k.fold.threads = 4, k.fold.seed = 73)

save(out_we, file = paste(outdir, "/", data_set_name, "_we-", n.chains, "-chain-", Sys.Date(), ".R", sep = ''))
  
#### model with main effects of distance to water and basal area
out_wb <- sfMsPGOcc(occ.formula = ~ dist_water_500 + basal_area,
                      det.formula = ~ shrubs + basal_area,
                      data = data_list, priors = prior.list, n.neighbors = n.neighbors,
                      cov.model = cov.model, NNGP = TRUE,
                      n.factors = n.factors, n.batch = n.batch, batch.length = batch.length,
                      n.burn = n.burn, accept.rate = 0.43, n.thin = n.thin, n.omp.threads = 7,
                      n.chains = n.chains, k.fold = 4, k.fold.threads = 4, k.fold.seed = 73)
  
save(out_wb, file = paste(outdir, "/", data_set_name, "_wb-", n.chains, "-chain-", Sys.Date(), ".R", sep = ''))

#### model with main effects of distance to forest edge and basal area
out_eb <- sfMsPGOcc(occ.formula = ~ dist_edge_500 + basal_area,
                    det.formula = ~ shrubs + basal_area,
                    data = data_list, priors = prior.list, n.neighbors = n.neighbors,
                    cov.model = cov.model, NNGP = TRUE,
                    n.factors = n.factors, n.batch = n.batch, batch.length = batch.length,
                    n.burn = n.burn, accept.rate = 0.43, n.thin = n.thin, n.omp.threads = 7,
                    n.chains = n.chains, k.fold = 4, k.fold.threads = 4, k.fold.seed = 73)

save(out_eb, file = paste(outdir, "/", data_set_name, "_eb-", n.chains, "-chain-", Sys.Date(), ".R", sep = ''))

### null model (intercept only)
out_null <- sfMsPGOcc(occ.formula = ~ 1,
                      det.formula = ~ 1,
                      data = data_list, priors = prior.list, n.neighbors = n.neighbors,
                      cov.model = cov.model, NNGP = TRUE,
                      n.factors = n.factors, n.batch = n.batch,
                      batch.length = batch.length, n.burn = n.burn,
                      accept.rate = 0.43, n.thin = n.thin,
                      n.omp.threads = 7, n.chains = n.chains,
                      k.fold = 4, k.fold.threads = 4, k.fold.seed = 73)

save(out_null, file = paste(outdir, "/", data_set_name, "_null-", n.chains, "-chain-", Sys.Date(), ".R", sep = ''))
}



# Run the models ----------------------------------------------------------

sapply(names(data_bundles), function(x) run_models(data_bundles[[x]], x,  outdir = 'results/spOccupancy_fits/julian_dry_rainy_spatial'), simplify = F)

