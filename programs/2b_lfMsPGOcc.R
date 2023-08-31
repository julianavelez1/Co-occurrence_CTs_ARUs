# 2b_lfMsPGOcc.R: this script runs a joint species distribution model
#                   with imperfect detection using lfMsPGOcc from the spOccupancy 
#                   package. This model accounts for (1) imperfect detection;
#                   and (2) residual species correlations.
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

#data_bundles <- data_bundles[c('rainy_ct', 'rainy_au')]

# Wrap model prep and execution into a single function ----------------------------------------------------------

run_models <- function(data_list, data_set_name, outdir = 'results/spOccupancy_fits/julian_dry_rainy_nonspatial', n.batch = 6000, batch.length = 25, n.thin = 50, n.chains = 4) {
  
  # create directory for outputs if needed
  suppressWarnings(dir.create(outdir, recursive = T))
  
  # center and scale occupancy and detection covariates
  data_list$occ.covs <- apply(data_list$occ.covs, 2, scale)
  
  data_list$det.covs <- lapply(data_list$det.covs, scale)
  
  # data_list$det.covs[1:(length(data_list$det.covs) - 1)] <- lapply(data_list$det.covs[1:(length(data_list$det.covs) - 1)], scale)
  # 
  # data_list$det.covs$precip <- t(apply(data_list$det.covs$precip, 1, scale))

  
    # Priors ------------------------------
    prior.list <- list(beta.comm.normal = list(mean = 0, var = 2.72),
                       alpha.comm.normal = list(mean = 0, var = 2.72),
                       tau.sq.beta.ig = list(a = 0.1, b = 0.1),
                       tau.sq.alpha.ig = list(a = 0.1, b = 0.1)
                       )
    
    # Run the model -----------------------------------------------------------
    n.factors <- 2
    n.samples <- n.batch * batch.length
    n.burn <- 0.6*(n.batch*batch.length)
    
#### model with main effects of cover, distance to water, distance to forest edge and basal area

out_cweb <- lfMsPGOcc(occ.formula = ~ cover + dist_water_500 + dist_edge_500 + basal_area,
                         det.formula = ~ shrubs + basal_area,
                         data = data_list, 
                         priors = prior.list, 
                         n.factors = n.factors, 
                         n.samples = n.samples,
                         n.omp.threads = 1,
                         verbose = TRUE,
                         n.report = 1000,
                         n.burn = n.burn,
                         n.thin = n.thin,
                         n.chains = n.chains, 
                         k.fold = 4,
                         k.fold.threads = 4, 
                         k.fold.seed = 73)
    
save(out_cweb, file = paste(outdir, "/", data_set_name, "_cweb-", n.chains, "-chain-", Sys.Date(), ".R", sep = ''))    
    
### model with main effects of cover, distance to water, and distance to forest edge
out_cwe <- lfMsPGOcc(occ.formula = ~ cover + dist_water_500 + dist_edge_500,
                     det.formula = ~ shrubs + basal_area,
                     data = data_list, 
                     priors = prior.list, 
                     n.factors = n.factors, 
                     n.samples = n.samples,
                     n.omp.threads = 1,
                     verbose = TRUE,
                     n.report = 1000,
                     n.burn = n.burn,
                     n.thin = n.thin,
                     n.chains = n.chains, 
                     k.fold = 4,
                     k.fold.threads = 4, 
                     k.fold.seed = 73)

save(out_cwe, file = paste(outdir, "/", data_set_name, "_cwe-", n.chains, "-chain-", Sys.Date(), ".R", sep = ''))

#### model with main effects of cover, distance to water, and basal area
out_cwb <- lfMsPGOcc(occ.formula = ~ cover + dist_water_500 + basal_area,
                     det.formula = ~ shrubs + basal_area,
                     data = data_list, 
                     priors = prior.list, 
                     n.factors = n.factors, 
                     n.samples = n.samples,
                     n.omp.threads = 1,
                     verbose = TRUE,
                     n.report = 1000,
                     n.burn = n.burn,
                     n.thin = n.thin,
                     n.chains = n.chains, 
                     k.fold = 4,
                     k.fold.threads = 4, 
                     k.fold.seed = 73)

save(out_cwb, file = paste(outdir, "/", data_set_name, "_cwb-", n.chains, "-chain-", Sys.Date(), ".R", sep = ''))

#### model with main effects of distance to water, distance to forest edge and basal area
out_web <- lfMsPGOcc(occ.formula = ~ dist_water_500 + dist_edge_500 + basal_area,
                     det.formula = ~ shrubs + basal_area,
                     data = data_list, 
                     priors = prior.list, 
                     n.factors = n.factors, 
                     n.samples = n.samples,
                     n.omp.threads = 1,
                     verbose = TRUE,
                     n.report = 1000,
                     n.burn = n.burn,
                     n.thin = n.thin,
                     n.chains = n.chains, 
                     k.fold = 4,
                     k.fold.threads = 4, 
                     k.fold.seed = 73)

save(out_web, file = paste(outdir, "/", data_set_name, "_web-", n.chains, "-chain-", Sys.Date(), ".R", sep = ''))

#### model with main effects of cover, distance to forest edge and basal area
out_ceb <- lfMsPGOcc(occ.formula = ~ cover + dist_edge_500 + basal_area,
                     det.formula = ~ shrubs + basal_area,
                     data = data_list, 
                     priors = prior.list, 
                     n.factors = n.factors, 
                     n.samples = n.samples,
                     n.omp.threads = 1,
                     verbose = TRUE,
                     n.report = 1000,
                     n.burn = n.burn,
                     n.thin = n.thin,
                     n.chains = n.chains, 
                     k.fold = 4,
                     k.fold.threads = 4, 
                     k.fold.seed = 73)

save(out_ceb, file = paste(outdir, "/", data_set_name, "_ceb-", n.chains, "-chain-", Sys.Date(), ".R", sep = ''))
   
### model with main effects of cover and distance to water
out_cw <- lfMsPGOcc(occ.formula = ~ cover + dist_water_500,
                    det.formula = ~ shrubs + basal_area,
                    data = data_list, 
                    priors = prior.list, 
                    n.factors = n.factors, 
                    n.samples = n.samples,
                    n.omp.threads = 1,
                    verbose = TRUE,
                    n.report = 1000,
                    n.burn = n.burn,
                    n.thin = n.thin,
                    n.chains = n.chains, 
                    k.fold = 4,
                    k.fold.threads = 4, 
                    k.fold.seed = 73)

save(out_cw, file = paste(outdir, "/", data_set_name, "_cw-", n.chains, "-chain-", Sys.Date(), ".R", sep = ''))
    
### model with main effects of cover and distance to forest edge
out_ce <- lfMsPGOcc(occ.formula = ~ cover + dist_edge_500,
                    det.formula = ~ shrubs + basal_area,
                    data = data_list, 
                    priors = prior.list, 
                    n.factors = n.factors, 
                    n.samples = n.samples,
                    n.omp.threads = 1,
                    verbose = TRUE,
                    n.report = 1000,
                    n.burn = n.burn,
                    n.thin = n.thin,
                    n.chains = n.chains, 
                    k.fold = 4,
                    k.fold.threads = 4, 
                    k.fold.seed = 73)

save(out_ce, file = paste(outdir, "/", data_set_name, "_ce-", n.chains, "-chain-", Sys.Date(), ".R", sep = ''))

#### model with main effects of cover and basal area
out_cb <- lfMsPGOcc(occ.formula = ~ cover + basal_area,
                    det.formula = ~ shrubs + basal_area,
                    data = data_list, 
                    priors = prior.list, 
                    n.factors = n.factors, 
                    n.samples = n.samples,
                    n.omp.threads = 1,
                    verbose = TRUE,
                    n.report = 1000,
                    n.burn = n.burn,
                    n.thin = n.thin,
                    n.chains = n.chains, 
                    k.fold = 4,
                    k.fold.threads = 4, 
                    k.fold.seed = 73)

save(out_cb, file = paste(outdir, "/", data_set_name, "_cb-", n.chains, "-chain-", Sys.Date(), ".R", sep = ''))

### model with main effects of distance to water and distance to forest edge
out_we <- lfMsPGOcc(occ.formula = ~ dist_water_500 + dist_edge_500,
                    det.formula = ~ shrubs + basal_area,
                    data = data_list, 
                    priors = prior.list, 
                    n.factors = n.factors, 
                    n.samples = n.samples,
                    n.omp.threads = 1,
                    verbose = TRUE,
                    n.report = 1000,
                    n.burn = n.burn,
                    n.thin = n.thin,
                    n.chains = n.chains, 
                    k.fold = 4,
                    k.fold.threads = 4, 
                    k.fold.seed = 73)

save(out_we, file = paste(outdir, "/", data_set_name, "_we-", n.chains, "-chain-", Sys.Date(), ".R", sep = ''))
  
#### model with main effects of distance to water and basal area
out_wb <- lfMsPGOcc(occ.formula = ~ dist_water_500 + basal_area,
                    det.formula = ~ shrubs + basal_area,
                    data = data_list, 
                    priors = prior.list, 
                    n.factors = n.factors, 
                    n.samples = n.samples,
                    n.omp.threads = 1,
                    verbose = TRUE,
                    n.report = 1000,
                    n.burn = n.burn,
                    n.thin = n.thin,
                    n.chains = n.chains, 
                    k.fold = 4,
                    k.fold.threads = 4, 
                    k.fold.seed = 73)
  
save(out_wb, file = paste(outdir, "/", data_set_name, "_wb-", n.chains, "-chain-", Sys.Date(), ".R", sep = ''))

#### model with main effects of distance to forest edge and basal area
out_eb <- lfMsPGOcc(occ.formula = ~ dist_edge_500 + basal_area,
                    det.formula = ~ shrubs + basal_area,
                    data = data_list, 
                    priors = prior.list, 
                    n.factors = n.factors, 
                    n.samples = n.samples,
                    n.omp.threads = 1,
                    verbose = TRUE,
                    n.report = 1000,
                    n.burn = n.burn,
                    n.thin = n.thin,
                    n.chains = n.chains, 
                    k.fold = 4,
                    k.fold.threads = 4, 
                    k.fold.seed = 73)

save(out_eb, file = paste(outdir, "/", data_set_name, "_eb-", n.chains, "-chain-", Sys.Date(), ".R", sep = ''))

#### null model (intercept only)
  out_null <- lfMsPGOcc(occ.formula = ~ 1,
                        det.formula = ~ 1,
                        data = data_list, 
                        priors = prior.list, 
                        n.factors = n.factors, 
                        n.samples = n.samples,
                        n.omp.threads = 1,
                        verbose = TRUE,
                        n.report = 1000,
                        n.burn = n.burn,
                        n.thin = n.thin,
                        n.chains = n.chains, 
                        k.fold = 4,
                        k.fold.threads = 4, 
                        k.fold.seed = 73)

  save(out_null, file = paste(outdir, "/", data_set_name, "_null-", n.chains, "-chain-", Sys.Date(), ".R", sep = ''))
}

# Run the models ----------------------------------------------------------

sapply(names(data_bundles), function(x) run_models(data_bundles[[x]], x,  outdir = 'results/spOccupancy_fits/julian_dry_rainy_nonspatial'), simplify = F)


