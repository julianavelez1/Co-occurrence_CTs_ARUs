This repository contains code and data associated with:

Vélez, J., McShea, W., Pukazhenthi, B., Stevenson, P. and J. Fieberg. 2023. Implications of the scale of detection for inferring co-occurrence patterns from paired camera traps and acoustic recorders. *Conservation Biology*.

Here we investigate the association between two measures of disturbance (poaching and livestock) and wild ungulates using data collected using camera traps and autonomous acoustic recording units. We quantified these associations using joint species distribution models (JSDMs) fit to data collected in multifunctional landscapes of the Orinoquía region of Colombia. We also evaluated the effect of the detection scale of different sampling methods (e.g., camera traps and acoustic recorders) for inferring co-occurrence patterns between wildlife and disturbance factors.

## Data

- cams_operation_site.RDS: camera traps operation matrix
- covariates.RDS: data frame of site covariates
- data_bundles.rdata: lists containing detection-nondetection arrays, occupancy and detection covariates and site coordinates.
- indep_recs.RDS: list containing independent records of the surveys with camera traps (1), acoustic recorders (2) and pooled records from camera traps and acoustic recorders (3)
- moth_operation_site.RDS: acoustic recorders operation matrix

## Programs

- 1_detection_array.R: creates data bundles to fit JSDMs. 
Input: data/indep_recs.RDS, data/cams_operation_site.RDS, data/moth_operation_site.RDS, data/covariates.RDS. 
Output: data/data_bundles.rdata
- 2a_sfMsPGOcc.R: fits multi-species spatial occupancy model with species correlations.
Input: data/indep_recs.RDS. 
Output: results/spOccupancy_fits/julian_dry_rainy_spatial
- 2b_lfMsPGOcc.R: fits multi-species non-spatial occupancy model with species correlations. 
Input: data/indep_recs.RDS. 
Output: results/spOccupancy_fits/julian_dry_rainy_nonspatial
- 3_waic_convergence: gets WAIC and assesses model convergence. 
Input: results/spOccupancy_fits. 
Output: results/waic/waic_df.RDS, results/best_mods/fits_best_models.RDS.
- 4_covariates_effects.R: gets and plots occupancy estimates and covariates effects. 
Input: results/best_mods/fits_best_models.RDS. 
Output: results/estimates, figures/
- 5_predictions.R: gets model predictions for the effect of covariates on occupancy probability. 
Input: data/data_bundles.rdata, results/best_mods/fits_best_models.RDS. 
Output: figures/marginal_*.RDS
- 6_predictons.R: gets species correlations
Input: data/data_bundles.rdata, results/best_mods/fits_best_models.RDS. 
Output: Output: results/estimates/ci_corr.RDS, figures/corr_*.RDS

