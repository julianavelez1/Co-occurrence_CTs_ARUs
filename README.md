## Co-occurrence_CTs_ARUs

This repository contains data, code and results associated with:

Vélez, J., McShea, W., Pukazhenthi, B., Stevenson, P. and J. Fieberg. 2023. Implications of the scale of detection for inferring co-occurrence patterns from paired camera traps and acoustic recorders. *Conservation Biology*.

The objective of this study was to investigate the association between two measures of disturbance (poaching and livestock) and wild ungulates using data collected using camera traps and autonomous acoustic recording units. We quantified these associations using joint species distribution models (JSDMs) fit to data collected in multifunctional landscapes of the Orinoquía region of Colombia. We also evaluated the effect of the detection scale of camera traps and acoustic recorders for inferring co-occurrence patterns between wildlife and disturbance factors.

## Data

- cams_operation_site.RDS: camera traps operation matrix
- covariates.RDS: data frame of site covariates
- data_bundles.rdata: lists containing detection-nondetection arrays, occupancy and detection covariates and site coordinates.
- indep_recs.RDS: list containing independent records of the surveys with camera traps and acoustic recorders.
- moth_operation_site.RDS: acoustic recorders operation matrix

## Programs

- 1_detection_array.R: creates data bundles to fit JSDMs. 
  - Input: data/indep_recs.RDS, data/cams_operation_site.RDS, data/moth_operation_site.RDS, data/covariates.RDS
  - Output: data/data_bundles.rdata
- 2a_sfMsPGOcc.R: fits multi-species spatial occupancy model with species correlations.
  - Input: data/indep_recs.RDS
  - Output: results/spOccupancy_fits/julian_dry_rainy_spatial
- 2b_lfMsPGOcc.R: fits multi-species non-spatial occupancy model with species correlations. 
  - Input: data/indep_recs.RDS
  - Output: results/spOccupancy_fits/julian_dry_rainy_nonspatial
- 3_waic_convergence.R: gets WAIC and assesses model convergence. 
  - Input: results/spOccupancy_fits
  - Output: results/waic/waic_df.RDS, results/best_mods/fits_best_models.RDS
- 4_covariates_effects.R: gets and plots occupancy estimates and covariates effects. 
  - Input: results/best_mods/fits_best_models.RDS
  - Output: results/estimates, figures/occ_\*.RDS, figures/det_\*.RDS
- 5_predictions.R: gets model predictions for the effect of covariates on occupancy probability. 
  - Input: data/data_bundles.rdata, results/best_mods/fits_best_models.RDS. 
  - Output: figures/marginal_\*.RDS
- 6_sp_correlation.R: gets species correlations
  - Input: data/data_bundles.rdata, results/best_mods/fits_best_models.RDS
  - Output: results/estimates/ci_corr.RDS, figures/corr_\*.RDS

## Session information

R version 4.2.1 (2022-06-23)
Platform: aarch64-apple-darwin20 (64-bit)
Running under: macOS Monterey 12.2.1

Matrix products: default
LAPACK: /Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/lib/libRlapack.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] ggdist_3.2.1      ggcorrplot_0.1.4  tidybayes_3.0.3   posterior_1.3.1   sf_1.0-9         
 [6] coda_0.19-4       spOccupancy_0.5.2 Matrix_1.5-3      here_1.0.1        forcats_1.0.0    
[11] stringr_1.5.0     dplyr_1.1.0       readr_2.1.4       tidyr_1.3.0       tibble_3.1.8     
[16] ggplot2_3.4.1     tidyverse_1.3.2   purrr_1.0.1      

loaded via a namespace (and not attached):
 [1] nlme_3.1-162          fs_1.6.1              lubridate_1.9.2       doParallel_1.0.17    
 [5] webshot_0.5.4         httr_1.4.4            rprojroot_2.0.3       tensorA_0.36.2       
 [9] tools_4.2.1           backports_1.4.1       utf8_1.2.3            R6_2.5.1             
[13] KernSmooth_2.23-20    DBI_1.1.3             colorspace_2.1-0      withr_2.5.0          
[17] tidyselect_1.2.0      compiler_4.2.1        cli_3.6.0             rvest_1.0.3          
[21] arrayhelpers_1.1-0    xml2_1.3.3            checkmate_2.1.0       scales_1.2.1         
[25] classInt_0.4-8        proxy_0.4-27          systemfonts_1.0.4     digest_0.6.31        
[29] minqa_1.2.5           rmarkdown_2.20        svglite_2.1.1         pkgconfig_2.0.3      
[33] htmltools_0.5.4       lme4_1.1-31           dbplyr_2.3.0          fastmap_1.1.0        
[37] rlang_1.0.6           readxl_1.4.2          rstudioapi_0.14       svUnit_1.0.6         
[41] farver_2.1.1          generics_0.1.3        jsonlite_1.8.4        distributional_0.3.1 
[45] googlesheets4_1.0.1   magrittr_2.0.3        kableExtra_1.3.4.9000 Rcpp_1.0.10          
[49] munsell_0.5.0         fansi_1.0.4           abind_1.4-5           lifecycle_1.0.3      
[53] stringi_1.7.12        yaml_2.3.7            MASS_7.3-58.2         grid_4.2.1           
[57] parallel_4.2.1        crayon_1.5.2          lattice_0.20-45       haven_2.5.1          
[61] splines_4.2.1         hms_1.1.2             knitr_1.42            pillar_1.8.1         
[65] boot_1.3-28.1         codetools_0.2-19      reprex_2.0.2          glue_1.6.2           
[69] evaluate_0.20         modelr_0.1.10         vctrs_0.5.2           nloptr_2.0.3         
[73] tzdb_0.3.0            foreach_1.5.2         cellranger_1.1.0      gtable_0.3.1         
[77] RANN_2.6.1            assertthat_0.2.1      xfun_0.40             broom_1.0.3          
[81] e1071_1.7-13          class_7.3-21          googledrive_2.0.0     viridisLite_0.4.1    
[85] gargle_1.3.0          iterators_1.0.14      units_0.8-1           timechange_0.2.0     
[89] ellipsis_0.3.2       