# 1_detection_array.R: this script creates an array containing species, sites and 
#                      replicates. This array is the detection-nondetection data 
#                      used to fit joint species distribution models using 
#                      spOccupancy 

rm(list=ls()) 

library(purrr)
library(tidyverse)
library(here)
library(Matrix)

source(here("functions", "prep_species_array.R"))

# Read independent records and operation matrices -----------------------------------

# Independent records
indep_recs <- readRDS(here("data", "indep_recs.RDS")) # contains independent records of the camera trap survey (1), acoustic survey with audiomoths (2) and pooled records from camera traps and audiomoths (3)

# Set unique values (1 and NA) in operation matrices. Rows are sites and columns are operational dates

cams_op <- readRDS(here("data", "cams_operation_site.RDS"))
cams_op[cams_op %in% 0] <- NA
cams_op[cams_op %in% 0.5] <- 1

unique(c(cams_op))

moths_op <- readRDS(here("data", "moth_operation_site.RDS"))
moths_op[moths_op %in% 0] <- NA
moths_op[moths_op %in% 2] <- 1

unique(c(moths_op))


# Remove days with no moths operating (i.e., columns with only NAs)
cols_na_count <- map(as_tibble(moths_op), ~sum(is.na(.))) %>% bind_rows() %>% t()
no_na <- names(cols_na_count[cols_na_count!=nrow(moths_op),])
moths_op <- moths_op[,no_na]

# Verify NA counts
cols_na_count <- map(as_tibble(moths_op), ~sum(is.na(.))) %>% bind_rows() %>% t()
no_na <- names(cols_na_count[cols_na_count==nrow(moths_op),])

# Format and select covariates --------------------------------------------------------------

covs_cams_moths <- readRDS(here("data", "covariates.RDS")) 
covs_cams_moths <- covs_cams_moths[order(covs_cams_moths$placename),]

coords_cams_moths <- covs_cams_moths %>% 
  select(-geometry) %>% 
  arrange(placename) %>%
  mutate(x = long/1000000, y = lat/1000000) %>% 
  select(x, y)

occ_covs_cams_moths <- covs_cams_moths %>% 
  select(core_500, dist_water_500, cover, 
         dist_edge_500, basal_area, shrubs)

# Trim operation matrices via joined operation of cameras and moths ---------------------------------------

cams_moths_op  <- cams_op[,colnames(cams_op) %in% colnames(moths_op)] # same operational dates for cameras and moths

cams_moths_op <- cams_moths_op*moths_op # propagate NAs from moths, we need both sources of NAs

cams_moths_op <- cams_moths_op[!rownames(cams_moths_op) %in% c("M00", "N18", "N31"), ] # remove sites where cameras or moths were stolen

# Get sampling effort  -----------------------------------------------------

operation_cams_moths <- rowSums(cams_moths_op, na.rm = TRUE) %>% as.data.frame() 
op_days_cams_moths <- cbind(rownames(operation_cams_moths), operation_cams_moths)
colnames(op_days_cams_moths) <- c("placename", "op_days")
op_days_cams_moths <- op_days_cams_moths %>% arrange(placename)

mean(op_days_cams_moths$op_days)

sum(op_days_cams_moths$op_days)


# Get covariates lists ----------------------------------------------------

det_covs_cams_moths <- list(basal_area = covs_cams_moths$basal_area,
                      shrubs = covs_cams_moths$shrubs,
                      cover = covs_cams_moths$cover)


covs_list_cams_moths <- list(det.covs = det_covs_cams_moths, 
                        occ.covs = occ_covs_cams_moths, 
                        coords = coords_cams_moths)


# Get long version of operation matrix -------------

op_names <- "cams_moths_op"

op_matrices <- list(cams_moths_op)

names(op_matrices) <- op_names

op_matrices_ix <- op_matrices %>% 
  map(~.x %>% 
        Matrix(., sparse = T) %>% 
        summary(.))

op_matrices_dimnames <- op_matrices %>% 
  map(~.x %>% 
        dimnames(.))

# Placename, date, and status of operation

op_matrices_l <- map2(
  .x = op_matrices_dimnames,
  .y = op_matrices_ix,
  .f = ~{
    data.frame(placename = .x[[1]][.y$i],
               date_string = .x[[2]][.y$j], 
               status = .y$x)
  }
)  


# Species codes -----------------------------------------------------------

sp.codes <-  c("Cattle", "Brocket Deer",
               "Collared Peccary", "Lowland Tapir",
               "Poaching", "White-lipped Peccary", 
               "White-tailed Deer")

# CT: join detections and operation matrix ----------

cams_recs <- indep_recs %>% pluck("indep_ct")

# Get daily detections
         
cams_recs <- cams_recs %>% 
  filter(species %in% sp.codes) %>%
  mutate(date_string = format(date_time, '%Y-%m-%d'), tmp = 1)  %>%
  select(placename, species, method, date_string, tmp) %>% 
  group_by(placename, species, date_string) %>% 
  summarise(tmp = sum(tmp)) # every detection gets a one and then we add them by day

# Long format of detection-nondetection data for all sampling days
cams_recs_l <- cams_recs %>%
  spread(date_string, tmp, fill = 0) %>%  
  gather(date_string, det, -placename, -species)

cams_recs_ls <- list(cams_recs_l)

# Get CTs operation matrices
op_matrices_l_ct <- names(op_matrices_l) %>% 
  str_detect('cams') %>% 
  keep(op_matrices_l, .)

# Matrix of placename, species and columns of dates/days with detection-nondetection data

recs_operation_ct <- map2(
  .x = cams_recs_ls,
  .y = op_matrices_l_ct,
  .f = ~{
    right_join(.x, .y) %>% 
      mutate(det = det*status, 
             year_day = format(strptime(date_string, format = '%Y-%m-%d'), 
                                '%Y_%j')) %>% 
      group_by(placename, species, year_day)  %>% 
      summarise(det = sum(det, na.rm = T)) %>% # add daily detections
      spread(year_day, det) %>% # move detections into columns with column names as julian day
      filter(!is.na(species))
  }
) 

names(recs_operation_ct) <- names(op_matrices_l_ct)


# Create array species*site*day with detection-nondetection data (0,1)

y_array_ct <- recs_operation_ct %>%
  map(~.x %>% 
        prep_species_array)

# Get raw detections
dets_ct <- apply(y_array_ct$cams_moths_op, 1, sum)
names(dets_ct) <- sp.codes

# AU: join detections and operation matrix ----------

# Filter poaching and cattle from audio

cams_recs <- indep_recs %>% pluck("indep_ct")
audio_recs <- indep_recs %>% pluck("indep_au") 

cams_recs <- cams_recs %>% filter(!species %in% c("Cattle", "Poaching")) # wildlife in cameras
audio_recs <- audio_recs %>% filter(species %in% c("Cattle", "Poaching")) # disturbance in audio

cams_au_recs <- rbind(cams_recs, audio_recs)

# Get daily detections

cams_au_recs <- cams_au_recs %>% 
  filter(species %in% sp.codes) %>%
  mutate(date_string = format(date_time, '%Y-%m-%d'), tmp = 1)  %>%
  select(placename, species, method, date_string, tmp) %>% 
  group_by(placename, species, date_string) %>% 
  summarise(tmp = sum(tmp)) # every detection gets a one and then we add them by day

# Long format of detection-nondetection data for all sampling days
cams_au_recs_l <- cams_au_recs %>%
  spread(date_string, tmp, fill = 0) %>%  
  gather(date_string, det, -placename, -species)

cams_au_recs_ls <- list(cams_au_recs_l)

# Get moths operation matrices
op_matrices_l_au <- names(op_matrices_l) %>% 
  str_detect('moths') %>% 
  keep(op_matrices_l, .)

# Matrix of placename, species and columns of dates/days with detection-nondetection data

recs_operation_au <- map2(
  .x = cams_au_recs_ls,
  .y = op_matrices_l_au,
  .f = ~{
    right_join(.x, .y) %>% 
      mutate(det = det*status, 
             year_day = format(strptime(date_string, format = '%Y-%m-%d'), 
                               '%Y_%j')) %>% 
      group_by(placename, species, year_day)  %>% 
      summarise(det = sum(det, na.rm = T)) %>% # add daily detections
      spread(year_day, det) %>% # move detections into columns with column names as julian day
      filter(!is.na(species))
  }
) 

names(recs_operation_au) <- names(op_matrices_l_au)


# Create array species*site*day with detection-nondetection data (0,1)

y_array_au <- recs_operation_au %>%
  map(~.x %>% 
        prep_species_array)

# Get raw detections
dets_au <- apply(y_array_au$cams_moths_op, 1, sum)

names(dets_au) <- sp.codes
dets_au <- dets_au[c("Cattle", "Poaching")]

# Join raw detections (Table 1 in manuscript) -----------------------------------------------------

dets_ct_au <-  data.frame(sp.codes = c("Cattle", "Poaching"), dets_au) %>% 
  full_join(., data.frame(sp.codes, dets_ct)) %>% 
  arrange(sp.codes) %>% select(sp.codes, dets_ct, dets_au)

cols <- c("Species", "CTs", "ARUs")

names(dets_ct_au) <- cols


# Save data ---------------------------------------------------------------

data_dry_rainy_ct <- list(y = y_array_ct$cams_moths_op, 
                      det.covs = det_covs_cams_moths, 
                      occ.covs = occ_covs_cams_moths, 
                      coords = coords_cams_moths)

data_dry_rainy_au <- list(y = y_array_au$cams_moths_op, 
                     det.covs = det_covs_cams_moths, 
                     occ.covs = occ_covs_cams_moths, 
                     coords = coords_cams_moths)

save(data_dry_rainy_ct,
     data_dry_rainy_au,
     sp.codes,
     file = here("data", "data_bundles.rdata"))

save(dets_ct_au, op_days_cams_moths, 
     file = here("results", "raw_dets", "raw_daily_dets.rdata"))
