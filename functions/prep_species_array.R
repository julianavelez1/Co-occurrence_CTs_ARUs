# Function to prepare species array

prep_species_array <- function(ops_df) {
  # # prepare empty data frame to explicitly represent sites with no detections for a given species
  # sum(x) - sum(x) to propagate NAs
  empty <- ops_df %>% group_by(placename) %>% 
    summarise(across(where(is.numeric), function(x) sum(x) - sum(x)))
  
  N <- length(unique(ops_df$species))
  K <- ncol(empty) - 1
  J <- length(unique(ops_df$placename))
  
  y <- array(NA, dim = c(N, J, K))
  
  for (i in 1:N) {
    tmp <- ops_df  %>% filter(species == sp.codes[i]) %>% 
      ungroup() %>% select(placename, where(is.numeric))
    tmp <- tmp %>% group_by(placename) %>% summarise_all(sum)
    tmp <- rbind(tmp, empty %>% 
                   filter(!placename %in% tmp$placename)) %>% 
      arrange(placename) %>% ungroup() %>% 
      select(where(is.numeric)) %>% 
      replace(is.na(.), 0)
    y[i, , ] <- as.matrix(tmp[, colnames(tmp)])
    y[i, , ] <- ifelse(y[i, , ] > 0, 1, 0)
  }
  return(y)
}