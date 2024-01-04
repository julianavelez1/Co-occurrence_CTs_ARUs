# Functions to get species correlations

get.corr <- function(idx, lambda_samp, N.spp, n.factors) {
  matr <- matrix(as.numeric(lambda_samp[idx,1:(N.spp*n.factors)]), ncol = n.factors) # get matrix of each mcmc sample
  lambdalambdaT <- tcrossprod(matr) / (nrow(matr) - 1) # covariance matrix
  cor_out <- cov2cor(lambdalambdaT) # from covariance to correlation
  return(cor_out)
}


get.corr.spOcc <- function(fit = out, N.spp, n.factors) {

  lambda_samp <- fit$lambda.samples %>% as_draws_df() %>% as.data.frame()
  
  corr_array <- sapply(1:nrow(lambda_samp), get.corr, lambda_samp, N.spp, n.factors, simplify = 'array')

  corr_p50 <- apply(corr_array, 1:2, median)
  
  corr_p055 <- apply(corr_array, 1:2, quantile, 0.055)
  
  corr_p945 <- apply(corr_array, 1:2, quantile, 0.945)
  
  corr_p50_sig <- corr_p50 * !(apply(corr_p055, 2, function(x) x < 0) & apply(corr_p945, 2, function(x) x > 0))
  
  return(list(correlation_array = corr_array, median_correlation_sig = corr_p50_sig))
}

