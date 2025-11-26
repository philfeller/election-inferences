# Analyze the effect of covariates on posterior standard deviations of transition probabilities
model_sds <- function(ei_model_with_covariates, ei_model_without_covariates, cell_name_pattern) {
  no_cov_sds <- c()
  cov_sds <- c()
  for (r in rownames(beta.sims.MD(ei_model_without_covariates, cell_name_pattern))) {
    for (c in colnames(beta.sims.MD(ei_model_without_covariates, cell_name_pattern))) {
      no_cov <- as.vector(beta.sims.MD(ei_model_without_covariates, cell_name_pattern)[r, c, ])
      cov <- as.vector(beta.sims.MD(ei_model_with_covariates, cell_name_pattern)[r, c, ])
      no_cov_sd <- sqrt(var(no_cov))
      cov_sd <- sqrt(var(cov))
      no_cov_sds <- c(no_cov_sds, no_cov_sd)
      cov_sds <- c(cov_sds, cov_sd)
      if (no_cov_sd > .1) {
        print(paste(r, "to", c, ":", as.character((no_cov_sd - cov_sd)/ no_cov_sd), "covariate improvement"))
      }
    }
  }
  list(no_cov_sds = no_cov_sds, cov_sds = cov_sds)
}

# Sumary of standard deviations with and without covariates
model_summary <- function(sd_list) {
  no_cov_sds <- sd_list$no_cov_sds
  cov_sds <- sd_list$cov_sds
  no_cov_quantiles <- quantile(no_cov_sds, probs = c(.05, .5, .95))
  cov_quantiles <- quantile(cov_sds, probs = c(.05, .5, .95))
  no_cov_median <- median(no_cov_sds)
  cov_median <- median(cov_sds)
  no_cov_mean <- no_cov_quantiles[2]
  cov_mean <- cov_quantiles[2]
  no_cov_5th <- no_cov_quantiles[1]
  cov_5th <- cov_quantiles[1]
  no_cov_95th <- no_cov_quantiles[3]
  cov_95th <- cov_quantiles[3]
  # Create a summary data frame
  summary_df <- data.frame(
    Statistic = c("5th Percentile", "Median", "Mean", "95th Percentile"),
    No_Covariates = c(no_cov_quantiles[1], no_cov_quantiles[2], no_cov_quantiles[2], no_cov_quantiles[3]),
    With_Covariates = c(cov_quantiles[1], cov_quantiles[2], cov_quantiles[2], cov_quantiles[3])
  )
  return(summary_df)
}
