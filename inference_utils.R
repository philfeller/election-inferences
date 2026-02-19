# Utility functions used when creating election inferences
## - calculate_residuals(ei.model, results, beg_yr, end_yr): Calculate residuals from an ei.MD.bayes model
## - regression_model(results, ei.model, beg_yr, end_yr, selection_method): Build stepwise regression models for each party
## - build_ei_model(beg_yr, end_yr, lambda1, lambda2, covariate, results_tibble): Tune and build ei.MD.bayes model for a given pair of years

library(eiPack)
source("./global.R")

library(glmnet)
library(purrr)

# Function to prepare spatial covariates properly
prepare_spatial_data <- function(data, col_names) {
  # Extract columns as matrix
  X <- as.matrix(data[, col_names, drop = FALSE])
  
  # Check for sufficient columns
  if (ncol(X) < 2) {
    stop("Need at least 2 predictor variables for glmnet")
  }
  
  # Remove any columns with all NAs or zero variance
  col_vars <- apply(X, 2, var, na.rm = TRUE)
  valid_cols <- !is.na(col_vars) & col_vars > 0
  X <- X[, valid_cols, drop = FALSE]
  
  if (ncol(X) < 2) {
    stop("After removing constant columns, need at least 2 predictors")
  }
  
  # Standardize to put all variables on same scale
  X_scaled <- scale(X)
  
  list(
    X = X_scaled,
    center = attr(X_scaled, "scaled:center"),
    scale = attr(X_scaled, "scaled:scale"),
    colnames = colnames(X)
  )
}

# Build elastic net regression models for each party's residuals; return the covariates
regression_model <- function(results, ei.model, beg_yr, end_yr, selection_method = "all", seed = 123) {
  # results: tibble with election results and covariates
  # ei.model: ei.MD.bayes model object
  # beg_yr: beginning year of election transition
  # end_yr: ending year of election transition
  # selection_method: "all" to return all covariates selected in any model;
  #                   "majority" to return only covariates selected in at least
  #                   half of the models
  
  set.seed(seed) # For reproducibility
  beg_party <- get(paste("p", substring(as.character(beg_yr), 3, 4), sep = ""))
  end_party <- get(paste("p", substring(as.character(end_yr), 3, 4), sep = ""))
  weight_col <- paste0("ELIG_", end_yr)
  resid <- calculate_residuals(ei.model, results, beg_yr, end_yr)
  data <- cbind(resid, results)
  
  # Build column vector (not formula string)
  base_cols <- colnames(results %>% select(lon:rr_dist, gini:pct_farm_1860))
  lag_cols <- c()
  for (party in beg_party) {
    lag_col <- paste("lag_", party, sep = "")
    lag_cols <- c(lag_cols, lag_col)
  }
  all_cols <- c(base_cols, lag_cols) # Character vector of column names
  
  # Determine covariates using an elastic net with spatial awareness
  
  elastic_models <- map(end_party, function(y) {
    spatial_data <- prepare_spatial_data(data, all_cols)
    X <- spatial_data$X
    y_vec <- data[[y]]
    weights <- data[[weight_col]]
    
    # Remove any NA values
    complete_cases <- complete.cases(X, y_vec, weights)
    X <- X[complete_cases, , drop = FALSE]
    y_vec <- y_vec[complete_cases]
    weights <- weights[complete_cases]
    
    # Stability selection: run elastic net multiple times
    n_stability_runs <- 50
    stability_results <- map(1:n_stability_runs, function(run) {
      # Try multiple alphas
      alpha_grid <- seq(0.1, 0.9, by = 0.1)
      
      cv_results <- map(alpha_grid, function(a) {
        cv_fit <- cv.glmnet(X, y_vec,
                            alpha = a, weights = weights,
                            nfolds = 10, type.measure = "mse"
        )
        list(
          alpha = a,
          lambda = cv_fit$lambda.min,
          mse = min(cv_fit$cvm)
        )
      })
      
      best_alpha_idx <- which.min(map_dbl(cv_results, "mse"))
      best_alpha <- cv_results[[best_alpha_idx]]$alpha
      
      cv_fit <- cv.glmnet(X, y_vec,
                          alpha = best_alpha, weights = weights,
                          nfolds = 10, type.measure = "mse"
      )
      
      final_fit <- glmnet(X, y_vec,
                          alpha = best_alpha, weights = weights,
                          lambda = cv_fit$lambda.min
      )
      
      coef_matrix <- coef(final_fit)
      selected <- rownames(coef_matrix)[which(coef_matrix != 0)]
      selected[selected != "(Intercept)"]
    })
    
    # Count how often each variable is selected
    all_selected <- unlist(stability_results)
    selection_freq <- table(all_selected) / n_stability_runs
    
    # Keep variables selected in at least 60% of runs
    stable_vars <- names(selection_freq)[selection_freq >= 0.6]
    
    # Fit one final model with median alpha to get coefficients
    alpha_grid <- seq(0.1, 0.9, by = 0.1)
    cv_results <- map(alpha_grid, function(a) {
      cv_fit <- cv.glmnet(X, y_vec,
                          alpha = a, weights = weights,
                          nfolds = 10, type.measure = "mse"
      )
      list(alpha = a, mse = min(cv_fit$cvm))
    })
    best_alpha <- cv_results[[which.min(map_dbl(cv_results, "mse"))]]$alpha
    
    cv_fit <- cv.glmnet(X, y_vec,
                        alpha = best_alpha, weights = weights,
                        nfolds = 10, type.measure = "mse"
    )
    final_fit <- glmnet(X, y_vec,
                        alpha = best_alpha, weights = weights,
                        lambda = cv_fit$lambda.min
    )
    
    list(
      model = final_fit,
      cv_model = cv_fit,
      selected_vars = stable_vars,
      selection_frequency = selection_freq,
      alpha = best_alpha,
      lambda = cv_fit$lambda.min,
      scaling = spatial_data
    )
  })
  
  selected_covariates_by_party <- map(elastic_models, "selected_vars")
  selected_covariates <- unique(unlist(selected_covariates_by_party))
  var_counts <- table(unlist(selected_covariates_by_party))
  majority_covars <- names(var_counts[var_counts >= length(end_party) / 2])
  
  if (selection_method == "all") {
    covars <- selected_covariates
  } else if (selection_method == "majority") {
    covars <- majority_covars
  } else {
    stop("selection_method must be 'all' or 'majority'")
  }
  return(covars)
}

estimate_alpha_params <- function(formula_eimd, data, total) {
  # formula_eimd: formula for ei.MD.bayes model
  # data: tibble with election results
  # total: vector with total eligible voters for each unit
  
  # Run tuning and model building with default hyperparameters; mean and
  # variance of the gamma posteriors will be used to set model factors
  tune <- tuneMD(
    formula_eimd,
    data = data, total = total, totaldraws = 10000, ntunes = 10
  )
  
  ei.model <- ei.MD.bayes(
    formula_eimd,
    data = data, total = total,
    sample = 2000, burnin = 1000, thin = 10,
    tune.list = tune
  )
  
  # Calculate mean and variance of Alpha posteriors from default model;
  # these will be used to set hyperparameters
  mu <- mean(unlist(list(ei.model$draws$Alpha)), na.rm = TRUE)
  sigma2 <- var(unlist(list(ei.model$draws$Alpha)), na.rm = TRUE)
  list(mu = mu, sigma2 = sigma2)
}

# Although the model factors for each pair of years have been chosen to produce MCMC results
# that converge fairly reliably, individual runs sometimes fail to converge. The calls to
# the ei.MD.bayes function are wrapped in a loop that tests for convergence to ensure that
# results are useable.

# Define function to tune and build ei.MD.bayes model for a given pair of years
build_ei_model <- function(beg_yr, end_yr, lambda1, lambda2, covariate = FALSE, cov_cols = c(), results_tibble = NULL) {
  # beg_yr: beginning year of election transition
  # end_yr: ending year of election transition
  # lambda1: shape parameter for gamma prior on Dirichlet concentration parameters
  # lambda2: scale parameter for gamma prior on Dirichlet concentration parameters
  # covariate: whether to include covariates in the model
  # cov_cols: optional vector of column names to use as covariates
  
  # results_tibble: optional tibble with election results and covariates; if not
  #                  provided, the function will look for a tibble named
  #                  results.<last two digits of end_yr>
  
  # Define regression formula for parties in the given years
  beg_party <- get(paste("p", substring(as.character(beg_yr), 3, 4), sep = ""))
  end_party <- get(paste("p", substring(as.character(end_yr), 3, 4), sep = ""))
  lhs <- paste(end_party, collapse = ", ")
  rhs <- paste(beg_party, collapse = ", ")
  formula_string <- paste("cbind(", lhs, ") ~ cbind(", rhs, ")")
  formula_eimd <- as.formula(formula_string)
  
  # Define data and total inputs
  if (missing(results_tibble)) {
    data <- get(paste("results.", substring(as.character(end_yr), 3, 4), sep = ""))
  } else {
    data <- results_tibble
  }
  total <- data[[paste("ELIG_", end_yr, sep = "")]]
  
  # Define model factors
  tune_size <- 10000
  sample <- 2000
  thin <- 10
  burnin <- 1000
  
  # Estimate alpha parameters from default model run
  estimated_alpha_params <- estimate_alpha_params(formula_eimd, data, total)
  mu <- estimated_alpha_params$mu
  sigma2 <- estimated_alpha_params$sigma2
  
  desired_max_ratio <- 4 # Posterior variance should not be more than 4 times prior variance
  sigma2_ratios <- data.frame(sigma2 = numeric(), ratio = numeric())
  
  # Re-tune and re-build model until variance of posterior exceeds
  # one-fourth of the prior variance; too restrictive a prior can
  # force the model into a restrictive region, missing the full
  # range of possible values
  while (TRUE) {
    # Calculate hyperparameters
    lambda2 <- mu / sigma2
    lambda1 <- mu * lambda2
    # Tune model with calculated hyperparameters
    tune <- tuneMD(
      formula_eimd,
      data = data, total = total, totaldraws = tune_size, ntunes = 10,
      lambda1 = lambda1, lambda2 = lambda2
    )
    
    # Build model; check for convergence
    h <- c(0, 1)
    while (sum(h) != length(h)) {
      ei.model <- ei.MD.bayes(
        formula_eimd,
        data = data, total = total,
        sample = sample, burnin = burnin, thin = thin,
        tune.list = tune, lambda1 = lambda1, lambda2 = lambda2
      )
      # Check for convergence
      h <- coda::heidel.diag(lambda.MD(ei.model, end_party))[, 1]
    }
    sigma2.post <- var(unlist(list(ei.model$draws$Alpha)), na.rm = TRUE)
    # Create an empty data frame to hold sigma2 values that have been tried and their associated ratios
    # Test whether the posterior variance is at most four times the prior variance
    if (sigma2.post >= desired_max_ratio * sigma2 & nrow(sigma2_ratios) < 4) {
      # Add the current sigma2 and ratio to the data frame
      ratio <- sigma2.post / sigma2
      sigma2_ratios <- rbind(sigma2_ratios, data.frame(sigma2 = sigma2, ratio = ratio))
      sigma2 <- sigma2 * 2
    } else {
      if (nrow(sigma2_ratios) >= 4) {
        # Select the variance with the lowest ratio
        best_row <- sigma2_ratios[which.min(sigma2_ratios$ratio), ]
        sigma2 <- best_row$sigma2
        lambda2 <- mu / sigma2
        lambda1 <- mu * lambda2
        print("Warning: Maximum number of hyperparameter adjustments reached; proceeding with best value.")
      }
      break
    }
  }
  
  if (covariate == TRUE & length(cov_cols) == 0) {
    # Determine covariates using regression model on residuals
    cov_cols <- regression_model(data, ei.model, beg_yr, end_yr, "all")
  }
  
  if (covariate == TRUE && length(cov_cols) > 0) {
    # Save the current model as the base model before adding covariates
    base.ei.model <- ei.model
    # Create covariate formula
    covariate_string <- paste("~", paste(cov_cols, collapse = " + "))
    covariate_formula <- as.formula(covariate_string)
    
    tune_cov <- tuneMD(
      formula_eimd,
      data = data, total = total, totaldraws = tune_size, ntunes = 10,
      lambda1 = lambda1, lambda2 = lambda2, covariate = covariate_formula
    )
    
    h <- c(0, 1)
    while (sum(h) != length(h)) {
      ei.model <- ei.MD.bayes(
        formula_eimd,
        data = data, total = total,
        sample = sample, burnin = burnin, thin = thin,
        tune.list = tune_cov, lambda1 = lambda1, lambda2 = lambda2,
        covariate = covariate_formula
      )
      # Check for convergence
      h <- coda::heidel.diag(lambda.MD(ei.model, end_party))[, 1]
    }
  }
  return(list(model = ei.model, cov = cov_cols))
}

quantify_ei_variability <- function(beg_yr, end_yr, n_runs = 10, covariate = TRUE) {
  results_list <- vector("list", n_runs)
  
  beg_yr_parties <- get(paste0("p", beg_yr - 1800))
  end_yr_parties <- get(paste0("p", end_yr - 1800))
  cov <- c()
  for (run in 1:n_runs) {
    cat("Running inference", run, "of", n_runs, "\n")
    
    ei_results <- build_ei_model(beg_yr, end_yr, covariate = covariate, cov_cols = cov)
    ei_model <- ei_results$model
    cov <- ei_results$cov
    
    # Save draws array for this run
    saveRDS(posterior::as_draws_array(ei_model$draws$Beta),
            file = paste0("ei_draws_run_", run, ".Rds"),
            compress = "xz")
    
    # Extract posterior means and credible intervals for betas
    betas <- beta.sims.MD(ei_model, end_yr_parties)
    
    n_source <- dim(betas)[1]
    n_target <- dim(betas)[2]
    n_iter <- dim(betas)[3]
    
    source_parties <- dimnames(betas)[[1]]
    target_parties <- dimnames(betas)[[2]]
    
    # Calculate summaries
    posterior_mean <- apply(betas, c(1, 2), mean)
    posterior_sd <- apply(betas, c(1, 2), sd)
    posterior_median <- apply(betas, c(1, 2), median)
    posterior_q025 <- apply(betas, c(1, 2), quantile, probs = 0.025)
    posterior_q975 <- apply(betas, c(1, 2), quantile, probs = 0.975)
    
    # Store summaries for this run
    results_list[[run]] <- list(
      n_draws_used = n_iter,
      covariates = cov,
      # covariate_effectiveness = cov_comparison,
      posterior_mean = posterior_mean,
      posterior_sd = posterior_sd,
      posterior_median = posterior_median,
      posterior_q025 = posterior_q025,
      posterior_q975 = posterior_q975
    )
  }
  
  # Calculate confidence interval variability metrics
  variability_analysis <- analyze_inference_variability(results_list)
  
  # Create a list of the draws arrays for all runs
  arr_list <- lapply(1:n_runs, function(i) {
    readRDS(paste0("ei_draws_run_", i, ".Rds"))
  })
  
  # Combine the draws arrays into a single draws object
  combined_draws <- posterior::bind_draws(
    lapply(arr_list, posterior::as_draws_array),
    along = "chain"
  )

  # Calculate Gelmen-Rubin R-hat values for all parameters across all runs
  rhat_values <- posterior::summarise_draws(combined_draws, "rhat")
  
  # Calculate posterior means for key parameters across all runs
  param_means <- lapply(arr_list, function(arr) {
    posterior::summarise_draws(posterior::as_draws_array(arr), "mean")$mean
  })
  
  # Find run closest to median
  param_matrix <- do.call(rbind, param_means)
  median_params <- apply(param_matrix, 2, median)
  
  # Calculate distance from median for each run
  distances <- apply(param_matrix, 1, function(row) {
    sqrt(sum((row - median_params)^2))
  })
  
  closest_run <- which.min(distances)
  median_betas <- arr_list[[closest_run]]
  vote_data <- prepare_vote_data(end_yr)
  
  all_flows <- expand.grid(
    from = beg_yr_parties,
    to = end_yr_parties,
    stringsAsFactors = FALSE
  ) %>%
    rowwise() %>%
    mutate(
      stats = list(flow_intervals(arr_list, from, to, vote_data))
    ) %>%
    unnest_wider(stats)
  
  return(list(
    median_run = median_betas,
    individual_runs = results_list,
    variability = variability_analysis,
    n_successful = length(results_list),
    rhat = rhat_values,
    flow_intervals = all_flows
  ))
}

analyze_inference_variability <- function(results_list) {
  n_runs <- length(results_list)
  
  # Extract posterior means from each run
  means_array <- array(
    dim = c(
      n_runs,
      nrow(results_list[[1]]$posterior_mean),
      ncol(results_list[[1]]$posterior_mean)
    )
  )
  
  for (i in 1:n_runs) {
    means_array[i, , ] <- results_list[[i]]$posterior_mean
  }
  
  # Calculate between-run statistics
  meta_mean <- apply(means_array, c(2, 3), mean)
  meta_sd <- apply(means_array, c(2, 3), sd)
  meta_min <- apply(means_array, c(2, 3), min)
  meta_max <- apply(means_array, c(2, 3), max)
  
  # Calculate coefficient of variation (CV) for each transition
  meta_cv <- meta_sd / meta_mean
  
  # Within-run vs between-run variance decomposition
  within_run_var <- mean(sapply(results_list, function(r) mean(r$posterior_sd^2)))
  between_run_var <- mean(meta_sd^2)
  total_var <- within_run_var + between_run_var
  
  # What proportion of uncertainty is due to model instability vs sampling?
  pct_between_run <- 100 * between_run_var / total_var
  
  # Check if 95% CIs from different runs overlap
  overlap_matrix <- matrix(TRUE, n_runs, n_runs)
  for (i in 1:(n_runs - 1)) {
    for (j in (i + 1):n_runs) {
      # Check if CIs overlap for all transitions
      ci_overlap <- all(
        results_list[[i]]$posterior_q975 >= results_list[[j]]$posterior_q025 &
          results_list[[i]]$posterior_q025 <= results_list[[j]]$posterior_q975
      )
      overlap_matrix[i, j] <- overlap_matrix[j, i] <- ci_overlap
    }
  }
  
  pct_runs_overlap <- 100 * sum(overlap_matrix) / (n_runs * (n_runs - 1))
  
  return(list(
    meta_mean = meta_mean,
    meta_sd = meta_sd,
    meta_cv = meta_cv,
    range = list(min = meta_min, max = meta_max),
    within_run_variance = within_run_var,
    between_run_variance = between_run_var,
    pct_uncertainty_from_instability = pct_between_run,
    pct_ci_overlap = pct_runs_overlap,
    assessment = assess_ci_stability(meta_cv, pct_between_run, pct_runs_overlap)
  ))
}

assess_ci_stability <- function(cv, pct_between, pct_overlap) {
  max_cv <- max(cv, na.rm = TRUE)
  
  if (max_cv < 0.05 && pct_between < 5 && pct_overlap > 95) {
    return("EXCELLENT: Results highly stable across runs")
  } else if (max_cv < 0.10 && pct_between < 10 && pct_overlap > 90) {
    return("GOOD: Results reasonably stable, minor variation")
  } else if (max_cv < 0.20 && pct_between < 20 && pct_overlap > 80) {
    return("FAIR: Noticeable instability, interpret with caution")
  } else {
    return("POOR: High instability, results may not be reliable")
  }
}

# Function to prepare vote data for a given year, for preparing the data that
# will be used to calculate weighted aggregate flows across towns
prepare_vote_data <- function(year) {
  yr_2digit <- sprintf("%02d", year %% 100)
  dataset_name <- paste0("results.", yr_2digit)
  
  if (!exists(dataset_name)) {
    stop(paste("Dataset", dataset_name, "not found"))
  }
  
  results_df <- get(dataset_name)
  prev_year <- year - 1
  
  results_df %>%
    select(town, ends_with(paste0("vote_in_", prev_year))) %>%
    mutate(town_id = row_number()) %>%
    rename_with(
      ~ gsub("_vote_in_", "_in_", .x),
      .cols = contains("_vote_in_")
    ) %>%
    rename_with(
      ~ gsub("nonvote_in_", "Abstaining_in_", .x),
      .cols = starts_with("nonvote_in_")
    ) %>%
    select(town_id, town, ends_with(paste0("_in_", prev_year)))
}

# Calculate weighted aggregate flow across all towns
calculate_aggregate_flow <- function(arr, from, to, vote_data) {
  draws_df <- posterior::as_draws_df(arr)
  
  # Get all town-specific transitions for this flow
  town_params <- grep(paste0("beta\\.", from, "\\.", to, "\\.\\d+"), 
                      colnames(draws_df), value = TRUE)
  
  # Extract town IDs from parameter names
  town_ids <- as.numeric(gsub(paste0(".*\\."), "", town_params))
  
  # Get weights (number of voters in 'from' category for each town)
  weights <- vote_data[[from]][match(town_ids, vote_data$town_id)]
  
  # Remove any NAs and align
  valid_idx <- !is.na(weights) & weights > 0
  weights <- weights[valid_idx]
  town_params <- town_params[valid_idx]
  
  # Normalize weights
  weights <- weights / sum(weights)
  
  # Check if we have valid data
  if (length(weights) == 0 || length(town_params) == 0) {
    return(rep(NA_real_, nrow(draws_df)))
  }
  
  # Weighted average across towns for each posterior draw
  aggregate_samples <- sapply(1:nrow(draws_df), function(i) {
    town_flows <- sapply(town_params, function(p) draws_df[[p]][i])
    
    # Make sure lengths match
    if (length(town_flows) != length(weights)) {
      return(NA_real_)
    }
    
    weighted.mean(town_flows, weights, na.rm = TRUE)
  })
  
  return(aggregate_samples)
}

# Calculate intervals across runs
flow_intervals <- function(arr_list, from, to, vote_data) {
  # Get mean aggregate flow from each run
  flow_estimates <- sapply(arr_list, function(arr) {
    samples <- calculate_aggregate_flow(arr, from, to, vote_data)
    mean(samples, na.rm = TRUE)
  })
  
  # Remove any NAs
  flow_estimates <- flow_estimates[!is.na(flow_estimates)]
  
  if (length(flow_estimates) == 0) {
    return(list(median = NA, mean = NA, min = NA, max = NA, range_pct = NA))
  }
  
  list(
    median = median(flow_estimates),
    mean = mean(flow_estimates),
    min = min(flow_estimates),
    max = max(flow_estimates),
    range_pct = max(flow_estimates) - min(flow_estimates)
  )
}
