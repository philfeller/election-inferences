# Utility functions used when creating election inferences
## - calculate_residuals(ei.model, results, beg_yr, end_yr): Calculate residuals from an ei.MD.bayes model
## - regression_model(results, ei.model, beg_yr, end_yr, selection_method): Build stepwise regression models for each party
## - build_ei_model(beg_yr, end_yr, lambda1, lambda2, covariate, results_tibble): Tune and build ei.MD.bayes model for a given pair of years

library(eiPack)
source("./global.R")
source("./model_evaluation.R", local = TRUE)

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

    # Compare model standard deviations with and without covariates
    if (covariate == TRUE && exists("base.ei.model")) {
      sd_comparison <- model_sds(ei.model, base.ei.model, end_party)
      cov_comparison <- model_summary(sd_comparison)
    } else {
      cov_comparison <- NULL
    }
  }

  return(list(model = ei.model, cov = cov_cols, cov_comparison = cov_comparison))
}

quantify_ei_variability <- function(beg_yr, end_yr, n_runs = 10) {
  results_list <- vector("list", n_runs)

  cov <- c()
  for (run in 1:n_runs) {
    cat("Running inference", run, "of", n_runs, "\n")

    end_yr_parties <- get(paste0("p", end_yr - 1800))
    ei_results <- build_ei_model(beg_yr, end_yr, covariate = TRUE, cov_cols = cov)
    ei_model <- ei_results$model
    cov <- ei_results$cov
    cov_comparison <- ei_results$cov_comparison

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
      covariate_effectiveness = cov_comparison,
      posterior_mean = posterior_mean,
      posterior_sd = posterior_sd,
      posterior_median = posterior_median,
      posterior_q025 = posterior_q025,
      posterior_q975 = posterior_q975
    )
  }

  # Calculate variability metrics
  variability_analysis <- analyze_inference_variability(results_list)

  return(list(
    individual_runs = results_list,
    variability = variability_analysis,
    n_successful = length(results_list)
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
    assessment = assess_stability(meta_cv, pct_between_run, pct_runs_overlap)
  ))
}

assess_stability <- function(cv, pct_between, pct_overlap) {
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

select_representative_run <- function(var_results, weight_covariates = FALSE) {
  n_runs <- length(var_results$individual_runs)
  
  # Calculate distance of each run from the meta_mean
  distances <- sapply(1:n_runs, function(i) {
    run_mean <- var_results$individual_runs[[i]]$posterior_mean
    meta_mean <- var_results$variability$meta_mean
    
    # Euclidean distance across all transitions
    sqrt(sum((run_mean - meta_mean)^2))
  })
  
  if (weight_covariates) {
    # Also consider covariate effectiveness
    effectiveness <- sapply(1:n_runs, function(i) {
      cov_eff <- var_results$individual_runs[[i]]$covariate_effectiveness
      if (is.list(cov_eff)) {
        mean(unlist(cov_eff), na.rm = TRUE)
      } else {
        mean(cov_eff, na.rm = TRUE)
      }
    })
    
    norm_dist <- (distances - min(distances)) / (max(distances) - min(distances))
    norm_eff <- (effectiveness - min(effectiveness)) / (max(effectiveness) - min(effectiveness))
    
    combined_score <- norm_dist - 0.3 * norm_eff
    
    median_run_idx <- which.min(combined_score)
    
    cat("Selected run", median_run_idx, "balancing centrality and covariate effectiveness\n")
  } else {
    median_run_idx <- which.min(distances)
    cat("Selected run", median_run_idx, "as most representative (distance only)\n")
  }
  
  cat("Distance from meta mean:", round(distances[median_run_idx], 4), "\n")
  cat("Seed used:", var_results$individual_runs[[median_run_idx]]$seed, "\n")  # This might not exist
  
  # Return just the matrices, not the whole run object
  return(list(
    run_index = median_run_idx,
    results = var_results$individual_runs[[median_run_idx]]$posterior_mean,      # JUST THE MATRIX
    covariates = var_results$individual_runs[[median_run_idx]]$covariates,
    covariate_effectiveness = var_results$individual_runs[[median_run_idx]]$covariate_effectiveness,
    ci_lower = var_results$individual_runs[[median_run_idx]]$posterior_q025,     # JUST THE MATRIX
    ci_upper = var_results$individual_runs[[median_run_idx]]$posterior_q975,     # JUST THE MATRIX
    distance = distances[median_run_idx],
    all_distances = distances
  ))
}

# Analyze covariates across runs
analyze_covariates_across_runs <- function(var_results, run_indices = NULL) {
  if (is.null(run_indices)) {
    run_indices <- 1:length(var_results$individual_runs)
  }

  n_runs <- length(run_indices)

  # Extract all covariates used
  all_covariates <- unique(unlist(lapply(run_indices, function(i) {
    var_results$individual_runs[[i]]$covariates
  })))

  # Count frequency of each covariate
  covariate_freq <- sapply(all_covariates, function(cov) {
    sum(sapply(run_indices, function(i) {
      cov %in% var_results$individual_runs[[i]]$covariates
    }))
  })

  covariate_freq_pct <- 100 * covariate_freq / n_runs

  # Categorize covariates
  core_covariates <- names(covariate_freq)[covariate_freq == n_runs]
  common_covariates <- names(covariate_freq)[covariate_freq >= n_runs * 0.5 &
    covariate_freq < n_runs]
  occasional_covariates <- names(covariate_freq)[covariate_freq < n_runs * 0.5]

  # Average effectiveness when covariate is used
  avg_effectiveness <- sapply(run_indices, function(i) {
    eff <- var_results$individual_runs[[i]]$covariate_effectiveness
    if (is.list(eff)) mean(unlist(eff), na.rm = TRUE) else mean(eff, na.rm = TRUE)
  })

  return(list(
    all_covariates = all_covariates,
    frequency = covariate_freq,
    frequency_pct = covariate_freq_pct,
    core = core_covariates,
    common = common_covariates,
    occasional = occasional_covariates,
    recommended_set = c(core_covariates, common_covariates)
  ))
}

create_consensus_estimate <- function(var_results, consistency_threshold = 0.75) {
  n_runs <- length(var_results$individual_runs)

  # Calculate pairwise similarities between runs
  similarity_matrix <- matrix(0, n_runs, n_runs)

  for (i in 1:n_runs) {
    for (j in 1:n_runs) {
      if (i == j) {
        similarity_matrix[i, j] <- 1
        next
      }

      mean_i <- var_results$individual_runs[[i]]$posterior_mean
      mean_j <- var_results$individual_runs[[j]]$posterior_mean

      # Calculate similarity (inverse of normalized distance)
      distance <- sqrt(sum((mean_i - mean_j)^2))
      max_distance <- sqrt(sum(pmax(mean_i, mean_j)^2)) # Maximum possible distance
      similarity_matrix[i, j] <- 1 - (distance / max_distance)
    }
  }

  # For each run, calculate its average similarity to all other runs
  avg_similarity <- rowMeans(similarity_matrix)

  # Select runs above the consistency threshold
  consistent_runs <- which(avg_similarity >= consistency_threshold)

  if (length(consistent_runs) < 3) {
    warning("Too few consistent runs, lowering threshold")
    consistent_runs <- order(avg_similarity, decreasing = TRUE)[1:min(5, n_runs)]
  }

  cat("Using", length(consistent_runs), "of", n_runs, "runs for consensus\n")
  cat("Selected runs:", consistent_runs, "\n")
  cat(
    "Average similarity of selected runs:",
    round(mean(avg_similarity[consistent_runs]), 3), "\n\n"
  )

  # Calculate weighted average (weight by similarity)
  weights <- avg_similarity[consistent_runs]
  weights <- weights / sum(weights)

  n_source <- nrow(var_results$individual_runs[[1]]$posterior_mean)
  n_target <- ncol(var_results$individual_runs[[1]]$posterior_mean)

  consensus_mean <- matrix(0, n_source, n_target)
  rownames(consensus_mean) <- rownames(var_results$individual_runs[[1]]$posterior_mean)
  colnames(consensus_mean) <- colnames(var_results$individual_runs[[1]]$posterior_mean)

  for (i in seq_along(consistent_runs)) {
    run_idx <- consistent_runs[i]
    consensus_mean <- consensus_mean +
      weights[i] * var_results$individual_runs[[run_idx]]$posterior_mean
  }

  # Calculate consensus credible intervals from selected runs
  consensus_q025 <- matrix(0, n_source, n_target,
    dimnames = list(
      rownames(consensus_mean),
      colnames(consensus_mean)
    )
  )
  consensus_q975 <- consensus_q025

  for (i in seq_along(consistent_runs)) {
    run_idx <- consistent_runs[i]
    consensus_q025 <- consensus_q025 +
      weights[i] * var_results$individual_runs[[run_idx]]$posterior_q025
    consensus_q975 <- consensus_q975 +
      weights[i] * var_results$individual_runs[[run_idx]]$posterior_q975
  }

  # Analyze covariates across all runs
  covariate_analysis <- analyze_covariates_across_runs(var_results)

  return(list(
    consensus_mean = consensus_mean,
    consensus_q025 = consensus_q025,
    consensus_q975 = consensus_q975,
    runs_used = consistent_runs,
    weights = weights,
    covariate_analysis = covariate_analysis,
    avg_similarity = avg_similarity
  ))
}

create_median_estimate <- function(var_results) {
  n_runs <- length(var_results$individual_runs)
  n_source <- nrow(var_results$individual_runs[[1]]$posterior_mean)
  n_target <- ncol(var_results$individual_runs[[1]]$posterior_mean)

  # Create array to hold all runs
  all_means <- array(dim = c(n_runs, n_source, n_target))

  for (i in 1:n_runs) {
    all_means[i, , ] <- var_results$individual_runs[[i]]$posterior_mean
  }

  # Take element-wise median
  median_estimate <- apply(all_means, c(2, 3), median)
  q025_estimate <- apply(all_means, c(2, 3), quantile, probs = 0.25) # IQR instead of CI
  q975_estimate <- apply(all_means, c(2, 3), quantile, probs = 0.75)

  rownames(median_estimate) <- rownames(var_results$individual_runs[[1]]$posterior_mean)
  colnames(median_estimate) <- colnames(var_results$individual_runs[[1]]$posterior_mean)

  # Analyze covariates across all runs
  covariate_analysis <- analyze_covariates_across_runs(var_results)

  return(list(
    median_estimate = median_estimate,
    q25 = q025_estimate,
    q75 = q975_estimate,
    covariate_analysis = covariate_analysis,
    interpretation = "Element-wise median across all runs"
  ))
}

# Choose representative run, consensus estimate, or median estimate,
# Depending on degree of maximum CV

choose_final_estimate <- function(var_results) {
  max_cv <- max(var_results$variability$meta_cv)

  if (max_cv < 0.10) {
    selected <- select_representative_run(var_results, weight_covariates = TRUE)

    chosen_estimate <- list(
      estimates = selected$results,
      ci_lower = selected$ci_lower,
      ci_upper = selected$ci_upper,
      covariates = selected$covariates,
      approach = "representative run",
      max_cv = max_cv
    )
  } else if (max_cv < 0.20) {
    consensus <- create_consensus_estimate(var_results, consistency_threshold = 0.75)

    chosen_estimate <- list(
      estimates = consensus$consensus_mean,
      ci_lower = consensus$consensus_q025,
      ci_upper = consensus$consensus_q975,
      covariates = consensus$covariate_analysis$recommended_set,
      approach = "weighted consensus",
      n_runs_used = length(consensus$runs_used),
      max_cv = max_cv
    )
  } else {
    median_est <- create_median_estimate(var_results)

    chosen_estimate <- list(
      estimates = median_est$median_estimate,
      ci_lower = median_est$q25,
      ci_upper = median_est$q75,
      covariates = median_est$covariate_analysis$recommended_set,
      approach = "element-wise median",
      max_cv = max_cv
    )
  }
  return(chosen_estimate)
}
