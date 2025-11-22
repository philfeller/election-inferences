# Utility functions used when creating election inferences
## - calculate_residuals(ei.model, results, beg_yr, end_yr): Calculate residuals from an ei.MD.bayes model
## - regression_model(results, ei.model, beg_yr, end_yr, selection_method): Build stepwise regression models for each party
## - build_ei_model(beg_yr, end_yr, lambda1, lambda2, covariate, results_tibble): Tune and build ei.MD.bayes model for a given pair of years

library(eiPack)
source("./global.R")

# Build stepwise regression models for each party's residuals; return the covariates
regression_model <- function(results, ei.model, beg_yr, end_yr, selection_method = "all") {
  # results: tibble with election results and covariates
  # ei.model: ei.MD.bayes model object
  # beg_yr: beginning year of election transition
  # end_yr: ending year of election transition
  # selection_method: "all" to return all covariates selected in any model;
  #                   "majority" to return only covariates selected in at least
  #                   half of the models

  beg_party <- get(paste("p", substring(as.character(beg_yr), 3, 4), sep = ""))
  end_party <- get(paste("p", substring(as.character(end_yr), 3, 4), sep = ""))
  resid <- calculate_residuals(ei.model, results, beg_yr, end_yr)
  data <- cbind(resid, results)
  base_cols <- paste(colnames(results %>% select(lon, lat, gini:pct_farm_1860)), sep = " + ")
  lag_cols <- c()
  for (party in beg_party) {
    lag_col <- paste("lag_", party, sep = "")
    lag_cols <- c(lag_cols, lag_col)
  }
  cols <- c(base_cols, lag_cols)
  test_cols <- paste(cols, collapse = " + ")

  step_models <- map(end_party, function(y) {
    form <- as.formula(paste(y, test_cols, sep = "~"))
    lm_full <- lm(form, data, weights = data$ELIG_1855)
    step(lm_full, direction = "both", trace = 0)
  })

  selected_covariates <- map(step_models, ~ names(coef(.x))[-1])
  var_counts <- table(unlist(selected_covariates))
  majority_covars <- names(var_counts[var_counts >= length(step_models) / 2])

  if (selection_method == "all") {
    covars <- unique(unlist(selected_covariates))
  } else if (selection_method == "majority") {
    covars <- majority_covars
  } else {
    stop("selection_method must be 'all' or 'majority'")
  }
  return(covars)
}

# Although the model factors for each pair of years have been chosen to produce MCMC results
# that converge fairly reliably, individual runs sometimes fail to converge. The calls to
# the ei.MD.bayes function are wrapped in a loop that tests for convergence to ensure that
# results are useable.

# Define function to tune and build ei.MD.bayes model for a given pair of years
build_ei_model <- function(beg_yr, end_yr, lambda1, lambda2, covariate = FALSE, results_tibble = NULL) {
  # beg_yr: beginning year of election transition
  # end_yr: ending year of election transition
  # lambda1: shape parameter for gamma prior on Dirichlet concentration parameters
  # lambda2: scale parameter for gamma prior on Dirichlet concentration parameters
  # covariate: whether to include covariates in the model
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

  # Run tuning and model building with default hyperparameters; mean and
  # variance of the gamma posteriors will be used to set model factors
  tune.default <- tuneMD(
    formula_eimd,
    data = data, total = total, totaldraws = tune_size, ntunes = 10
  )

  ei.model <- ei.MD.bayes(
    formula_eimd,
    data = data, total = total,
    sample = sample, burnin = burnin, thin = thin,
    tune.list = tune.default
  )

  # Calculate mean and variance of Alpha posteriors from default model;
  # these will be used to set hyperparameters
  mu <- mean(unlist(list(ei.model$draws$Alpha)), na.rm = TRUE)
  sigma2 <- var(unlist(list(ei.model$draws$Alpha)), na.rm = TRUE)
  desired_max_ratio <- 4  # Posterior variance should not be more than 4 times prior variance
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

  cols <- regression_model(data, ei.model, beg_yr, end_yr, "majority")
  if (covariate == TRUE && length(cols) > 0) {
    covariate_string <- paste("~", paste(cols, collapse = " + "))
    print(paste("Using covariates for", beg_yr, "to", end_yr, ":", covariate_string))
    save(covariate_string, file = paste(end_yr, "_", beg_yr, "_covariates.Rda", sep = ""))
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
  return(ei.model)
}
