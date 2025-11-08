# Utility functions used when creating election inferences

library(coda)
library(eiPack)
library(purrr, warn.conflicts = FALSE)

# Build stepwise regression models for each party's residuals; return the covariates
regression_model <- function(results, ei.model, beg_yr, end_yr, selection_method = "all") {
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
  
  tune <- tuneMD(
    formula_eimd,
    data = data, total = total, totaldraws = tune_size, ntunes = 10,
    lambda1 = lambda1, lambda2 = lambda2
  )
  
  h <- c(0, 1)
  while (sum(h) != length(h)) {
    ei.model <- ei.MD.bayes(
      formula_eimd,
      data = data, total = total,
      sample = sample, burnin = burnin, thin = thin,
      tune.list = tune, lambda1 = lambda1, lambda2 = lambda2
    )
    # Check for convergence
    h <- heidel.diag(lambda.MD(ei.model, end_party))[, 1]
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
      h <- heidel.diag(lambda.MD(ei.model, end_party))[, 1]
    }
  }
  return(ei.model)
}

