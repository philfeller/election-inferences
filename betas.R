# Extract individual simulations of beta values, either statewide or for a town
# The output is particularly useful for plotting beta distributions

source("./global.R")

# Extract individual simulations of beta values, either statewide or for a town
beta.sims.MD <- function(ei.model, cols, town.id) {
  # ei.model: ei.MD.bayes model object
  # cols: vector of column names corresponding to the marginals for which
  #       beta simulations are desired
  # town.id: optional; if provided, extract beta simulations for the
  #          specified town only (by name)

  # Define functions used to extract information from the MCMC model
  getY <- function(x) x[[1]]$y
  get2 <- function(x) x[2]

  # Use Heidelberger and Welsh to determine the point at which all betas converge.
  min_conv <- max(coda::heidel.diag(eiPack::lambda.MD(ei.model, cols))[, 2])
  num_sims <- length(ei.model$draws$Beta[, 1])

  # The model contains the MCMC-created betas by district for each row/column combination;
  # the lambda.MD function returns overall betas
  if (missing(town.id)) {
    x <- eiPack::lambda.MD(ei.model, cols)
    tnames <- strsplit(colnames(x), "lambda.")
  } else {
    # The model contains the MCMC-created betas by district for each row/column combination;
    # work with the raw betas in order to extract the betas for an individual district
    x <- ei.model$draws$Beta
    tnames <- strsplit(colnames(x), "beta.")
  }
  idx <- strsplit(sapply(tnames, get2), ".", fixed = TRUE)
  idx <- as.list(as.data.frame(matrix(unlist(idx),
    byrow = TRUE,
    nrow = length(idx), ncol =
      length(idx[[1]])
  )))
  idx <- lapply(idx, as.character)
  idx <- lapply(idx, unique)
  if (missing(town.id)) {
    x_dim <- c(sapply(idx, length), nrow(x))
    x_dim_names <- list(
      rows = gsub("_", " ", substring(idx[[1]], 1, nchar(idx[[1]]) - 8)),
      columns = gsub("_", " ", substring(idx[[2]], 1, nchar(idx[[2]]) - 8)),
      simulations = 1:nrow(x)
    )
  } else {
    x_dim <- c(sapply(idx, length), nrow(x))
    x_dim_names <- ldimnames <- list(
      rows = gsub("_", " ", substring(idx[[1]], 1, nchar(idx[[1]]) - 8)),
      columns = gsub("_", " ", substring(idx[[2]], 1, nchar(idx[[2]]) - 8)),
      towns = idx[[3]],
      simulations = 1:nrow(x)
    )
  }
  # Change "Abstaining" to "Non-voting" in dimension names
  x_dim_names <- lapply(x_dim_names, function(n) gsub("Abstaining", "Non-voting", n))
  if (missing(town.id)) {
    x <- array(t(x), dim = x_dim, dimnames = x_dim_names)
  } else {
    x <- array(t(x), dim = x_dim, dimnames = x_dim_names)[, , town.id, ]
  }
  # Drop all simulations before all betas have converged
  return(x[, , min_conv:num_sims])
}

# Combine the output of beta.sims.MD() into a matrix, by the rows and columns
# from which marginals were taken, of the mean beta values.
betas.MD <- function(beta_sims) {
  # beta_sims: output of beta.sims.MD()

  # Function to calculate rounded mean
  rmed <- function(x) round(mean(x), 3)

  r_names <- rownames(beta_sims)
  c_names <- colnames(beta_sims)
  vote_trans <- array(dim = c(length(r_names), length(c_names)), dimnames = list(r_names, c_names))
  # Iterate through the parties in the columns
  for (op in 1:dim(beta_sims)[1]) {
    p <- t(beta_sims[op, , ])
    value <- apply(p, 2, rmed)
    vote_trans[op, ] <- value
  }
  return(vote_trans)
}

# Combine the output of beta.sims.MD() into a matrix, by the rows and columns
# from which marginals were taken, of the 2.5 percentile beta values.
betas_lo.MD <- function(beta_sims) {
  # beta_sims: output of beta.sims.MD()
  
  # Function to calculate rounded mean
  rmed <- function(x) round(quantile(x, 0.025), 3)
  
  r_names <- rownames(beta_sims)
  c_names <- colnames(beta_sims)
  vote_trans <- array(dim = c(length(r_names), length(c_names)), dimnames = list(r_names, c_names))
  # Iterate through the parties in the columns
  for (op in 1:dim(beta_sims)[1]) {
    p <- t(beta_sims[op, , ])
    value <- apply(p, 2, rmed)
    vote_trans[op, ] <- value
  }
  return(vote_trans)
}

# Combine the output of beta.sims.MD() into a matrix, by the rows and columns
# from which marginals were taken, of the .5 percentile beta values.
betas_lo.MD <- function(beta_sims) {
  # beta_sims: output of beta.sims.MD()
  
  # Function to calculate rounded mean
  rmed <- function(x) round(quantile(x, 0.025), 3)
  
  r_names <- rownames(beta_sims)
  c_names <- colnames(beta_sims)
  vote_trans <- array(dim = c(length(r_names), length(c_names)), dimnames = list(r_names, c_names))
  # Iterate through the parties in the columns
  for (op in 1:dim(beta_sims)[1]) {
    p <- t(beta_sims[op, , ])
    value <- apply(p, 2, rmed)
    vote_trans[op, ] <- value
  }
  return(vote_trans)
}

# Create single-row data frame with the results for a given year
get_shares <- function(results, yr, town.id) {
  # results: tibble with election results
  # yr: year for which to extract shares
  # town.id: optional; if provided, extract shares for the specified town only (by row index)

  vote_cols <- paste("vote_in_", yr, sep = "")
  elig_col <- paste("G_", yr, sep = "")
  if (missing(town.id)) {
    totals <- apply(results %>%
      select(ends_with(vote_cols), ends_with(elig_col)), 2, sum)
  } else {
    totals <- apply(results[town.id, ] %>%
      select(ends_with(vote_cols), ends_with(elig_col)), 2, sum)
  }
  elig <- totals[length(totals)]
  votes <- totals[-1 * length(totals)]
  return(votes / elig)
}

# Get the names of the vote columns for a given year
get_names <- function(results, yr) {
  # results: tibble with election results
  # yr: year for which to extract column names

  vote_cols <- paste("vote_in_", yr, sep = "")
  return(colnames(results %>% select(ends_with(vote_cols))))
}

# Get the names of the share columns for a given year
get_share_names <- function(results, yr) {
  # results: tibble with election results
  # yr: year for which to extract share column names

  vote_cols <- paste("vote_in_", yr, sep = "")
  in_cols <- paste("in_", yr, sep = "")
  return(colnames(results %>% select(ends_with(in_cols) & ! starts_with("lag_") & !ends_with(vote_cols))))
}

# Calculate residuals for the model
calculate_residuals <- function(ei.model, results, beg_yr, end_yr) {
  # ei.model: ei.MD.bayes model object
  # results: tibble with election results
  # beg_yr: beginning year of the transition
  # end_yr: ending year of the transition

  col_names <- get_names(results, end_yr)
  share_names <- get_share_names(results, end_yr)
  resid <- data.frame(matrix(ncol = length(col_names), nrow = 0))
  weights <- as.double(unlist(sqrt(results %>%
    select(starts_with("ELIG_") & ends_with(as.character(end_yr))))))
  for (i in 1:nrow(results)) {
    betas <- betas.MD(beta.sims.MD(ei.model, share_names))
    row_marginals <- get_shares(results, beg_yr, i)
    col_marginals <- as.vector(get_shares(results, end_yr, i))
    predicted <- as.vector(apply(betas * row_marginals, 2, sum))
    resid <- rbind(resid, weights[i] * (predicted - col_marginals))
  }
  colnames(resid) <- col_names
  return(resid)
}

# Use Shapiro-Wilk test to evaluate normality of residuals
evaluate_residuals <- function(residuals) {
  # residuals: data frame of residuals, as output by calculate_residuals()

  suppressWarnings({
    for (i in 1:ncol(residuals)) {
      # Shapiro-Wilk test requires between 3 and 5000 samples
      if (nrow(residuals) >= 3 && nrow(residuals) <= 5000) {
        sw <- shapiro.test(residuals[, i])
        if (sw$p.value < 0.05) {
          print(paste("Residuals for", colnames(residuals)[i], "are not normally distributed."))
        }
      } else {
        print(paste("Not enough (or too many) observations for Shapiro-Wilk test on", colnames(residuals)[i]))
      }
    }
  })
}
