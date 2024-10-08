# Extract individual simulations of beta values, either statewide or for a town
# The output is particularly useful for plotting beta distributions

library(coda)

beta.sims.MD <- function(ei.model, cols, town.id) {
  # Define functions used to extract information from the MCMC model
  getY <- function(x) x[[1]]$y
  get2 <- function(x) x[2]

# Use Heidelberger and Welsh to determine the point at which all betas converge.
  min_conv <- max(heidel.diag(lambda.MD(ei.model, cols))[, 2])
  num_sims <- length(ei.model$draws$Beta[, 1])

  # The model contains the MCMC-created betas by district for each row/column combination;
  # the lambda.MD function returns overall betas
  if (missing(town.id)) {
    x <- lambda.MD(ei.model, cols)
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
# Ideally, only the confirmed converged simulations should be used—
# heidel.diag(lambda.MD(ei.model, cols))[,2]:1000—but the MCMC output in this
# project shows no substantial change in the estimated distribution, not enough
# to justify the added complexity
betas.MD <- function(beta_sims) {
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

# Create single-row data frame with the results for a given year
get_shares <- function(results, yr, town.id) {
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

get_names <- function(results, yr) {
  vote_cols <- paste("vote_in_", yr, sep = "")
  return(colnames(results %>% select(ends_with(vote_cols))))
}

get_share_names <- function(results, yr) {
  vote_cols <- paste("vote_in_", yr, sep = "")
  in_cols <- paste("in_", yr, sep = "")
  return(colnames(results %>% select(ends_with(in_cols) & !ends_with(vote_cols))))
}

calculate_residuals <- function(ei.model, results, beg_yr, end_yr) {
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

evaluate_residuals <- function(residuals) {
  suppressWarnings({
    for (i in 1:ncol(residuals)) {
      if (ks.test(residuals[, i], "pnorm", 0, sd(residuals[, i]))$p.value < 0.05) {
        print(paste("Residuals for", colnames(residuals)[i], "are not normally distributed."))
      }
    }
  })
}
