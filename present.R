# Functions for presenting the results of the MD inference
# construct_contingency(), boxplotMD(), ridgelineMD(), densityMD(),
# create_map(), corr_matrix_by_party()

source("./global.R")

# Define functions useful for extracting tabular data from betas and for
# presenting the information

# Ensure that rounded percentages total 100
pct_round <- function(x, sig_fig) {
  # x: vector of percentages summing to 1
  # sig_fig: number of significant figures to round to

  factor <- 100 * 10^sig_fig
  x <- x * factor # multiply by a factor that reflects the desired number of significant digits
  y <- floor(x) # truncate to nearest tenth of a percent
  # allocate the total remainder, beginning with the greatest difference
  indices <- tail(order(x - y), round(sum(x)) - sum(y))
  y[indices] <- y[indices] + 1
  y / factor
}

# Ensure that rounded betas add up to rounded row and column marginals
beta_round <- function(matrix, beg_round, end_round, sig_fig) {
  # matrix: matrix of betas to be rounded
  # beg_round: vector of rounded row marginals
  # end_round: vector of rounded column marginals
  # sig_fig: number of significant figures to round to

  factor <- 100 * 10^sig_fig
  # Return TRUE if the beta has an allocatable remainder for both column and row
  remain <- function(floor, x, y) {
    col_remain <- (beg_round - rowSums(floor))[y]
    row_remain <- (end_round - colSums(floor))[x]
    if (col_remain > 0 && row_remain > 0) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  }
  # Calculate the number of remainders yet to be allocated to a beta
  num_remain <- function() {
    return(sum(beg_round - rowSums(floor)) + sum(end_round - colSums(floor)))
  }
  matrix <- matrix * factor
  beg_round <- beg_round * factor
  end_round <- end_round * factor
  floor <- floor(matrix)
  col <- ncol(matrix)
  row <- nrow(matrix)
  # Repeat until all remainders are allocated
  while (num_remain() > 0) {
    for (i in order(matrix - floor, decreasing = TRUE)) {
      d <- arrayInd(i, dim(matrix))
      x <- d[, 2]
      y <- d[, 1]
      # Only if beta has both row and column remainders
      if (remain(floor, x, y)) {
        floor[i] <- floor[i] + 1
      }
      if (num_remain() == 0) {
        return(floor / factor)
      }
    }
  }
}

# Create single-row data frame with the results for a given year
get_shares <- function(results, yr, townID) {
  # results: tibble with election results
  # yr: year for which to extract shares
  # townID: optional; if provided, extract shares for the specified town only (by row index)

  vote_cols <- paste("vote_in_", yr, sep = "")
  elig_col <- paste("G_", yr, sep = "")
  if (missing(townID)) {
    totals <- apply(results %>%
      select(ends_with(vote_cols), ends_with(elig_col)), 2, sum)
  } else {
    totals <- apply(results[townID, ] %>%
      select(ends_with(vote_cols), ends_with(elig_col)), 2, sum)
  }
  elig <- totals[length(totals)]
  votes <- totals[-1 * length(totals)]
  # Rearrange votes to match order of parties in betas
  party_order <- str_replace(gsub("vote_in_", "", names(votes)), "non18", "Abstaining_18")
  betas_parties <- str_replace(get(paste("p", substring(as.character(yr), 3, 4), sep = "")), "_in_", "_")
  votes <- votes[match(betas_parties, party_order)]
  return(votes / elig)
}

# Create a data frame with the contingency table for voter transitions
# between two years, using the results of an MD inference
construct_contingency <- function(results, betas, beg_yr, end_yr, townID, sig_fig) {
  # results: tibble with election results
  # betas: matrix of betas from the MD inference
  # beg_yr: beginning year of the transition
  # end_yr: ending year of the transition
  # townID: optional; if provided, extract shares for the specified town only (by row index)
  # sig_fig: number of significant figures to round to

  if (missing(sig_fig)) {
    sig_fig <- 1
  }
  beg_shares <- as.vector(get_shares(results, beg_yr, townID))
  c_names <- colnames(betas)
  c_name <- paste("share in", beg_yr, sep = " ")
  r_names <- rownames(betas)
  r_name <- paste("share in", end_yr, sep = " ")
  end_shares <- as.vector(get_shares(results, end_yr, townID))
  beta_shares <- beta_round(
    betas * beg_shares, pct_round(beg_shares, sig_fig),
    pct_round(end_shares, sig_fig), sig_fig
  ) * 100
  beg_shares <- pct_round(beg_shares, sig_fig) * 100
  end_shares <- pct_round(end_shares, sig_fig) * 100
  matrix <- cbind(rbind(beta_shares, end_shares), c(beg_shares, 100))
  colnames(matrix) <- c(c_names, c_name)
  rownames(matrix) <- c(r_names, r_name)
  return(matrix)
}

# Transform betas tibble for presenting in a single chart
combine_betas <- function(betas) {
  # betas: matrix of betas from the MD inference

  beta_df <- data.frame(matrix(ncol = 4, nrow = 0))
  colnames(beta_df) <- c("party", "value", "from", "to")
  for (op in 1:dim(betas)[1]) {
    p <- t(betas[op, , ])
    dep_party <- rownames(betas)[op]
    df <- as.data.frame(p) %>%
      pivot_longer(colnames(p), names_to = "party") %>%
      mutate(value = 100 * value) %>%
      mutate(to = party) %>%
      mutate(party = paste(dep_party, "to", party)) %>%
      mutate(from = dep_party)
    beta_df <- rbind(beta_df, df)
  }
  return(beta_df)
}

# Plot the distributions of the estimated betas
boxplotMD <- function(betas, col_yr, row_yr, town, covariate = FALSE) {
  # betas: matrix of betas from the MD inference
  # col_yr: beginning year of the transition
  # row_yr: ending year of the transition
  # town: optional; if provided, include town name in title

  x_label <- ifelse(col_yr > row_yr, "Party composition (percent)", "Vote shift (percent)")
  pre_title <- ""
  if (!missing(town)) {
    pre_title <- paste(town, ": ", sep = "")
  }
  beta_df <- combine_betas(betas)
  title <- paste("Voter transitions between", col_yr, "and", row_yr, ifelse(covariate, "with covariates", "without covariates"))
  print(ggplot(data = beta_df, aes(y = from, x = value, fill = to)) +
    stat_boxplot(geom = "errorbar") +
    geom_boxplot(outlier.shape = NA) +
    scale_fill_manual(values = party_colors) +
    guides(fill = guide_legend(title = paste(row_yr, "vote"))) +
    labs(title = title, x = x_label, y = paste(col_yr, "vote")) +
    theme_classic())
}

# Plot the ridgeline distributions of the estimated betas
ridgelineMD <- function(betas, col_yr, row_yr, town, covariate = FALSE) {
  # betas: matrix of betas from the MD inference
  # col_yr: beginning year of the transition
  # row_yr: ending year of the transition
  # town: optional; if provided, include town name in title

  x_label <- ifelse(col_yr > row_yr, "Party composition (percent)", "Vote shift (percent)")
  pre_title <- ""
  if (!missing(town)) {
    pre_title <- paste(town, ": ", sep = "")
  }
  beta_df <- combine_betas(betas)
  title <- paste("Voter transitions between", col_yr, "and", row_yr, ifelse(covariate, "with covariates", "without covariates"))
  suppressMessages(
    print(ggplot(data = beta_df, aes(y = party, x = value, fill = to)) +
      ggridges::stat_density_ridges(quantile_lines = FALSE, scale = 4, alpha = 1) +
      ggridges::theme_ridges(center_axis_labels = TRUE, font_size = 12) +
      scale_color_manual(values = party_colors, aesthetics = c("colour", "fill")) +
      guides(fill = guide_legend(title = paste(row_yr, "party"))) +
      labs(title = title, x = x_label, y = ""))
  )
}

# Plot the density distributions of the estimated betas
densityMD <- function(betas, col_yr, row_yr, town) {
  # betas: matrix of betas from the MD inference
  # col_yr: beginning year of the transition
  # row_yr: ending year of the transition
  # town: optional; if provided, include town name in title

  pre_title <- ""
  if (!missing(town)) {
    pre_title <- paste(town, ": ", sep = "")
  }
  for (op in 1:dim(betas)[1]) {
    p <- t(betas[op, , ])
    dep_party <- rownames(betas)[op]
    df <- as.data.frame(p) %>%
      pivot_longer(colnames(p), names_to = "party") %>%
      rename(vote_share = value)
    title <- paste(ifelse(dep_party != "Abstaining", "Voted ", ""), dep_party, " in ", col_yr, sep = "")
    print(ggplot(data = df, aes(x = vote_share, y = after_stat(density))) +
      geom_freqpoly(mapping = aes(color = party), binwidth = 0.02) +
      labs(title = title, x = "Vote share", color = paste(row_yr, "party")))
  }
}

# Create a map with filled values
create_map <- function(yr, map, data) {
  # yr: year for which to create the map
  # map: sf object with the map data
  # data: vector of data values to fill the map

  map$data <- data
  ggplot() +
    geom_sf(data = map, aes(fill = data), col = NA) +
    theme_bw() +
    scale_fill_viridis_c()
}

# Create a correlation matrix for statewide office results

corr_matrix_by_party <- function(governor_results, lt_governor_results, secretary_results, treasurer_results, probate_results, party) {
  # Set office order and display names
  office_names <- c("Governor", "Lieutenant Governor", "Secretary of the State", "Treasurer", "Judge of Probate")
  # Get the correct vote column name for the selected party
  vote_col <- paste0(party, "_votes")

  # Prepare a named list so we can loop programmatically, ordered for display
  dfs <- list(
    "Governor" = governor_results,
    "Lieutenant Governor" = lt_governor_results,
    "Secretary of the State" = secretary_results,
    "Treasurer" = treasurer_results,
    "Judge of Probate" = probate_results
  )

  # Build empty correlation matrix with readable row/col order
  correlation_table <- matrix(NA, nrow = 5, ncol = 5)
  rownames(correlation_table) <- office_names
  colnames(correlation_table) <- office_names

  # Compute pairwise correlations for the given party and fill the matrix symmetrically
  for (i in 1:5) {
    for (j in i:5) {
      cor_val <- NA
      # Only compute if columns exist
      if (vote_col %in% colnames(dfs[[office_names[i]]]) && vote_col %in% colnames(dfs[[office_names[j]]])) {
        cor_val <- round(cor(dfs[[office_names[i]]][[vote_col]], dfs[[office_names[j]]][[vote_col]], use = "pairwise.complete.obs"), 3)
      }
      correlation_table[i, j] <- cor_val
      correlation_table[j, i] <- cor_val
    }
  }

  # Turn to dataframe for printing with knitr::kable
  correlation_df <- as.data.frame(correlation_table)
  return(correlation_df)
}
