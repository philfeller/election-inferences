# Functions for presenting the results of the MD inference
# construct_contingency(), boxplotMD(), ridgelineMD(), densityMD(),
# create_map(), corr_matrix_by_party(), create_rhatio_heatmap(),
# create_transition_uncertainty_heatmap(), create_rhat_town_map()

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

create_uncertainty_decomposition_heatmap <- function(var_results) {
  # Extract party names from row/column names
  source_parties <- rownames(var_results$individual_runs[[1]]$posterior_mean)
  target_parties <- colnames(var_results$individual_runs[[1]]$posterior_mean)
  
  n_runs <- length(var_results$individual_runs)
  n_source <- length(source_parties)
  n_target <- length(target_parties)
  
  # Calculate within-run and between-run variance for each transition
  within_run_var <- matrix(0, n_source, n_target)
  between_run_var <- matrix(0, n_source, n_target)
  pct_between <- matrix(0, n_source, n_target)
  
  rownames(within_run_var) <- rownames(between_run_var) <- rownames(pct_between) <- source_parties
  colnames(within_run_var) <- colnames(between_run_var) <- colnames(pct_between) <- target_parties
  
  for (i in 1:n_source) {
    for (j in 1:n_target) {
      within_var <- mean(sapply(var_results$individual_runs, function(run) {
        run$posterior_sd[i, j]^2
      }))
      
      means <- sapply(var_results$individual_runs, function(run) {
        run$posterior_mean[i, j]
      })
      between_var <- var(means)
      
      within_run_var[i, j] <- within_var
      between_run_var[i, j] <- between_var
      
      total_var <- within_var + between_var
      pct_between[i, j] <- if (total_var > 0) {
        100 * between_var / total_var
      } else {
        0
      }
    }
  }
  
  # Create plot data
  plot_df <- as.data.frame(pct_between) %>%
    rownames_to_column("Source") %>%
    pivot_longer(-Source, names_to = "Target", values_to = "PctBetween") %>%
    mutate(
      Source = factor(Source, levels = source_parties),
      Target = factor(Target, levels = target_parties)
    )
  
  # Create heatmap with granular scale
  ggplot(plot_df, aes(x = Target, y = Source, fill = PctBetween)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.0f%%", PctBetween)), 
              size = 3, color = "white", fontface = "bold") +
    scale_fill_gradientn(
      colors = c(
        "#1a9850",  # Dark green: 0-20% (excellent)
        "#91cf60",  # Light green: 20-30% (good)
        "#d9ef8b",  # Yellow-green: 30-40% (acceptable)
        "#ffffbf",  # Yellow: 40-50% (borderline)
        "#fee08b",  # Light orange: 50-60% (concerning)
        "#fc8d59",  # Orange: 60-70% (problematic)
        "#d73027",  # Red: 70-80% (severe)
        "#a50026",  # Dark red: 80-90% (extreme)
        "#67001f"   # Very dark red: 90-100% (critical)
      ),
      values = scales::rescale(c(0, 20, 30, 40, 50, 60, 70, 80, 90, 100)),
      name = "Between-Run\nUncertainty\n(%)",
      limits = c(0, 100),
      breaks = c(0, 25, 50, 75, 100),
      labels = c("0%", "25%", "50%", "75%", "100%")
    ) +
    labs(
      title = "Model Specification Uncertainty by Transition",
      subtitle = "Percentage of total variance due to model specification (vs. MCMC sampling)",
      x = "Target Party",
      y = "Source Party",
      caption = "Higher values indicate greater sensitivity to model specification"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(size = 14, face = "bold"),
      panel.grid = element_blank()
    )
}

# Updated overlap heatmap (with your inverted scale)
create_overlap_heatmap <- function(var_results) {
  source_parties <- rownames(var_results$individual_runs[[1]]$posterior_mean)
  target_parties <- colnames(var_results$individual_runs[[1]]$posterior_mean)
  
  n_runs <- length(var_results$individual_runs)
  
  overlap_pct <- matrix(0, length(source_parties), length(target_parties))
  rownames(overlap_pct) <- source_parties
  colnames(overlap_pct) <- target_parties
  
  for (i in seq_along(source_parties)) {
    for (j in seq_along(target_parties)) {
      overlaps <- 0
      comparisons <- 0
      
      for (r1 in 1:(n_runs-1)) {
        for (r2 in (r1+1):n_runs) {
          ci1_lower <- var_results$individual_runs[[r1]]$posterior_q025[i, j]
          ci1_upper <- var_results$individual_runs[[r1]]$posterior_q975[i, j]
          ci2_lower <- var_results$individual_runs[[r2]]$posterior_q025[i, j]
          ci2_upper <- var_results$individual_runs[[r2]]$posterior_q975[i, j]
          
          if ((ci1_upper >= ci2_lower) && (ci1_lower <= ci2_upper)) {
            overlaps <- overlaps + 1
          }
          comparisons <- comparisons + 1
        }
      }
      
      # Inverted: 100 - overlap% to show disagreement
      overlap_pct[i, j] <- 100 - 100 * overlaps / comparisons
    }
  }
  
  # Create plot
  overlap_df <- as.data.frame(overlap_pct) %>%
    rownames_to_column("Source") %>%
    pivot_longer(-Source, names_to = "Target", values_to = "Disagreement") %>%
    mutate(
      Source = factor(Source, levels = source_parties),
      Target = factor(Target, levels = target_parties)
    )
  
  ggplot(overlap_df, aes(x = Target, y = Source, fill = Disagreement)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.0f%%", Disagreement)), 
              size = 3, color = "white", fontface = "bold") +
    scale_fill_gradientn(
      colors = c(
        "#1a9850",  # Dark green: 0-5% (excellent agreement)
        "#91cf60",  # Light green: 5-10% (good)
        "#d9ef8b",  # Yellow-green: 10-15% (acceptable)
        "#ffffbf",  # Yellow: 15-20% (borderline)
        "#fee08b",  # Light orange: 20-25% (concerning)
        "#fc8d59",  # Orange: 25-30% (problematic)
        "#d73027",  # Red: 30-40% (severe)
        "#a50026",  # Dark red: 40-50% (extreme)
        "#67001f"   # Very dark red: 50%+ (critical)
      ),
      values = scales::rescale(c(0, 5, 10, 15, 20, 25, 30, 40, 50)),
      name = "CI\nDisagreement\n(%)",
      limits = c(0, 50),
      breaks = c(0, 10, 20, 30, 40, 50),
      labels = c("0%", "10%", "20%", "30%", "40%", "50%+")
    ) +
    labs(
      title = "Credible Interval Disagreement Across Model Runs",
      subtitle = "Percentage of pairwise run comparisons with non-overlapping 95% CIs",
      x = "Target Party",
      y = "Source Party",
      caption = "Higher values indicate greater model instability"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(size = 14, face = "bold"),
      panel.grid = element_blank()
    )
}

# Updated combined heatmap with CV % scale
create_combined_uncertainty_heatmap <- function(var_results, beg_yr, end_yr) {
  source_parties <- rownames(var_results$individual_runs[[1]]$posterior_mean)
  target_parties <- colnames(var_results$individual_runs[[1]]$posterior_mean)
  
  n_runs <- length(var_results$individual_runs)
  n_source <- length(source_parties)
  n_target <- length(target_parties)
  
  # Calculate metrics for each transition
  uncertainty_data <- tibble()
  
  for (i in 1:n_source) {
    for (j in 1:n_target) {
      # Between-run percentage
      within_var <- mean(sapply(var_results$individual_runs, function(run) {
        run$posterior_sd[i, j]^2
      }))
      means <- sapply(var_results$individual_runs, function(run) {
        run$posterior_mean[i, j]
      })
      between_var <- var(means)
      total_var <- within_var + between_var
      pct_between <- if (total_var > 0) {
        100 * between_var / total_var
      } else {
        0
      }
      
      # CI disagreement (inverted)
      overlaps <- 0
      comparisons <- 0
      for (r1 in 1:(n_runs-1)) {
        for (r2 in (r1+1):n_runs) {
          ci1_lower <- var_results$individual_runs[[r1]]$posterior_q025[i, j]
          ci1_upper <- var_results$individual_runs[[r1]]$posterior_q975[i, j]
          ci2_lower <- var_results$individual_runs[[r2]]$posterior_q025[i, j]
          ci2_upper <- var_results$individual_runs[[r2]]$posterior_q975[i, j]
          
          if ((ci1_upper >= ci2_lower) && (ci1_lower <= ci2_upper)) {
            overlaps <- overlaps + 1
          }
          comparisons <- comparisons + 1
        }
      }
      ci_disagreement <- 100 - 100 * overlaps / comparisons
      
      # CV as percentage
      cv_pct <- var_results$variability$meta_cv[i, j] * 100
      
      uncertainty_data <- bind_rows(uncertainty_data, tibble(
        Source = source_parties[i],
        Target = target_parties[j],
        `Between-Run %` = pct_between,
        `CI Disagreement %` = ci_disagreement,
        `CV %` = cv_pct
      ))
    }
  }
  
  # Find the actual max CV value
  max_cv <- max(uncertainty_data$`CV %`, na.rm = TRUE)
  
  # Reshape for faceting
  plot_data <- uncertainty_data %>%
    pivot_longer(cols = c(`Between-Run %`, `CI Disagreement %`, `CV %`),
                 names_to = "Metric", values_to = "Value") %>%
    mutate(
      Source = factor(Source, levels = source_parties),
      Target = factor(Target, levels = target_parties),
      Metric = factor(Metric, levels = c("Between-Run %", "CI Disagreement %", "CV %"))
    )
  
  # Cap extremely high values for color scale but show actual in text
  plot_data <- plot_data %>%
    mutate(
      # Cap at 150 for color purposes
      Value_for_color = if_else(Value > 150, 150, Value),
      # Show actual value with indicator if >150
      Label = if_else(Value > 150, 
                      sprintf("%.0f*", Value),
                      sprintf("%.0f", Value))
    )
  
  # Create faceted plot with extended CV scale
  ggplot(plot_data, aes(x = Target, y = Source, fill = Value_for_color)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = Label), 
              size = 2.5, color = "white", fontface = "bold") +
    facet_wrap(~Metric, ncol = 3) +
    scale_fill_gradientn(
      colors = c(
        "#1a9850",  # Dark green: 0-15
        "#91cf60",  # Light green: 15-25
        "#d9ef8b",  # Yellow-green: 25-40
        "#ffffbf",  # Yellow: 40-50
        "#fee08b",  # Light orange: 50-65
        "#fc8d59",  # Orange: 65-80
        "#d73027",  # Red: 80-95
        "#a50026",  # Dark red: 95-120
        "#67001f"   # Very dark red: 120-150+
      ),
      values = scales::rescale(c(0, 15, 25, 40, 50, 65, 80, 95, 120, 150)),
      name = "Value",
      limits = c(0, 150),  # Hard cap at 150
      oob = scales::squish,  # Force values >150 to use the highest color
      breaks = c(0, 30, 60, 90, 120, 150),
      labels = c("0", "30", "60", "90", "120", "150+")
    ) +
    labs(
      title = paste("Uncertainty Decomposition between", beg_yr, "and", end_yr),
      x = paste(as.character(end_yr), "Party"),
      y = paste(as.character(beg_yr), "Party"),
      caption = "All metrics: Higher = worse. *Values >150% shown with asterisk and capped at darkest red."
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = element_text(size = 8),
      plot.title = element_text(size = 14, face = "bold"),
      strip.text = element_text(size = 10, face = "bold"),
      panel.grid = element_blank()
    )
}

# And the summary table
create_uncertainty_summary_table <- function(var_results) {
  source_parties <- rownames(var_results$individual_runs[[1]]$posterior_mean)
  target_parties <- colnames(var_results$individual_runs[[1]]$posterior_mean)
  
  n_runs <- length(var_results$individual_runs)
  
  summary_table <- tibble()
  
  for (i in seq_along(source_parties)) {
    for (j in seq_along(target_parties)) {
      # Calculate metrics
      within_var <- mean(sapply(var_results$individual_runs, function(run) {
        run$posterior_sd[i, j]^2
      }))
      means <- sapply(var_results$individual_runs, function(run) {
        run$posterior_mean[i, j]
      })
      between_var <- var(means)
      total_var <- within_var + between_var
      pct_between <- if (total_var > 0) {
        100 * between_var / total_var
      } else {
        0
      }
      
      summary_table <- bind_rows(summary_table, tibble(
        Transition = paste(source_parties[i], "→", target_parties[j]),
        `Mean Est. (%)` = var_results$variability$meta_mean[i, j] * 100,
        `CV` = var_results$variability$meta_cv[i, j],
        `Between-Run (%)` = pct_between,
        `Assessment` = case_when(
          pct_between > 60 ~ "Very High Instability",
          pct_between > 40 ~ "High Instability",
          pct_between > 20 ~ "Moderate Instability",
          TRUE ~ "Low Instability"
        )
      ))
    }
  }
  
  # Sort by between-run percentage
  summary_table <- summary_table %>%
    arrange(desc(`Between-Run (%)`))
  
  return(summary_table)
}

# Create a heatmap of the range of uncertainty across all runs for each transition
create_transition_uncertainty_heatmap <- function(all_flows) {
  # all_flows: data frame with columns 'from', 'to', and 'range_pct' (percentage points)
  
  # Prepare data for heatmap
  heatmap_data <- all_flows %>%
    mutate(
      from_label = gsub("_", " ", gsub("_in_185.", "", from)),
      to_label = gsub("_", " ", gsub("_in_185.", "", to)),
      from_label = ifelse(from_label == "Abstaining", "Non-voting", from_label),
      to_label = ifelse(to_label == "Abstaining", "Non-voting", to_label),
      from_label = fct_relevel(from_label, party_sort),
      to_label = fct_relevel(to_label, party_sort),
      from_year = gsub(".*_in_(185.)", "\\1", from),
      to_year = gsub(".*_in_(185.)", "\\1", to),
      range_pct_display = range_pct * 100  # Convert to percentage points
    )
  
  from_year <- unique(heatmap_data$from_year)
  to_year <- unique(heatmap_data$to_year)
  
  # Create heatmap of uncertainty ranges
  ggplot(heatmap_data, aes(x = to_label, y = from_label)) +
    geom_tile(aes(fill = range_pct_display), color = "white", size = 1) +
    geom_text(aes(label = sprintf("±%.1f", range_pct_display / 2)), 
              size = 3) +
    scale_fill_gradient2(
      low = "#2166ac",
      mid = "#ffffbf",
      high = "#d73027",
      midpoint = 5,  # Adjust based on your data
      name = "Range\n(pp)",
      labels = function(x) sprintf("±%.0f", x/2)
    ) +
    scale_x_discrete(expand = expansion(add = 0.5)) +
    scale_y_discrete(expand = expansion(add = 0.5)) +
    labs(
      title = "Estimate Uncertainty",
      subtitle = "Range across 10 model runs",
      x = paste0("Vote in ", to_year),
      y = paste0("Vote in ", from_year)
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold"),
      axis.title.y = element_blank()
    )
}

# Helper function to extract R-hat values by transition type and summarize across towns 
rhat_by_transition <- function(rhat_values) {
  # rhat_values: data frame with columns 'variable' (beta variable) and 'rhat' (R-hat value)
  
  rhat_by_transition <- rhat_values %>%
    mutate(
      # Extract transition information from variable names
      # e.g., "beta.Democrat_in_1851.Democrat_in_1852.1" 
      variable_clean = gsub("beta\\.", "", variable),
      # Extract from, to, and town_id
      from = sub("\\..*", "", variable_clean),
      rest = sub("^[^.]+\\.", "", variable_clean),
      to = sub("\\..*", "", rest),
      from_year = gsub(".*_in_(185.)", "\\1", from),
      to_year = gsub(".*_in_(185.)", "\\1", to),
      town_id = as.numeric(sub(".*\\.", "", rest)),
      # Clean up labels
      from_label = gsub("_in_\\d+", "", from),
      to_label = gsub("_in_\\d+", "", to),
      from_label = gsub("_", " ", from_label),
      to_label = gsub("_", " ", to_label),
      from_label = ifelse(from_label == "Abstaining", "Non-voting", from_label),
      to_label = ifelse(to_label == "Abstaining", "Non-voting", to_label),
      from_label = fct_relevel(from_label, party_sort),
      to_label = fct_relevel(to_label, party_sort)
    ) %>%
    group_by(from_label, to_label) %>%
    summarize(
      max_rhat = max(rhat, na.rm = TRUE),
      mean_rhat = mean(rhat, na.rm = TRUE),
      from_year = unique(from_year),
      to_year = unique(to_year),
      .groups = "drop"
    )
  
  return(rhat_by_transition)
}

# Create a heatmap of the R-hat values for each transition across all runs
create_rhat_heatmap <- function(rhat_values) {
  # rhat_values: data frame with columns 'from', 'to', and 'rhat' (R-hat value)
  
  # Prepare data for heatmap
  rhat_by_transition <- rhat_by_transition(rhat_values)
  
  from_year <- unique(rhat_by_transition$from_year)
  to_year <- unique(rhat_by_transition$to_year)
  
  # Show max R-hat
  ggplot(rhat_by_transition, aes(x = to_label, y = from_label)) +
    geom_tile(aes(fill = max_rhat), color = "white", size = 1) +
    geom_text(aes(label = sprintf("%.2f", max_rhat)), 
              size = 3) +
    scale_fill_gradient2(
      low = "#2166ac",
      mid = "#ffffbf",
      high = "#d73027",
      midpoint = 1.05,
      name = "Max R-hat"
    ) +
    scale_x_discrete(expand = expansion(add = 0.5)) +
    scale_y_discrete(expand = expansion(add = 0.5)) +
    labs(
      title = "Maximum Convergence",
      subtitle = "Max R-hat across towns",
      x = paste0("Vote in ", to_year),
      y = paste0("Vote in ", from_year)
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold")
    )
}

# Create a town-level map of the maximum R-hat value across all transitions
create_rhat_town_map <- function(rhat_values, map_data) {
  # rhat_values: data frame with columns 'variable' (beta variable) and 'rhat' (R-hat value)
  # map_data: sf object with town geometries and a 'town_id' column
  
  rhat_by_town <- rhat_values %>%
    mutate(
      # Extract transition information from variable names
      # e.g., "beta.Democrat_in_1851.Democrat_in_1852.1" 
      variable_clean = gsub("beta\\.", "", variable),
      # Extract from, to, and town_id
      from = sub("\\..*", "", variable_clean),
      rest = sub("^[^.]+\\.", "", variable_clean),
      to = sub("\\..*", "", rest),
      from_year = gsub(".*_in_(185.)", "\\1", from),
      to_year = gsub(".*_in_(185.)", "\\1", to),
      town_id = as.numeric(sub(".*\\.", "", rest)),
      # Clean up labels
      from_label = gsub("_in_\\d+", "", from),
      to_label = gsub("_in_\\d+", "", to),
      from_label = gsub("_", " ", from_label),
      to_label = gsub("_", " ", to_label),
      from_label = ifelse(from_label == "Abstaining", "Non-voting", from_label),
      to_label = ifelse(to_label == "Abstaining", "Non-voting", to_label),
      from_label = fct_relevel(from_label, party_sort),
      to_label = fct_relevel(to_label, party_sort)
    ) %>%
    group_by(town_id) %>%
    summarize(
      max_rhat = max(rhat, na.rm = TRUE),
      mean_rhat = mean(rhat, na.rm = TRUE),
      from_year = unique(from_year),
      to_year = unique(to_year),
      .groups = "drop"
    )
  
  from_year <- unique(rhat_by_town$from_year)
  to_year <- unique(rhat_by_town$to_year)
  
  map_data <- map_data %>%
    bind_cols(rhat_by_town)
  
  ggplot(map_data) +
    geom_sf(aes(fill = max_rhat), color = "white") +
    scale_fill_gradient2(
      low = "#2166ac",
      mid = "#ffffbf",
      high = "#d73027",
      midpoint = 1.05,
      name = "Max R-hat"
    ) +
    labs(
      title = (paste0("Maximum Convergence by Town (", from_year, "→", to_year, ")")),
      subtitle = "Maximum R-hat across all transitions",
      fill = "Max R-hat"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold")
    )
}
