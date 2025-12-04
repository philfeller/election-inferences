# Sankey / alluvial pipeline for cov.ei.* models
# Assumes models cov.ei.52 .. cov.ei.57 exist in your R session and each has:
#   model$draws$Cell.counts  (matrix: draws x n_cells) with colnames like
#   "ccount.Source_in_YYYY.Target_in_YYYY"
#
# Requires: dplyr, tidyr, purrr, networkD3, ggalluvial, ggplot2
#
# Run the whole script; see comments for options.

library(dplyr)
library(tidyr)
library(purrr)
library(networkD3)
library(ggalluvial)
library(ggplot2)

# ---- 1. collect models (adjust names if needed) ----
model_names <- c("cov.ei.52", "cov.ei.53", "cov.ei.54", "cov.ei.55", "cov.ei.56", "cov.ei.57")
models <- mget(model_names, ifnotfound = NA, envir = .GlobalEnv)
if (any(is.na(models))) stop("One or more model names not found in global env. Check model_names.")

# ---- 2. helper: extract per-draw aggregated cell counts and parse cell names ----
extract_cell_draws <- function(em) {
  cc <- as.matrix(em$draws$Cell.counts)  # draws x n_cells
  cell_names <- colnames(cc)
  if (is.null(cell_names)) stop("No colnames on Cell.counts; cannot parse cells.")
  parts <- strsplit(cell_names, "\\.")
  # Expect parts[[i]] like c("ccount", "Free_Soil_in_1854", "Know_Nothing_in_1855")
  source <- sapply(parts, `[`, 2)
  target <- sapply(parts, `[`, 3)
  df <- as.data.frame(cc)
  names(df) <- cell_names
  # tidy: draws x cell -> long with draw index
  draws <- seq_len(nrow(cc))
  long <- as_tibble(df) %>%
    mutate(.draw = draws) %>%
    pivot_longer(-.draw, names_to = "cell", values_to = "count") %>%
    mutate(source = source[match(cell, cell_names)],
           target = target[match(cell, cell_names)])
  return(list(long = long, cell_names = cell_names, source = source, target = target,
              draws = nrow(cc), n_cells = length(cell_names)))
}

# extract for each model, compute per-draw conditional probs p(target | source)
pair_data <- map(models, function(em) {
  ex <- extract_cell_draws(em)
  long <- ex$long
  # compute conditional p(target | source) per draw: normalize counts within each draw and source
  long <- long %>%
    group_by(.draw, source) %>%
    mutate(source_total = sum(count),
           prop = ifelse(source_total > 0, count / source_total, 0)) %>%
    ungroup()
  # Summarize mean and CI across draws per source->target cell
  cell_summary <- long %>%
    group_by(source, target) %>%
    summarize(mean_count = mean(count),
              mean_prop  = mean(prop),
              prop_lo = quantile(prop, 0.025),
              prop_hi = quantile(prop, 0.975),
              .groups = "drop")
  list(long = long, summary = cell_summary)
})

names(pair_data) <- names(models)

# ---- 3. build nodes & links for Sankey (mean flows) ----
# nodes: unique party-year labels across all pairs
all_nodes <- map(pair_data, ~.$summary %>% select(source, target)) %>%
  bind_rows() %>%
  pivot_longer(everything(), values_to = "node") %>%
  distinct(node) %>%
  arrange(node)

nodes <- tibble(name = all_nodes$node)

# links: combine summaries from each pair
links_df <- map_dfr(pair_data, ~.$summary %>% select(source, target, mean_prop)) %>%
  mutate(source_id = match(source, nodes$name) - 1,
         target_id = match(target, nodes$name) - 1,
         value = mean_prop)

# optionally filter tiny links to reduce clutter:
min_prop_display <- 0.005  # adjust
links_df_filtered <- links_df %>% filter(value >= min_prop_display)

# networkD3 sankey (interactive)
sankey <- sankeyNetwork(Links = links_df_filtered %>% select(source = source_id, target = target_id, value),
                        Nodes = nodes, Source = "source", Target = "target",
                        Value = "value", NodeID = "name",
                        fontSize = 12, nodeWidth = 30)
# Print sankey in RStudio viewer or save htmlwidget:
print(sankey)

# ---- 4. approximate multi-year alluvial by chaining pairwise mean conditional probs ----
# This approximation multiplies conditional means across pairs to get path weights (fast but ignores posterior dependence).
# Steps: create list of dataframes with columns source/year and target/year and mean_prop; then join sequentially.

# First, we need to order pairs by source-year. We'll parse year from label assuming pattern like "Party_in_YYYY"
extract_year <- function(label) {
  m <- regmatches(label, regexpr("\\d{4}", label))
  if (length(m) == 0) return(NA_character_)
  return(m)
}

# Build a small helper to make pair_df with columns: source_label, target_label, mean_prop, source_year, target_year
pair_summary_byname <- map2(names(pair_data), pair_data, ~ {
  df <- .y$summary %>% rename(source_label = source, target_label = target, mean = mean_prop)
  df <- df %>% mutate(source_year = extract_year(source_label),
                      target_year = extract_year(target_label),
                      pair_name = .x)
  df
})

# Order pairs by source_year so we chain correctly
pair_summary_byname <- pair_summary_byname %>% 
  set_names(names(pair_data)) %>%
  enframe(name = "pair_name", value = "df") %>%
  mutate(src_year = map_chr(df, ~ unique(.x$source_year)[1])) %>%
  arrange(src_year)

# If the chain is contiguous (e.g., 1851->1852,1852->1854,1854->1855), adapt accordingly.
# For an example 4-step chain: A -> B -> C -> D, we can do inner joins to form A-B-C-D table.

# Extract the dfs in order
dfs_ordered <- pair_summary_byname$df

# If you have exactly 3+ pairs, you can chain them. Here's generic join for N pairs:
# We'll build successive joins to form multiple columns A,B,C,... representing party at each step.
chain_paths <- function(dfs) {
  # rename columns to Level1, Level2,... for successive joins
  renamed <- map2(dfs, seq_along(dfs), function(df, k) {
    df %>% rename(!!paste0("L", k) := source_label,
                  !!paste0("L", k+1) := target_label,
                  !!paste0("p", k) := mean) %>%
      select(starts_with("L"), starts_with("p"))
  })
  # iterative inner join on overlapping L columns (L2 of first with L2 of second, etc.)
  res <- renamed[[1]]
  if (length(renamed) > 1) {
    for (k in 2:length(renamed)) {
      # join res (has Lk) to renamed[[k]] on column Lk
      join_col <- paste0("L", k)
      res <- inner_join(res, renamed[[k]], by = join_col)
    }
  }
  # compute path weight as product of p1*p2*...
  pcols <- grep("^p", names(res), value = TRUE)
  res <- res %>% mutate(weight = apply(select(., all_of(pcols)), 1, prod))
  return(res)
}

# Build approximate chained paths
paths_approx <- chain_paths(dfs_ordered)

# create ggalluvial plot if there are <= 6 levels (aesthetics)
num_levels <- length(dfs_ordered) + 1
if (num_levels <= 6) {
  # rename to human-friendly axis labels
  axis_labels <- str_replace(names(select(paths_approx, starts_with("L"))), "L", "Year")
  ggplot(paths_approx,
         aes_string(axis1 = "L1", axis2 = "L2",
                    x = NULL)) # we'll add aesthetics below
  # build alluvium aesthetics dynamically:
  plot_data <- paths_approx %>% mutate(weight = weight)
  # Build formula for n levels
  aes_args <- setNames(paste0("L", 1:num_levels), paste0("axis", 1:num_levels))
  ggplot(plot_data,
         aes_string(x = "1", alluvium = paste0("interaction(", paste0("L", 1:num_levels, collapse = ","), ")"),
                    y = "weight", fill = "L1")) +
    geom_flow(stat = "identity", lode.guidance = "frontback", alpha = 0.8) +
    geom_stratum(width = 0.2) +
    theme_minimal() +
    ggtitle("Approximate multi-year flows (product of pairwise conditional means)") +
    theme(axis.text.y = element_blank(), axis.ticks = element_blank())
} else {
  message("Too many levels for simple ggalluvial plotting; consider networkD3 or subsetting parties.")
}

# ---- 5. Posterior-chained full paths (more accurate) ----
# This computes path weights per draw by multiplying per-draw conditional probabilities across pairs,
# then summarizes across draws to get mean and CI for each full path.
# WARNING: number of possible full paths = product(#parties at each year). It can explode.
#
# Steps:
# - For each pair, create an array draws x source x target of per-draw p(target|source)
# - Build list of unique parties per year (year sequence inferred from pairs)
# - Enumerate all possible paths (cartesian product). If total paths > 20000, the enumeration may be slow.
# - For each draw s, compute path weight as product of the appropriate p's for that path.
# - Summarize path weights across draws.

build_per_draw_cond_array <- function(pair_long_df) {
  # pair_long_df: from extract_cell_draws, with .draw, source, target, prop
  draws <- max(pair_long_df$.draw)
  sources <- unique(pair_long_df$source)
  targets <- unique(pair_long_df$target)
  # create array draws x source x target where missing combinations get 0
  arr <- array(0, dim = c(draws, length(sources), length(targets)),
               dimnames = list(draw = seq_len(draws), source = sources, target = targets))
  for (d in seq_len(draws)) {
    dd <- pair_long_df %>% filter(.draw == d)
    for (r in seq_len(nrow(dd))) {
      s <- dd$source[r]; t <- dd$target[r]; p <- dd$prop[r]
      arr[d, as.character(s), as.character(t)] <- p
    }
  }
  return(list(arr = arr, sources = sources, targets = targets))
}

# build per-draw arrays for each pair in order
pair_arrays <- map(pair_data, ~ build_per_draw_cond_array(.x$long))

# Determine year ordering and species per timepoint
# For each pair, extract source_year and target_year from labels
get_year_label <- function(pair_df) {
  # pair_df = .$summary with source and target
  all_sources <- unique(pair_df$source)
  yrs <- sort(unique(extract_year(all_sources)))
  return(yrs)
}
# Build a list of unique party labels per year by scanning all pairs
party_by_year <- list()
for (i in seq_along(pair_arrays)) {
  sources <- dimnames(pair_arrays[[i]]$arr)$source
  targets <- dimnames(pair_arrays[[i]]$arr)$target
  # extract year strings
  syears <- unique(extract_year(sources))
  tyears <- unique(extract_year(targets))
  if (!is.na(syears)) party_by_year[[syears]] <- unique(c(party_by_year[[syears]], sources))
  if (!is.na(tyears)) party_by_year[[tyears]] <- unique(c(party_by_year[[tyears]], targets))
}
# sort years
years_sorted <- sort(names(party_by_year))
message("Years found (inferred): ", paste(years_sorted, collapse = ", "))

# build list of party vectors in chronological order
party_lists <- map(years_sorted, ~ party_by_year[[.x]])

# compute number of paths
num_paths <- prod(map_int(party_lists, length))
if (num_paths > 20000) {
  message("WARNING: number of full paths = ", num_paths, " (may be slow). Consider subsetting parties.")
}

# enumerate all paths (cartesian product)
paths_df <- expand.grid(!!!set_names(party_lists, paste0("Y", seq_along(party_lists))),
                        stringsAsFactors = FALSE) %>%
  as_tibble()

# For each draw, compute path weight:
n_draws <- dim(pair_arrays[[1]]$arr)[1]
compute_path_weights_per_draw <- function(draw_idx) {
  # draw_idx: scalar
  # for each pair k (between year k and k+1), get per-draw p(target|source)
  # we need to map Lk (party at year k) -> Lk+1 (party at year k+1) probabilities
  pvec <- numeric(nrow(paths_df))
  for (rowi in seq_len(nrow(paths_df))) {
    prob <- 1
    for (k in seq_len(length(pair_arrays))) {
      # pair_arrays[[k]] corresponds to pair k -> k+1 in the ordered sequence of pairs we gave earlier (models order)
      arr <- pair_arrays[[k]]$arr
      src_label <- paths_df[[paste0("Y", k)]][rowi]
      tgt_label <- paths_df[[paste0("Y", k+1)]][rowi]
      # if src or tgt not present in this pair's dimnames, treat p=0
      if (!(src_label %in% dimnames(arr)$source) || !(tgt_label %in% dimnames(arr)$target)) {
        prob <- 0; break
      } else {
        prob <- prob * arr[draw_idx, src_label, tgt_label]
      }
    }
    pvec[rowi] <- prob
  }
  return(pvec)
}

# run across draws (careful: this loop can be slow)
message("Computing per-draw path weights for ", n_draws, " draws and ", nrow(paths_df), " paths...")
path_weights_by_draw <- replicate(n_draws, numeric(nrow(paths_df)))  # will be overwritten; we reassign below

for (s in seq_len(n_draws)) {
  path_weights_by_draw[, s] <- compute_path_weights_per_draw(s)
}

# summarize across draws: mean and 95% CI
path_mean <- rowMeans(path_weights_by_draw)
path_lo   <- apply(path_weights_by_draw, 1, quantile, 0.025)
path_hi   <- apply(path_weights_by_draw, 1, quantile, 0.975)
paths_df$mean <- path_mean
paths_df$lo <- path_lo
paths_df$hi <- path_hi

# Sort and show top paths
top_paths <- paths_df %>% arrange(desc(mean)) %>% slice_head(n = 30)
print(top_paths)

# Optionally, build an alluvial plot from the posterior-mean paths (this is accurate under chaining)
# Use only top-k paths to avoid clutter:
topk <- 50
plot_paths <- top_paths %>% slice_head(n = topk)

# For plotting need columns in order:
if (ncol(plot_paths) >= 4) {
  # convert to long format for ggalluvial
  plot_long <- plot_paths %>%
    select(starts_with("Y"), mean) %>%
    mutate(id = row_number()) %>%
    pivot_longer(cols = starts_with("Y"), names_to = "year_idx", values_to = "party") %>%
    mutate(year = as.integer(gsub("Y", "", year_idx)))
  ggplot(plot_long,
         aes(x = factor(year), stratum = party, alluvium = id, y = plot_paths$mean[match(id, plot_paths$row_number)],
             fill = party)) +
    geom_flow(stat = "identity") +
    geom_stratum() +
    theme_minimal() +
    ggtitle("Top posterior-mean full paths (chained per-draw multiplications)")
} else {
  message("Not enough levels for multi-year alluvial plot.")
}

# ---- DONE ----
# The script produced:
# - sankey interactive widget (print(sankey))
# - approximate chained alluvial (if small)
# - posterior-chained full paths summarized and top paths printed
#
# Interpretations & caveats:
# - The 'approximate' approach multiplies pairwise conditional MEANS. That is fast but biases path means when the p's are variable (E[prod] != prod(E)).
# - The 'posterior-chained' approach multiplies per-draw conditional probabilities then averages: it estimates E[ product_k p_k ] and is coherent with the posterior if you assume independence across pairs' draws or align draws as an approximation.
# - If models were fitted jointly over multiple years you'd get coherent draws; if fit separately, the per-draw multiplication mixes independent posteriors — still useful as an exploratory summary but interpret cautiously.
#
# If you want, I can:
# - adjust plotting aesthetics, colors, and node ordering,
# - produce HTML output of the sankey,
# - reduce dimensionality by collapsing small parties into "Other",
# - or produce an animation over sampled posterior draws to show uncertainty visually.