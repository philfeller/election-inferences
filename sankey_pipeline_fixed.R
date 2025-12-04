library(dplyr)
library(tidyr)
library(purrr)
library(networkD3)
library(ggalluvial)
library(ggplot2)
library(stringr)

# ---- 1. collect models (adjust names if needed) ----
model_names <- c("cov.ei.52", "cov.ei.53", "cov.ei.54", "cov.ei.55", "cov.ei.56", "cov.ei.57")
models <- mget(model_names, ifnotfound = NA, envir = .GlobalEnv)
if (any(is.na(models))) stop("One or more model names not found in global env. Check model_names.")

# ---- helper: clean party label (strip '_in_185x' suffix and convert underscores to spaces) ----
clean_party_label <- function(full_label) {
  # Remove trailing "_in_YYYY" (any 4-digit year) and replace underscores with spaces
  lbl <- sub("_in_\\d{4}$", "", full_label)
  lbl <- gsub("_", " ", lbl)
  trimws(lbl)
}

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
# nodes: unique party-year labels across all pairs (keep full label internally for uniqueness)
all_nodes <- map(pair_data, ~.$summary %>% select(source, target)) %>%
  bind_rows() %>%
  pivot_longer(everything(), values_to = "node") %>%
  distinct(node) %>%
  arrange(node)

nodes <- tibble(full_name = all_nodes$node)
# Add display label (just party name) and group (party) for colors
nodes <- nodes %>%
  mutate(party = clean_party_label(full_name),
         label = party)  # label will be used for display

# If the user provided a named vector party_colors (names = party), use it; otherwise create a default palette
# Example expected format: party_colors <- c("Democrat" = "#1f77b4", "Whig" = "#ff7f0e", ...)
if (exists("party_colors") && is.character(party_colors) && !is.null(names(party_colors))) {
  # ensure all parties have a color; if missing, assign a default
  missing_parties <- setdiff(unique(nodes$party), names(party_colors))
  if (length(missing_parties) > 0) {
    default_palette <- grDevices::rainbow(length(missing_parties))
    names(default_palette) <- missing_parties
    party_colors <- c(party_colors, default_palette)
  }
} else {
  unique_parties <- unique(nodes$party)
  assigned <- grDevices::rainbow(length(unique_parties))
  names(assigned) <- unique_parties
  party_colors <- assigned
}

# attach group (party) to nodes for colour mapping
nodes <- nodes %>% mutate(group = party)

# links: combine summaries from each pair
links_df <- map_dfr(pair_data, ~.$summary %>% select(source, target, mean_prop)) %>%
  mutate(source_id = match(source, nodes$full_name) - 1,
         target_id = match(target, nodes$full_name) - 1,
         value = mean_prop)

# optionally filter tiny links to reduce clutter:
min_prop_display <- 0.005  # adjust
links_df_filtered <- links_df %>% filter(value >= min_prop_display)

# Build JS colourScale for networkD3 using party_colors (domain = party names, range = hex colors)
party_names_js <- jsonlite::toJSON(names(party_colors))
party_colors_js <- jsonlite::toJSON(unname(party_colors))
colourScale <- sprintf("d3.scaleOrdinal().domain(%s).range(%s)", party_names_js, party_colors_js)

# networkD3 sankey (interactive)
# node data must include 'name' column for networkD3; we'll supply 'label' for display, and NodeGroup = "group"
nodes_for_widget <- nodes %>% mutate(name = full_name)  # keep full_name as unique internal id
# networkD3 expects the Nodes data frame to have a column that matches NodeID; we'll use 'label' so labels show as party names
sankey <- sankeyNetwork(
  Links = links_df_filtered %>% select(source = source_id, target = target_id, value),
  Nodes = nodes_for_widget,
  Source = "source", Target = "target",
  Value = "value", NodeID = "label", NodeGroup = "group",
  fontSize = 12, nodeWidth = 30,
  colourScale = htmlwidgets::JS(colourScale)
)
# Print sankey in RStudio viewer or save htmlwidget:
print(sankey)

# ---- 4. corrected chain_paths implementation (fix: compute numeric product robustly) ----
chain_paths <- function(dfs) {
  # rename columns to Level1, Level2,... and create p1, p2, ...
  renamed <- map2(dfs, seq_along(dfs), function(df, k) {
    df %>%
      rename(!!paste0("L", k) := source_label,
             !!paste0("L", k+1) := target_label,
             !!paste0("p", k) := mean) %>%
      select(starts_with("L"), starts_with("p"))
  })
  # iterative inner join on overlapping L columns (L2 of first with L2 of second, etc.)
  res <- renamed[[1]]
  if (length(renamed) > 1) {
    for (k in 2:length(renamed)) {
      join_col <- paste0("L", k)
      res <- inner_join(res, renamed[[k]], by = join_col)
    }
  }
  # robust numeric product across p columns: use purrr::pmap_dbl to avoid character coercion issues
  pcols <- grep("^p", names(res), value = TRUE)
  if (length(pcols) == 0) stop("No p columns found in chain_paths result.")
  res <- res %>%
    mutate(across(all_of(pcols), ~ as.numeric(.))) %>%               # ensure numeric
    mutate(weight = pmap_dbl(select(., all_of(pcols)), function(...) prod(unlist(list(...)), na.rm = FALSE)))
  return(res)
}

# Example use (after building pair_summary_byname as earlier in the pipeline):
# dfs_ordered <- pair_summary_byname$df
# paths_approx <- chain_paths(dfs_ordered)
# (the rest of the plotting/processing continues unchanged)