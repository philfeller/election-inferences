# Sankey pipeline: use mean counts and preserve unique nodes while displaying only the party name.
# Run in the same session as cov.ei.52..cov.ei.57
#
# Requirements: dplyr, tidyr, purrr, networkD3, ggplot2, stringr, jsonlite, htmlwidgets
#
library(dplyr)
library(tidyr)
library(purrr)
library(networkD3)
library(ggplot2)
library(stringr)
library(jsonlite)
library(htmlwidgets)

# ---- collect models ----
model_names <- c("cov.ei.52", "cov.ei.53", "cov.ei.54", "cov.ei.55", "cov.ei.56", "cov.ei.57")
models <- mget(model_names, ifnotfound = NA, envir = .GlobalEnv)
if (any(is.na(models))) stop("One or more model names not found in global env. Check model_names.")

# ---- helper: clean party label (strip '_in_185x' suffix and convert underscores to spaces) ----
clean_party_label <- function(full_label) {
  lbl <- sub("_in_\\d{4}$", "", full_label)
  lbl <- gsub("_", " ", lbl)
  trimws(lbl)
}

# ---- helper: extract per-draw aggregated cell counts and parse cell names ----
extract_cell_draws <- function(em) {
  cc <- as.matrix(em$draws$Cell.counts)  # draws x n_cells
  cell_names <- colnames(cc)
  if (is.null(cell_names)) stop("No colnames on Cell.counts; cannot parse cells.")
  parts <- strsplit(cell_names, "\\.")
  source <- sapply(parts, `[`, 2)
  target <- sapply(parts, `[`, 3)
  df <- as.data.frame(cc)
  names(df) <- cell_names
  draws <- seq_len(nrow(cc))
  long <- as_tibble(df) %>%
    mutate(.draw = draws) %>%
    pivot_longer(-.draw, names_to = "cell", values_to = "count") %>%
    mutate(source = source[match(cell, cell_names)],
           target = target[match(cell, cell_names)])
  list(long = long, cell_names = cell_names, source = source, target = target,
       draws = nrow(cc), n_cells = length(cell_names))
}

# ---- extract and summarize (mean_count, mean_prop, CI) for each pair ----
pair_data <- map(models, function(em) {
  ex <- extract_cell_draws(em)
  long <- ex$long
  long <- long %>%
    group_by(.draw, source) %>%
    mutate(source_total = sum(count),
           prop = ifelse(source_total > 0, count / source_total, 0)) %>%
    ungroup()
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

# ---- build nodes ----
all_nodes <- map(pair_data, ~.$summary %>% select(source, target)) %>%
  bind_rows() %>%
  pivot_longer(everything(), values_to = "node") %>%
  distinct(node) %>%
  arrange(node)

nodes <- tibble(full_name = all_nodes$node) %>%
  mutate(party = clean_party_label(full_name))

# Create a visually identical display label for each node (party only) but keep NodeIDs unique.
# We append a zero-width space + index to NodeID. The zero-width space is invisible but makes string unique.
zw <- "\u200B"
nodes <- nodes %>%
  mutate(name_id = paste0(party, zw, seq_len(n())),
         label_display = party,
         group = party)  # group used for colour mapping

# ---- ensure party_colors maps exactly to nodes$group ----
if (exists("party_colors") && is.character(party_colors) && !is.null(names(party_colors))) {
  # normalize names for matching
  norm <- function(x) trimws(gsub("\\s+", " ", as.character(x)))
  det <- unique(nodes$group); det_norm <- norm(det)
  user_names <- names(party_colors); user_norm <- norm(user_names)
  map_idx <- match(det_norm, user_norm)
  matched <- party_colors[map_idx]
  names(matched) <- det
  missing <- which(is.na(matched) | matched == "")
  if (length(missing) > 0) {
    default_palette <- grDevices::rainbow(length(missing))
    matched[missing] <- default_palette
  }
  party_colors_matched <- matched
} else {
  unique_parties <- unique(nodes$group)
  assigned <- grDevices::rainbow(length(unique_parties))
  party_colors_matched <- setNames(assigned, unique_parties)
}

# Build JS colourScale for D3
party_names_js <- jsonlite::toJSON(names(party_colors_matched))
party_colors_js <- jsonlite::toJSON(unname(party_colors_matched))
colourScale <- sprintf("d3.scaleOrdinal().domain(%s).range(%s)", party_names_js, party_colors_js)

# ---- build links using mean_count so node widths reflect vote counts ----
links_df_raw <- map_dfr(pair_data, ~.$summary %>% select(source, target, mean_count)) %>%
  mutate(source_id = match(source, nodes$full_name) - 1,
         target_id = match(target, nodes$full_name) - 1,
         value = mean_count)

# If you want relative proportions of total votes instead of raw counts, normalize:
# total_votes_all <- sum(links_df_raw$value)
# links_df_raw <- links_df_raw %>% mutate(value = value / total_votes_all)

# optional: filter very tiny links relative to total (use fraction of total counts)
total_votes_all <- sum(links_df_raw$value)
min_frac_display <- 0.001  # show links >= 0.1% of total votes; adjust
links_df_filtered <- links_df_raw %>% filter(value >= min_frac_display * total_votes_all)

# ---- prepare nodes_for_widget and links_for_widget ----
nodes_for_widget <- nodes %>% transmute(name = name_id, label = label_display, group = group, full_name = full_name)
links_for_widget <- links_df_filtered %>% transmute(source = source_id, target = target_id, value = value)

# networkD3: use name (unique id) as NodeID but display label will be the visible text.
# networkD3 uses the column specified as NodeID to print labels; because we used zero-width characters
# the label looks like the party name but internal IDs are unique.
sankey_widget <- sankeyNetwork(
  Links = links_for_widget,
  Nodes = nodes_for_widget,
  Source = "source", Target = "target",
  Value = "value", NodeID = "name", NodeGroup = "group",
  fontSize = 12, nodeWidth = 30,
  colourScale = htmlwidgets::JS(colourScale)
)

# Inject a small JavaScript snippet to replace visible node text with the 'label' (party) without the zero-width suffix.
# This step is optional but makes node labels nicer (visual only).
sankey_widget <- htmlwidgets::onRender(
  sankey_widget,
  '
  function(el,x) {
    // find the node labels (text elements) and replace their textContent with the "label" field from nodes data (if present)
    // nodes data is available as x.nodes
    d3.select(el).selectAll(".node text").each(function(d,i) {
      // x.nodes[i] should correspond; use the label property if present otherwise fallback to current text
      if (x.nodes && x.nodes[i] && x.nodes[i].label) {
        d3.select(this).text(x.nodes[i].label);
      } else {
        // fallback: remove zero-width characters from current label
        var txt = d3.select(this).text();
        d3.select(this).text(txt.replace(/\\u200B/g,""));
      }
    });
  }
  '
)

# Print sankey
print(sankey_widget)

# ---- Quick diagnostics printed for verification ----
cat("Total votes (sum mean_count across all pairs):", total_votes_all, "\n")
cat("Top links (by mean_count):\n")
print(links_df_raw %>% arrange(desc(value)) %>% slice_head(n = 10) %>% select(source, target, value))
cat("\nNode parties and colors used:\n")
print(party_colors_matched)

# ---- Additionally: if you want the initial-year shares for a specific model (e.g., cov.ei.52), compute: ----
# Example for cov.ei.52:
if ("cov.ei.52" %in% names(models)) {
  df52 <- pair_data[["cov.ei.52"]]$summary
  source_totals_52 <- df52 %>% group_by(source) %>% summarize(total = sum(mean_count))
  # convert source label to cleaned party
  source_totals_52 <- source_totals_52 %>% mutate(party = clean_party_label(source))
  source_totals_52 <- source_totals_52 %>% mutate(share = total / sum(total))
  cat("\nInferred source shares for cov.ei.52 (should match get_shares(results.52,..) reasonably):\n")
  print(source_totals_52 %>% select(party, total, share))
}