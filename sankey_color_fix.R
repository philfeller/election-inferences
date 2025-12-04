# Fix for mismatched party colours in the sankey pipeline
#
# Run this after you run the pipeline that created `nodes` (the tibble with columns full_name, party, label).
# If you didn't keep `nodes`, re-run the node construction from the pipeline or run extract_cell_draws() +
# nodes <- tibble(full_name = all_nodes$node) %>% mutate(party = clean_party_label(full_name), label = party)
#
# This script:
# 1) shows the detected party names and the names supplied in `party_colors`
# 2) produces a matched party_colors vector (exactly matching the detected party names)
# 3) fills in any missing colours with a default palette
# 4) rebuilds the JS colourScale and returns a sankey widget using the corrected colours
#
# Usage:
# - If you have `nodes`, `links_df_filtered` and `links_df` in your session from the previous pipeline, just source() this file.
# - Otherwise re-create `nodes` as described and then run the functions below.

library(jsonlite)
library(htmlwidgets)

# Inspect current detected parties and provided party_colors (if any)
inspect_party_colors <- function(nodes, party_colors) {
  detected <- unique(nodes$party)
  cat("Detected parties (nodes$party):\n")
  print(detected)
  cat("\nNames in provided party_colors (if present):\n")
  if (exists("party_colors") && is.character(party_colors) && !is.null(names(party_colors))) {
    print(names(party_colors))
  } else {
    cat("(no party_colors found or not named)\n")
  }
  invisible(detected)
}

# Build a party_colors vector that exactly matches the detected parties (names = party)
# - party_colors_in: user-supplied named vector (may have different names/format)
# - nodes: tibble with nodes$party (cleaned party names)
# Returns: named character vector with names = detected parties
match_party_colors <- function(nodes, party_colors = NULL) {
  detected <- unique(nodes$party)
  # Normalize both sides to comparable keys (trim, collapse multiple spaces)
  norm <- function(x) trimws(gsub("\\s+", " ", x))
  detected_norm <- norm(detected)
  if (!is.null(party_colors) && is.character(party_colors) && !is.null(names(party_colors))) {
    user_names <- names(party_colors)
    user_norm <- norm(user_names)
    # attempt to match user-supplied names to detected parties
    map_idx <- match(detected_norm, user_norm)
    matched <- party_colors[map_idx]
    names(matched) <- detected
    # For unmatched (NA) assign defaults
    missing <- which(is.na(matched) | matched == "")
    if (length(missing) > 0) {
      # generate distinct colours for missing parties
      default_cols <- grDevices::rainbow(length(missing))
      matched[missing] <- default_cols
    }
  } else {
    # No user colours provided; create a palette
    cols <- grDevices::rainbow(length(detected))
    matched <- setNames(cols, detected)
  }
  return(matched)
}

# Rebuild the JS colourScale from a named vector of colours (names = parties)
build_colourScale_js <- function(party_colors_named) {
  party_names_js <- jsonlite::toJSON(names(party_colors_named))
  party_colors_js <- jsonlite::toJSON(unname(party_colors_named))
  sprintf("d3.scaleOrdinal().domain(%s).range(%s)", party_names_js, party_colors_js)
}

# Recreate sankey with corrected colours
# nodes_for_widget: a data.frame/tibble with columns full_name, party, label, name (name can be same as full_name)
# links_for_widget: data.frame with columns source (0-based), target (0-based), value
rebuild_sankey_with_colors <- function(nodes_for_widget, links_for_widget, party_colors_named) {
  colourScale <- build_colourScale_js(party_colors_named)
  sankey <- networkD3::sankeyNetwork(
    Links = links_for_widget,
    Nodes = nodes_for_widget,
    Source = "source", Target = "target",
    Value = "value", NodeID = "label", NodeGroup = "group",
    fontSize = 12, nodeWidth = 30,
    colourScale = htmlwidgets::JS(colourScale)
  )
  sankey
}

# Example flow (uncomment and run after the pipeline created nodes and links_df_filtered):
# detected <- inspect_party_colors(nodes, party_colors = if(exists("party_colors")) party_colors else NULL)
# party_colors_matched <- match_party_colors(nodes, if(exists("party_colors")) party_colors else NULL)
# # Recreate nodes_for_widget and links_for_widget to match expected columns:
# nodes_for_widget <- nodes %>% mutate(name = full_name)  # ensure 'name' present
# links_for_widget <- links_df_filtered %>% select(source = source_id, target = target_id, value)
# sankey_fixed <- rebuild_sankey_with_colors(nodes_for_widget, links_for_widget, party_colors_matched)
# print(sankey_fixed)
#
# If you still see mismatched colours after this, paste:
# - head(nodes)
# - names(party_colors) and party_colors (print)
# and I will check exact string mismatches (e.g., trailing spaces, capitalization differences).