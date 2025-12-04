# Fix & rebuild sankey with exact party -> color mapping
# Run after the pipeline that created `nodes` and `links_df_filtered`.
# Expects: nodes (tibble with columns full_name, party, label, group),
#          links_df_filtered (data.frame with source_id, target_id, value),
#          party_colors (named char vector, optional)

library(jsonlite)
library(htmlwidgets)
library(networkD3)
library(dplyr)

# Diagnostic
cat("Detected node parties:\n"); print(unique(nodes$party))
cat("\nNames in provided party_colors (if present):\n")
if (exists("party_colors") && is.character(party_colors) && !is.null(names(party_colors))) {
  print(names(party_colors))
} else {
  cat("(no party_colors found or not a named vector)\n")
}

# Normalizer: replace underscores with spaces, collapse multi-space, trim
norm <- function(x) trimws(gsub("\\s+", " ", gsub("_", " ", as.character(x))))

detected <- unique(nodes$party)
det_norm <- norm(detected)

# If no user party_colors, create a default palette
if (!(exists("party_colors") && is.character(party_colors) && !is.null(names(party_colors)))) {
  message("No named party_colors found; creating default palette.")
  party_colors <- setNames(grDevices::rainbow(length(det_norm)), det_norm)
}

user_names <- names(party_colors)
user_norm <- norm(user_names)

# Map user colors to detected parties using normalized names
matched_colors <- rep(NA_character_, length(det_norm))
names(matched_colors) <- detected
for (i in seq_along(det_norm)) {
  mi <- match(det_norm[i], user_norm)
  if (!is.na(mi)) {
    matched_colors[i] <- party_colors[mi]
  } else {
    matched_colors[i] <- NA_character_
  }
}

# Report unmatched
if (any(is.na(matched_colors))) {
  missing_parties <- detected[is.na(matched_colors)]
  message("No matching color for detected parties (will assign defaults): ", paste(missing_parties, collapse = ", "))
  # assign defaults for missing
  defaults <- grDevices::rainbow(sum(is.na(matched_colors)))
  matched_colors[is.na(matched_colors)] <- defaults
}

# Ensure names are exactly the detected party strings
names(matched_colors) <- detected

# Rebuild JS colourScale from matched_colors
party_names_js <- jsonlite::toJSON(names(matched_colors))
party_colors_js <- jsonlite::toJSON(unname(matched_colors))
colourScale <- sprintf("d3.scaleOrdinal().domain(%s).range(%s)", party_names_js, party_colors_js)

# Rebuild sankey widget
nodes_for_widget <- nodes %>% mutate(name = full_name)  # keep internal unique ids
links_for_widget <- links_df_filtered %>% transmute(source = source_id, target = target_id, value = value)

sankey_fixed <- sankeyNetwork(
  Links = links_for_widget,
  Nodes = nodes_for_widget,
  Source = "source", Target = "target",
  Value = "value", NodeID = "label", NodeGroup = "group",
  fontSize = 12, nodeWidth = 30,
  colourScale = htmlwidgets::JS(colourScale)
)

print(sankey_fixed)

# Print final mapping for verification
cat("\nFinal party -> color mapping used:\n")
print(matched_colors)