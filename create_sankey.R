# Create a sankey chart showing the inferred party transitions
# Requires two packages that aren't installed in the Docker image,
# including sankey3Dplus: https://schmidtpaul.github.io/sankeyD3plus/index.html

# Install sankey3Dplus from GitHub
# remotes::install_github('SchmidtPaul/sankeyD3plus')

library(coda)
library(htmltools)
library(jsonlite)
library(sankeyD3plus)

source("./global.R")
source("./betas.R", local = TRUE)
source("./present.R", local = TRUE)
source("./results_utils.R", local = TRUE)
source("./inference_utils.R", local = TRUE)
source("./census_utils.R", local = TRUE)

# Load inferences with covariates
load("1857_covariate_inference.Rda")
load("1856_covariate_inference.Rda")
load("1855_covariate_inference.Rda")
load("1854_covariate_inference.Rda")
load("1853_covariate_inference.Rda")
load("1852_covariate_inference.Rda")
load("results.Rda")

# Collect all variability results into a list
model_names <- c("variability_results.52", "variability_results.53", "variability_results.54", 
                 "variability_results.55", "variability_results.56", "variability_results.57")
models <- mget(model_names, ifnotfound = NA, envir = .GlobalEnv)
if (any(is.na(models))) stop("One or more model names not found in global env. Check model_names.")

# Helper function to extract party name (strip '_in_185x' suffix and convert underscores to spaces)
clean_party_label <- function(full_label) {
  lbl <- gsub("_", " ", full_label)
  lbl <- sub(" in \\d{4}$", "", lbl)
  trimws(lbl)
}

# Helper function to extract election year
node_year <- function(full_label) {
  in_year <- sub("^.*_in_", "", full_label)
  trimws(in_year)
}

choose_method_for_year <- function(var_result) {
  max_cv <- max(var_result$variability$meta_cv)
  ci_overlap <- 100 - var_result$variability$pct_ci_overlap  # Convert to disagreement
  
  if (max_cv < 0.10 || ci_overlap < 10) {
    return("representative")
  } else if (max_cv < 0.25 && ci_overlap < 25) {
    return("consensus")
  } else {
    return("median")
  }
}

# Helper function to select appropriate estimate from variability results
# Options: "representative", "consensus", "median"
select_estimate <- function(var_result, method = "representative") {
  if (method == "representative") {
    # Use the representative run (closest to meta_mean)
    rep_run <- select_representative_run(var_result)
    return(list(
      mean = rep_run$results,
      ci_lower = rep_run$ci_lower,
      ci_upper = rep_run$ci_upper
    ))
  } else if (method == "consensus") {
    # Use weighted consensus
    consensus <- create_consensus_estimate(var_result)
    return(list(
      mean = consensus$consensus_mean,
      ci_lower = consensus$consensus_q025,
      ci_upper = consensus$consensus_q975
    ))
  } else if (method == "median") {
    # Use element-wise median
    median_est <- create_median_estimate(var_result)
    return(list(
      mean = median_est$median_estimate,
      ci_lower = median_est$q25,
      ci_upper = median_est$q75
    ))
  } else {
    stop("method must be 'representative', 'consensus', or 'median'")
  }
}

# Extract cell summaries from variability results
extract_cell_summary <- function(var_result, method = "representative") {
  estimates <- select_estimate(var_result, method)
  
  source_parties <- rownames(estimates$mean)
  target_parties <- colnames(estimates$mean)
  
  cell_summary <- tibble()
  
  for (i in seq_along(source_parties)) {
    for (j in seq_along(target_parties)) {
      source <- source_parties[i]
      target <- target_parties[j]
      
      cell_summary <- bind_rows(cell_summary, tibble(
        source = source,
        target = target,
        mean_share = estimates$mean[i, j],
        share_lo = estimates$ci_lower[i, j],
        share_hi = estimates$ci_upper[i, j],
      ))
    }
  }
  
  return(cell_summary)
}

# Extract data for each transition
pair_data <- imap(models, function(var_result, model_name) {
  method <- choose_method_for_year(var_result)
  
  cell_summary <- extract_cell_summary(var_result, method = method)
  
  list(summary = cell_summary, method_used = method)
})
names(pair_data) <- names(models)

pair_data <- imap(pair_data, function(data, model_name) {
  # Extract target year from name (e.g., "variability_results.52" -> "52" -> 1852)
  year_suffix <- sub(".*\\.", "", model_name)  # Gets "52"
  target_year <- as.integer(paste0("18", year_suffix))
  source_year <- target_year - 1
  
  cat("Processing", model_name, ": source =", source_year, ", target =", target_year, "\n")
  
  # Add year suffixes to party names
  data$summary <- data$summary %>%
    mutate(
      source = gsub(" ", "_", source),  # Replace spaces with underscores
      source = paste0(source, "_in_", source_year),
      target = gsub(" ", "_", target),  # Replace spaces with underscores
      target = paste0(target, "_in_", target_year),
      yr = as.character(source_year)  # Add yr column for compatibility
    )
  
  # Add year metadata to the list
  data$source_year <- source_year
  data$target_year <- target_year
  
  return(data)
})

# Create sankey nodes
all_nodes <- map(pair_data, ~ .$summary %>% select(source, target)) %>%
  bind_rows() %>%
  pivot_longer(everything(), values_to = "node") %>%
  mutate(node = gsub(" ", "_", node)) %>%
  distinct(node) %>%
  arrange(node)

nodes <- tibble(full_name = all_nodes$node) %>%
  mutate(
    party = clean_party_label(full_name),
    label = sub("_in_", " ", full_name),
    label = gsub("_", " ", label),
    full_name = gsub(" ", "_", full_name)  # Replace spaces with underscores
  )

# Node values based on actual vote shares
shares <- data.frame(get_shares(results.52, 1851), fix.empty.names = FALSE) %>% 
  rbind(data.frame(get_shares(results.52, 1852), fix.empty.names = FALSE)) %>%
  rbind(data.frame(get_shares(results.53, 1853), fix.empty.names = FALSE)) %>%
  rbind(data.frame(get_shares(results.54, 1854), fix.empty.names = FALSE)) %>%
  rbind(data.frame(get_shares(results.55, 1855), fix.empty.names = FALSE)) %>%
  rbind(data.frame(get_shares(results.56, 1856), fix.empty.names = FALSE)) %>%
  rbind(data.frame(get_shares(results.57, 1857), fix.empty.names = FALSE))
names(shares) <- c("share")
row.names(shares) <- gsub("non", "Non-voting_", gsub("vote_","", row.names(shares)))

# Create sankey links with vote share flows
links_df_raw <- map_dfr(pair_data, ~ .$summary %>% 
                          select(source, target, mean_share, share_lo, share_hi)) %>%
  mutate(
    source_id = match(source, nodes$full_name) - 1,
    target_id = match(target, nodes$full_name) - 1,
    # Get the source party's vote share
    source_vote_share = shares[source, "share"]
  ) %>%
  # Multiply transition probabilities by source vote share to get actual flows
  mutate(
    flow_mean = mean_share * source_vote_share,
    flow_lo = share_lo * source_vote_share,
    flow_hi = share_hi * source_vote_share
  )

# Filter out insignificant transitions
min_frac_display <- 0.004
links_df_filtered <- links_df_raw %>% filter(flow_mean >= min_frac_display)

# Prepare the nodes and links for use in the sankey widget
nodes_for_widget <- nodes %>%
  transmute(id = full_name, label = label, group = party, full_name = full_name) %>%
  mutate(
    share = shares[full_name,], 
    nodecolor = party_colors[group], 
    xpos = as.integer(node_year(full_name)) - 1851
  )

links_for_widget <- links_df_filtered %>%
  transmute(
    source = source_id,
    target = target_id,
    value = as.numeric(flow_mean),      # Use actual flow, not just probability
    share_lo = as.numeric(flow_lo),     # Use actual flow CI
    share_hi = as.numeric(flow_hi)      # Use actual flow CI
  )

links_df <- as.data.frame(links_for_widget, stringsAsFactors = FALSE) %>%
  mutate(link_id = paste0("L", seq_len(n())))

# Create sankey widget
sankey_widget <- sankeyNetwork(
  Links = links_df,
  Nodes = nodes_for_widget,
  NodePosX = "xpos", dragY = FALSE,
  NodeValue = "share", Source = "source", Target = "target",
  Value = "value", NodeID = "label", NodeGroup = "group",
  NodeColor = "nodecolor", showNodeValues = FALSE, align = "none",
  units = "vote share", numberFormat = ".1%",
  fontSize = 12, nodeWidth = 30, height = 500, width = 1100,
  zoom = TRUE
)

# Add onRender JavaScript (same as before)
js <- '
function(el, x) {
  try {
    // Build linksArr robustly (prefer column-wise x.links, then array x.links)
    var linksArr = [];
    if (x && x.links && Array.isArray(x.links.source)) {
      var n = x.links.source.length;
      for (var i = 0; i < n; i++) {
        linksArr[i] = {
          link_id: x.links.link_id && x.links.link_id[i],
          source: x.links.source[i],
          target: x.links.target && x.links.target[i],
          value:  x.links.value  && x.links.value[i],
          share_lo: x.links.share_lo && x.links.share_lo[i],
          share_hi: x.links.share_hi && x.links.share_hi[i]
        };
      }
    } else if (x && Array.isArray(x.links)) {
      linksArr = x.links;
    }

    // Build map by link_id for direct lookup
    var mapById = {};
    for (var k = 0; k < linksArr.length; k++) {
      var L = linksArr[k];
      if (L && L.link_id !== undefined && L.link_id !== null) {
        mapById[String(L.link_id)] = L;
      }
    }

    // Formatting & unit
    var numberFormat = (x && x.options && x.options.numberFormat) ? x.options.numberFormat : ".1%";
    var fmt = d3.format(numberFormat);
    var unit = (x && x.options && (x.options.units || x.options.unit)) ? (x.options.units || x.options.unit) : "";

    var nl = String.fromCharCode(10);

    // Helper: visible name resolution
    function visibleName(nodeRef) {
      if (!x || !x.nodes) return "";
      if (nodeRef && typeof nodeRef === "object") return nodeRef.name || nodeRef.label || "";
      var idx = +nodeRef;
      if (!isNaN(idx) && x.nodes[idx]) return x.nodes[idx].label || x.nodes[idx].name || "";
      return "";
    }

    // 1) Update tooltips for all links
    d3.select(el).selectAll(".link").each(function(d, i) {
      var bound = d || {};
      var lid = bound.link_id || (bound && bound.__data__ && bound.__data__.link_id) || null;
      var linkObj = null;
      if (lid && mapById[String(lid)]) {
        linkObj = mapById[String(lid)];
      } else if (linksArr && linksArr[i]) {
        linkObj = linksArr[i];
      } else {
        linkObj = bound;
      }

      var mean = (linkObj && linkObj.value !== undefined) ? +linkObj.value : (d && d.value !== undefined ? +d.value : NaN);
      var lo   = (linkObj && linkObj.share_lo !== undefined) ? +linkObj.share_lo : (d && d.share_lo !== undefined ? +d.share_lo : NaN);
      var hi   = (linkObj && linkObj.share_hi !== undefined) ? +linkObj.share_hi : (d && d.share_hi !== undefined ? +d.share_hi : NaN);

      var srcName = (d && d.source && (d.source.name || d.source.label)) ? (d.source.name || d.source.label) : visibleName(linkObj && linkObj.source);
      var tgtName = (d && d.target && (d.target.name || d.target.label)) ? (d.target.name || d.target.label) : visibleName(linkObj && linkObj.target);

      var meanTxt = isFinite(mean) ? fmt(mean) : "-";
      var loTxt   = isFinite(lo)   ? fmt(lo)   : "-";
      var hiTxt   = isFinite(hi)   ? fmt(hi)   : "-";

      var title = (srcName || "") + " \u2192 " + (tgtName || "") + nl +
                  meanTxt + (unit ? (" " + unit) : "") + nl +
                  "95% credible interval: [" + loTxt + ", " + hiTxt + "]";

      d3.select(this).selectAll("title").remove();
      d3.select(this).append("title").text(title);
    });

    // 2) Node hover highlight integration
    function resolveNodeIndex(nodeDatum) {
      if (nodeDatum === null || nodeDatum === undefined) return null;
      if (typeof nodeDatum === "number") return +nodeDatum;
      if (nodeDatum && (nodeDatum.index !== undefined)) return +nodeDatum.index;
      if ((nodeDatum.name || nodeDatum.label) && x && Array.isArray(x.nodes)) {
        var key = nodeDatum.name || nodeDatum.label;
        for (var j = 0; j < x.nodes.length; j++) {
          if ((x.nodes[j].name && x.nodes[j].name === key) || (x.nodes[j].label && x.nodes[j].label === key)) return j;
        }
      }
      return null;
    }

    var linkSel = d3.select(el).selectAll(".link");
    var linkNodes = linkSel.nodes();
    var linksMeta = new Array(linkNodes.length);

    linkSel.each(function(d, i) {
      var nodeEl = this;
      var meta = { source: nodeEl.__data__.source.name, target: nodeEl.__data__.target.name, node: nodeEl };

      var computed = d3.select(nodeEl).style("stroke-opacity");
      meta.origOpacity = (computed === null || computed === "" ) ? null : computed;
      var width = d3.select(nodeEl).style("stroke-width");
      meta.origWidth = (width === null || width === "") ? null : width;

      linksMeta[i] = meta;
    });

    function setHighlight(txt) {
      linkSel.each(function(d,i) {
        var m = linksMeta[i];
        var el = d3.select(this);
        var isConnected = (m && (m.source === txt || m.target === txt));
        if (isConnected) {
          el.style("stroke-opacity", .8);
        } else {
          el.style("stroke-opacity", 0.2);
        }
      });
    }

    function clearHighlight() {
      linkSel.each(function(d,i) {
        var m = linksMeta[i];
        var el = d3.select(this);
        if (m && m.origOpacity !== null && m.origOpacity !== undefined) el.style("stroke-opacity", m.origOpacity);
        else el.style("stroke-opacity", null);
        if (m && m.origWidth !== null && m.origWidth !== undefined) el.style("stroke-width", m.origWidth);
        else el.style("stroke-width", null);
      });
    }

    var nodeSel = d3.select(el).selectAll(".node");
    nodeSel.on("mouseover.highlightLinks", function(d,i) {
      var txt = d3.select(this).select("text").text();
      setHighlight(txt);
    });

    nodeSel.on("mouseout.highlightLinks", function(d,i) {
      clearHighlight();
    });

    d3.select(el).on("mouseleave.highlightLinks", function() { clearHighlight(); });
  } catch (err) {
    console.error("sankey onRender (tooltip+highlight) error:", err && err.stack ? err.stack : err);
  }
}
'

# Replace links with full data and add JavaScript
sankey_widget$x$links <- jsonlite::fromJSON(jsonlite::toJSON(links_df, dataframe = "rows", digits = 12))
sankey_widget <- htmlwidgets::onRender(sankey_widget, htmlwidgets::JS(js))

# Save HTML widget
save(sankey_widget, file = "sankey.Rda")
