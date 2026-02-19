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
model_names <- c(
  "variability_results.52", "variability_results.53", "variability_results.54",
  "variability_results.55", "variability_results.56", "variability_results.57"
)
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

# Extract data for each transition
pair_data <- imap(models, function(var_result, model_name) {
  year_suffix <- sub(".*\\.", "", model_name)
  target_year <- as.integer(paste0("18", year_suffix))
  source_year <- target_year - 1

  cell_summary <- transition_matrix_from_flow(var_result$flow_intervals)
  rhat <- var_result$rhat

  list(summary = cell_summary, rhat = rhat, target_year = target_year, source_year = source_year)
})
names(pair_data) <- names(models)

# Create sankey nodes
all_nodes <- map_df(pair_data, function(pair) {
  tibble(
    node = c(
      paste(rownames(pair$summary), "in", pair$source_year, sep = " "),
      paste(colnames(pair$summary), "in", pair$target_year, sep = " ")
    )
  )
}) %>%
  mutate(node = gsub(" ", "_", node)) %>%
  distinct() %>%
  arrange()

nodes <- tibble(full_name = all_nodes$node) %>%
  mutate(
    party = clean_party_label(full_name),
    label = sub("_in_", " ", full_name),
    label = gsub("_", " ", label),
    full_name = gsub(" ", "_", full_name) # Replace spaces with underscores
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
row.names(shares) <- gsub("non", "Non-voting_", gsub("vote_", "", row.names(shares)))

# Create sankey links with vote share flows
links_df_raw <- map_df(names(pair_data), function(pair_name) {
  # Extract the year from the pair name (e.g., "variability_results.52" -> 52)
  year_suffix <- sub("variability_results\\.", "", pair_name)
  target_year <- as.integer(paste0("20", year_suffix))  # or 1800 + as.integer(year_suffix) depending on your century
  
  rhat_by_transition(pair_data[[pair_name]]$rhat) %>%
    mutate(
      source = gsub(" ", "_", paste0(from_label, "_in_", from_year)),
      source_id = match(source, nodes$full_name) - 1,
      target = gsub(" ", "_", paste0(to_label, "_in_", to_year)),
      target_id = match(target, nodes$full_name) - 1,
      source_vote_share = shares[source, "share"]
    ) %>%
    mutate(
      beta = map2_dbl(from_label, to_label, function(from, to) {
        pair_data[[pair_name]]$summary[as.character(from), as.character(to)]
      }),
      flow_mean = source_vote_share * beta
    ) %>%
    select(source_id, target_id, flow_mean, max_rhat)
})

# Filter out insignificant transitions
min_frac_display <- 0.003
links_df_filtered <- links_df_raw %>% filter(flow_mean >= min_frac_display)

# Prepare the nodes and links for use in the sankey widget
nodes_for_widget <- nodes %>%
  transmute(id = full_name, label = label, group = party, full_name = full_name) %>%
  mutate(
    share = shares[full_name, ],
    nodecolor = party_colors[group],
    xpos = as.integer(node_year(full_name)) - 1851
  )

links_for_widget <- links_df_filtered %>%
  transmute(
    source = source_id,
    target = target_id,
    value = as.numeric(flow_mean), # Use actual flow, not just probability
    max_rhat = as.numeric(max_rhat)
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
          rhat: x.links.max_rhat && x.links.max_rhat[i]
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
    var rhatFormat = ".2f";
    var fmtRhat = d3.format(rhatFormat);
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
      var rhat   = (linkObj && linkObj.rhat !== undefined) ? +linkObj.rhat : (d && d.rhat !== undefined ? +d.rhat : NaN);

      var srcName = (d && d.source && (d.source.name || d.source.label)) ? (d.source.name || d.source.label) : visibleName(linkObj && linkObj.source);
      var tgtName = (d && d.target && (d.target.name || d.target.label)) ? (d.target.name || d.target.label) : visibleName(linkObj && linkObj.target);

      var meanTxt = isFinite(mean) ? fmt(mean) : "-";
      var rhatTxt   = isFinite(rhat)   ? fmtRhat(rhat)   : "-";

      var title = (srcName || "") + " \u2192 " + (tgtName || "") + nl +
                  meanTxt + (unit ? (" " + unit) : "") + nl +
                  "Maximum R-hat: " + rhatTxt;

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
