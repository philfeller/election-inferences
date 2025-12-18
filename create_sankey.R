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
source("./census_utils.R", local = TRUE)

# Load inferences with covariates
load("1857_covariate_inference.Rda")
load("1856_covariate_inference.Rda")
load("1855_covariate_inference.Rda")
load("1854_covariate_inference.Rda")
load("1853_covariate_inference.Rda")
load("1852_covariate_inference.Rda")
load("results.Rda")

# Collect all models into a list
model_names <- c("cov.ei.52", "cov.ei.53", "cov.ei.54", "cov.ei.55", "cov.ei.56", "cov.ei.57")
models <- mget(model_names, ifnotfound = NA, envir = .GlobalEnv)
if (any(is.na(models))) stop("One or more model names not found in global env. Check model_names.")

# Helper function to extract party name (strip '_in_185x' suffix and convert underscores to spaces)
clean_party_label <- function(full_label) {
  lbl <- sub("_in_\\d{4}$", "", full_label)
  lbl <- gsub("_", " ", lbl)
  trimws(lbl)
}

# Helper function to extract election year
node_year <- function(full_label) {
  in_year <- sub("^.*_in_", "", full_label)
  trimws(in_year)
}

# Helper function to extract per-draw aggregated cell counts and parse cell names
extract_cell_draws <- function(em) {
  cc <- as.matrix(em$draws$Cell.counts) # draws x n_cells
  cell_names <- colnames(cc)
  if (is.null(cell_names)) stop("No colnames on Cell.counts; cannot parse cells.")
  parts <- strsplit(cell_names, "\\.")
  source <- sapply(parts, `[`, 2)
  target <- sapply(parts, `[`, 3)
  df <- as.data.frame(cc)
  names(df) <- cell_names
  # Use Heidelberger and Welsh to determine the point at which all betas converge.
  min_conv <- max(coda::heidel.diag(eiPack::lambda.MD(em, unique(target)))[, 2])
  draws <- seq_len(nrow(cc) - min_conv + 1)
  long <- as_tibble(df) %>%
    slice(min_conv:nrow(cc)) %>%
    mutate(.draw = draws) %>%
    pivot_longer(-.draw, names_to = "cell", values_to = "count") %>%
    mutate(
      source = source[match(cell, cell_names)],
      target = target[match(cell, cell_names)]
    )
  list(
    long = long, cell_names = cell_names, source = source, target = target,
    draws = nrow(cc), n_cells = length(cell_names)
  )
}

pair_data <- map(models, function(em) {
  ex <- extract_cell_draws(em)
  long <- ex$long
  draws <- max(long$.draw)
  long <- long %>%
    mutate(yr = node_year(source)) %>%
    group_by(yr) %>%
    mutate(yr_count = sum(count) / draws) %>%
    ungroup() %>%
    group_by(yr, .draw, source, target) %>%
    mutate(share = sum(count) / yr_count) %>%
    ungroup() %>%
    group_by(.draw, source) %>%
    mutate(
      source_total = sum(count),
      prop = ifelse(source_total > 0, count / source_total, 0)
    ) %>%
    ungroup()

  cell_summary <- long %>%
    group_by(yr, source, target) %>%
    summarize(
      mean_share = mean(share),
      share_lo = quantile(share, 0.025),
      share_hi = quantile(share, 0.975),
      .groups = "drop"
    )
  list(long = long, summary = cell_summary)
})
names(pair_data) <- names(models)

# Create sankey nodes
all_nodes <- map(pair_data, ~ .$summary %>% select(source, target)) %>%
  bind_rows() %>%
  pivot_longer(everything(), values_to = "node") %>%
  distinct(node) %>%
  arrange(node)

nodes <- tibble(full_name = all_nodes$node) %>%
  mutate(
    party = clean_party_label(full_name),
    yr = node_year(full_name),
    label = paste(party, yr)
  )

# Create sankey links
links_df_raw <- map_dfr(pair_data, ~ .$summary %>% select(source, target, mean_share, share_lo, share_hi)) %>%
  mutate(
    source_id = match(source, nodes$full_name) - 1,
    target_id = match(target, nodes$full_name) - 1
  )

# Node values are calculated based on the maximum of the sum of incoming and outgoing links,
# but rounding can cause this to differ from the actual vote share in the later years.
shares <- data.frame(get_shares(results.52, 1851), fix.empty.names = FALSE) %>% 
  rbind(data.frame(get_shares(results.52, 1852), fix.empty.names = FALSE)) %>%
  rbind(data.frame(get_shares(results.53, 1853), fix.empty.names = FALSE)) %>%
  rbind(data.frame(get_shares(results.54, 1854), fix.empty.names = FALSE)) %>%
  rbind(data.frame(get_shares(results.55, 1855), fix.empty.names = FALSE)) %>%
  rbind(data.frame(get_shares(results.56, 1856), fix.empty.names = FALSE)) %>%
  rbind(data.frame(get_shares(results.57, 1857), fix.empty.names = FALSE))
names(shares) <- c("share")
row.names(shares) <- gsub("non", "Abstaining_", gsub("vote_","", row.names(shares)))

# Assigning node values directly also allows filtering out of less significant transitions
min_frac_display <- 0.004
links_df_filtered <- links_df_raw %>% filter(mean_share >= min_frac_display)

# Prepare the nodes and links for use in the sankey widget
nodes_for_widget <- nodes %>%
  transmute(id = full_name, label = label, group = party, full_name = full_name) %>%
  mutate(share = shares[full_name,], nodecolor = party_colors[group], xpos = as.integer(node_year(full_name)) - 1851)
links_for_widget <- links_df_filtered %>%
  transmute(source = source_id,
            target = target_id,
            value = as.numeric(mean_share),
            share_lo = as.numeric(share_lo),
            share_hi = as.numeric(share_hi))
links_df <- as.data.frame(links_for_widget, stringsAsFactors = FALSE) %>%
  mutate(link_id = paste0("L", seq_len(n())))
                                                                                                                           
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

# Add onRender JavaScript that replaces the default link tooltip with one
# that also provides information about the 95% confidence interval
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

    // 1) Update tooltips for all links (use index-aligned linksArr when available)
    d3.select(el).selectAll(".link").each(function(d, i) {
      var bound = d || {};
      // Prefer link_id lookup, else index-based mapping from linksArr
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

      // Names: prefer bound datum node names, else resolve from linkObj indices
      var srcName = (d && d.source && (d.source.name || d.source.label)) ? (d.source.name || d.source.label) : visibleName(linkObj && linkObj.source);
      var tgtName = (d && d.target && (d.target.name || d.target.label)) ? (d.target.name || d.target.label) : visibleName(linkObj && linkObj.target);

      var meanTxt = isFinite(mean) ? fmt(mean) : "-";
      var loTxt   = isFinite(lo)   ? fmt(lo)   : "-";
      var hiTxt   = isFinite(hi)   ? fmt(hi)   : "-";

      var title = (srcName || "") + " \u2192 " + (tgtName || "") + nl +
                  meanTxt + (unit ? (" " + unit) : "") + nl +
                  "95% confidence interval: [" + loTxt + ", " + hiTxt + "]";

      d3.select(this).selectAll("title").remove();
      d3.select(this).append("title").text(title);
    });

    // 2) Node hover highlight integration
    // Helper to resolve node index from various datum shapes
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

    // Precompute metadata for each link element to speed hover work
    linkSel.each(function(d, i) {
      var nodeEl = this;
      var meta = { source: nodeEl.__data__.source.name, target: nodeEl.__data__.target.name, node: nodeEl };

      // store original styles
      var computed = d3.select(nodeEl).style("stroke-opacity");
      meta.origOpacity = (computed === null || computed === "" ) ? null : computed;
      var width = d3.select(nodeEl).style("stroke-width");
      meta.origWidth = (width === null || width === "") ? null : width;

      linksMeta[i] = meta;
    });

    // Highlight/dim functions
    function setHighlight(txt) {
      linkSel.each(function(d,i) {
        var m = linksMeta[i];
        var el = d3.select(this);
        var isConnected = (m && (m.source === txt || m.target === txt));
        if (isConnected) {
          // bring to front if possible and highlight
          // try { this.parentNode.appendChild(this); } catch(e) {}
          el.style("stroke-opacity", .8);
        } else {
          // dim non-connected links
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

    // Attach handlers to nodes (replace previous handlers in the same namespace)
    var nodeSel = d3.select(el).selectAll(".node");
    nodeSel.on("mouseover.highlightLinks", function(d,i) {
      var txt = d3.select(this).select("text").text();
      setHighlight(txt);
    });

    nodeSel.on("mouseout.highlightLinks", function(d,i) {
      clearHighlight();
    });

    // Also clear highlights if mouse leaves the widget area
    d3.select(el).on("mouseleave.highlightLinks", function() { clearHighlight(); });
  } catch (err) {
    console.error("sankey onRender (tooltip+highlight) error:", err && err.stack ? err.stack : err);
  }
}
'

# The links object in the sankey widget will only have the minimal data; need
# to replace it with the full data

sankey_widget$x$links <- jsonlite::fromJSON(jsonlite::toJSON(links_df, dataframe = "rows", digits = 12))
sankey_widget <- htmlwidgets::onRender(sankey_widget, htmlwidgets::JS(js))

# page <- tagList(
#   tags$h1("Sankey of party transitions (1851→1857)"),
#   tags$p("This chart shows inferred vote-share flows; link transparency indicates uncertainty (CI)."),
#   sankey_widget,
#   # tags$footer("Generated on ", Sys.Date())
# )

# Save HTML file
save(sankey_widget, file = "sankey.Rda")
rmarkdown::render("sankey.Rmd", output_file = "sankey.html", output_dir = "./html")
