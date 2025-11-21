# Test the whether code works by running a simple regression and generating
# plots and maps

source("./global.R")
source("./betas.R", local = TRUE)
source("./present.R", local = TRUE)
source("./ct_demographics.R", local = TRUE)
source("./prepare_results.R", local = TRUE)
source("./inference_utils.R", local = TRUE)

data.52 <- create_results(raw_results, map.1855, 1851, 1852, eligible_pct, factors)
results.52 <- data.52$results
map.52 <- data.52$shp

ei.52 <- build_ei_model(1851, 1852)
table.52 <- construct_contingency(results.52, betas.MD(beta.sims.MD(ei.52, p52)), 1851, 1852)
print(table.52)

boxplot <- boxplotMD(beta.sims.MD(ei.52, p52), 1851, 1852)
ridge <- ridgelineMD(beta.sims.MD(ei.52, p52), 1851, 1852)
map <- create_map(1851, map.52, (results.1851 %>% dplyr::filter(!stringr::str_detect(town, "^Weighted")))$Free_Soil)
