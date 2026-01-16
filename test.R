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

variability_results <- quantify_ei_variability(1851, 1852, n_runs = 2)
estimate <- choose_final_estimate(variability_results)
betas.52 <- estimate$estimates
table.52 <- construct_contingency(results.52, betas.52, 1851, 1852)
print(table.52)

map <- create_map(1851, map.52, (results.1851 %>% dplyr::filter(!stringr::str_detect(town, "^Weighted")))$Free_Soil)
