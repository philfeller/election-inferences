# Test the whether code works by running a simple regression and generating
# plots and maps

library(sf)
library(stringr)
source("betas.R")
source("present.R")
source("prepare_results.R")
source("inference_utils.R")

data.52 <- create_results(raw_results, map.1855, 1851, 1852, eligible_pct, factors)
results.52 <- data.52$results
map.52 <- data.52$shp

lambda2 <- 5 / 32.5
lambda1 <- 5 * lambda2
ei.52 <- build_ei_model(1851, 1852, lambda1, lambda2)
table.52 <- construct_contingency(results.52, betas.MD(beta.sims.MD(ei.52, p52)), 1851, 1852)
print(table.52)

boxplot <- boxplotMD(beta.sims.MD(ei.52, p52), 1851, 1852)
ridge <- ridgelineMD(beta.sims.MD(ei.52, p52), 1851, 1852)
map <- create_map(1851, map.52, (results.1851 %>% filter(!str_detect(town, "^Weighted")))$Free_Soil)
