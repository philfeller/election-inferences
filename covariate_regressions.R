# Calculate ecological inferences of voters transitions using the multinomial
# Dirichlet method, with covariates.
#
# The eiPack pacakge is used to perform ecological inferences:
# https://www.rdocumentation.org/packages/eiPack/versions/0.1-7

source("./global.R")
source("./betas.R", local = TRUE)
source("./present.R", local = TRUE)
source("./ct_demographics.R", local = TRUE)
source("./prepare_results.R", local = TRUE)
source("./inference_utils.R", local = TRUE)

# Prepare a tibble with percentage election results for each pair of years. A
# separate tibble is created for each pair of years in order to have the most
# granular data possible, without unnecessarily combining towns that split before
# the period to be analyzed:

results <- create_results(raw_results, map.1855, 1851, 1852, eligible_pct, factors)
results.52 <- results$results
map.52 <- results$shp
results <- create_results(raw_results, map.1855, 1852, 1853, eligible_pct, factors)
results.53 <- results$results
map.53 <- results$shp
results <- create_results(raw_results, map.1855, 1853, 1854, eligible_pct, factors)
results.54 <- results$results
map.54 <- results$shp
results <- create_results(raw_results, map.1855, 1854, 1855, eligible_pct, factors)
results.55 <- results$results
map.55 <- results$shp
results <- create_results(raw_results, map.1857, 1855, 1856, eligible_pct, factors)
results.56 <- results$results
map.56 <- results$shp
results <- create_results(raw_results, map.1857, 1856, 1857, eligible_pct, factors)
results.57 <- results$results
map.57 <- results$shp
results <- create_results(raw_results, map.1857, 1851, 1857, eligible_pct, factors)
results.51_57 <- results$results
map.51_57 <- results$shp

save(
  raw_results, results.52, results.53, results.54, results.55, results.56, results.57, results.51_57,
  file = "results.Rda"
)

save(
  map.52, map.53, map.54, map.55, map.56, map.57, map.51_57,
  file = "maps.Rda"
)

# Create models for transitions
# Set parameters for gamma-distribution based on observed alpha results
# when ei.MD.bayes is run with default lambda1 and lambda2 values

variability_results.52 <- quantify_ei_variability(1851, 1852, n_runs = 10)
estimate <- choose_final_estimate(variability_results.52)
betas.52 <- estimate$estimates
table.52 <- construct_contingency(results.52, betas.52, 1851, 1852)
print(table.52)
save(variability_results.52, file = "1852_covariate_inference.Rda")

variability_results.53 <- quantify_ei_variability(1852, 1853, n_runs = 10)
estimate <- choose_final_estimate(variability_results.53)
betas.53 <- estimate$estimates
table.53 <- construct_contingency(results.53, betas.53, 1852, 1853)
print(table.53)
save(variability_results.53, file = "1853_covariate_inference.Rda")

variability_results.54 <- quantify_ei_variability(1853, 1854, n_runs = 10)
estimate <- choose_final_estimate(variability_results.54)
betas.54 <- estimate$estimates
table.54 <- construct_contingency(results.54, betas.54, 1853, 1854)
print(table.54)
save(variability_results.54, file = "1854_covariate_inference.Rda")

variability_results.55 <- quantify_ei_variability(1854, 1855, n_runs = 10)
estimate <- choose_final_estimate(variability_results.55)
betas.55 <- estimate$estimates
table.55 <- construct_contingency(results.55, betas.55, 1854, 1855)
print(table.55)
save(variability_results.55, file = "1855_covariate_inference.Rda")

variability_results.56 <- quantify_ei_variability(1855, 1856, n_runs = 10)
estimate <- choose_final_estimate(variability_results.56)
betas.56 <- estimate$estimates
table.56 <- construct_contingency(results.56, betas.56, 1855, 1856)
print(table.56)
save(variability_results.56, file = "1856_covariate_inference.Rda")

variability_results.57 <- quantify_ei_variability(1856, 1857, n_runs = 10)
estimate <- choose_final_estimate(variability_results.57)
betas.57 <- estimate$estimates
table.57 <- construct_contingency(results.57, betas.57, 1856, 1857)
print(table.57)
save(variability_results.57, file = "1857_covariate_inference.Rda")

variability_results.51-57 <- quantify_ei_variability(1851, 1857, n_runs = 10)
estimate <- choose_final_estimate(variability_results.51-57)
betas.51_57 <- estimate$estimates
table.51_57 <- construct_contingency(results.51_57, betas.51_57, 1851, 1857)
print(table.51_57)
save(variability_results.51-57, file = "1851_1857_covariate_inference.Rda")
