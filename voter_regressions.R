# Calculate ecological inferences of voters transitions using the multinomial
# Dirichlet method, with covariates.
#
# The eiPack pacakge is used to perform ecological inferences:
# https://www.rdocumentation.org/packages/eiPack/versions/0.1-7

library(sf)
source("betas.R")
source("present.R")
source("prepare_results.R")
source("inference_utils.R")

# Prepare a tibble with percentage election results for each pair of years. A
# separate tibble is created for each pair of years in order to have the most
# granular data possible, without unnecessarily combining towns that split before
# the period to be analyzed:

results.52 <- create_results(raw_results, map.1855, 1851, 1852)$results
results.53 <- create_results(raw_results, map.1855, 1852, 1853)$results
results.54 <- create_results(raw_results, map.1855, 1853, 1854)$results
results.55 <- create_results(raw_results, map.1855, 1854, 1855)$results
results.56 <- create_results(raw_results, map.1857, 1855, 1856)$results
results.57 <- create_results(raw_results, map.1857, 1856, 1857)$results
results.51_57 <- create_results(raw_results, map.1857, 1851, 1857)$results

save(
  raw_results, results.52, results.53, results.54, results.55, results.56, results.51_57,
  file = "results.Rda"
)

# Although the model factors for each pair of years have been chosen to produce MCMC results
# that converge fairly reliably, individual runs sometimes fail to converge. The calls to
# the ei.MD.bayes function are wrapped in a loop that tests for convergence to ensure that
# results are useable.

# Create models for transitions
# Set parameters for gamma-distribution based on observed alpha results
# when ei.MD.bayes is run with default lambda1 and lambda2 values

lambda2 <- 5 / 32.5
lambda1 <- 5 * lambda2
ei.52 <- build_ei_model(1851, 1852, lambda1, lambda2)
table.52 <- construct_contingency(results.52, betas.MD(beta.sims.MD(ei.52, p52)), 1851, 1852)
print(table.52)
save(ei.52, file = "1852_inference.Rda")

lambda2 <- 4.55 / 56
lambda1 <- 5 * lambda2
ei.53 <- build_ei_model(1852, 1853, lambda1, lambda2)
table.53 <- construct_contingency(results.53, betas.MD(beta.sims.MD(ei.53, p53)), 1852, 1853)
print(table.53)
save(ei.53, file = "1853_inference.Rda")

lambda2 <- 1.75 / 13.5
lambda1 <- 1.75 * lambda2
ei.54 <- build_ei_model(1853, 1854, lambda1, lambda2)
table.54 <- construct_contingency(results.54, betas.MD(beta.sims.MD(ei.54, p54)), 1853, 1854)
print(table.54)
save(ei.54, file = "1854_inference.Rda")

lambda2 <- 2.5 / 16
lambda1 <- 2.5 * lambda2
ei.55 <- build_ei_model(1854, 1855, lambda1, lambda2)
table.55 <- construct_contingency(results.55, betas.MD(beta.sims.MD(ei.55, p55)), 1854, 1855)
print(table.55)
save(ei.55, file = "1855_inference.Rda")

lambda2 <- 2 / 16
lambda1 <- 2 * lambda2
ei.56 <- build_ei_model(1855, 1856, lambda1, lambda2)
table.56 <- construct_contingency(results.56, betas.MD(beta.sims.MD(ei.56, p56)), 1855, 1856)
print(table.56)
save(ei.56, file = "1856_inference.Rda")

lambda2 <- 2 / 4
lambda1 <- 2 * lambda2
ei.57 <- build_ei_model(1856, 1857, lambda1, lambda2)
table.57 <- construct_contingency(results.57, betas.MD(beta.sims.MD(ei.57, p57)), 1856, 1857)
print(table.57)
save(ei.57, file = "1857_inference.Rda")

lambda2 <- 2 / 5.5
lambda1 <- 2 * lambda2
ei.51_57 <- build_ei_model(1851, 1857, lambda1, lambda2, results_tibble = results.51_57)
table.51_57 <- construct_contingency(results.51_57, betas.MD(beta.sims.MD(ei.51_57, p57)), 1851, 1857)
print(table.51_57)
save(ei.51_57, file = "1851_1857_inference.Rda")
