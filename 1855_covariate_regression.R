# Calculate ecological inferences of voters transitions using the multinomial
# Dirichlet method, without covariates. If necessary, repeat with a greater
# burn-in factor until the results show convergence. Calculate the residuals
# and test them for normality. For any that aren't normally distributed, test
# the significance of various factors.
#
# The eiPack pacakge is used to perform ecological inferences:
# https://www.rdocumentation.org/packages/eiPack/versions/0.1-7

library(coda)
library(eiPack)
library(sf)
source("betas.R")
source("present.R")
source("prepare_results.R")

# Prepare a tibble with percentage election results for 1854 and 1855.

results.55 <- create_results(1854, 1855)

# Run an inference using longitude as a covariate

# Define model factors
tune_size <- 10000
sample <- 2000
thin <- 10
burnin <- 1000
covariate <- "~ LON"

lambda2 <- 2.5 / 16
lambda1 <- 2.5 * lambda2

tune.55 <- tuneMD(
  cbind(Democrat_in_1855, Whig_in_1855, Know_Nothing_in_1855, Abstaining_in_1855)
  ~ cbind(Democrat_in_1854, Whig_in_1854, Free_Soil_in_1854, Temperance_in_1854, Abstaining_in_1854),
  data = results.55, total = results.55$ELIG_1855, totaldraws = tune_size, ntunes = 10,
  lambda1 = lambda1, lambda2 = lambda2, covariate = covariate
)

h <- c(0, 1)
while (sum(h) != length(h)) {
  cov.ei.55 <- ei.MD.bayes(
    cbind(Democrat_in_1855, Whig_in_1855, Know_Nothing_in_1855, Abstaining_in_1855)
    ~ cbind(Democrat_in_1854, Whig_in_1854, Free_Soil_in_1854, Temperance_in_1854, Abstaining_in_1854),
    data = results.55, total = results.55$ELIG_1855,
    sample = sample, burnin = burnin, thin = thin,
    tune.list = tune.55, lambda1 = lambda1, lambda2 = lambda2, covariate = covariate
  )
  h <- heidel.diag(lambda.MD(cov.ei.55, p55))[, 1]
}
table.55 <- construct_contingency(results.55, betas.MD(beta.sims.MD(cov.ei.55, p55)), 1854, 1855)
print(table.55)

combined_betas <- combine_betas(beta.sims.MD(cov.ei.55, p55))

for (from_party in pull(distinct(combined_betas %>% select(from)))) {
  for (to_party in pull(distinct(combined_betas %>% select(to)))) {
    transition <- combined_betas %>% filter(from == from_party) %>% filter(to == to_party) %>% select(value)
    beta_sd <- sd(pull(transition))
    if (beta_sd > 2.5) {
      print(paste(from_party, " to ", to_party, ": ", as.character(beta_sd)))
    }
  }
}

save(cov.ei.55, file = "1855_covariate_inference.Rda")
