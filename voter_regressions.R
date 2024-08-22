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

# Prepare a tibble with percentage election results for each pair of years. A
# separate tibble is created for each pair of years in order to have the most
# granular data possible, without unnecessarily combining towns that split before
# the period to be analyzed:

results.52 <- create_results(1851, 1852)
results.53 <- create_results(1852, 1853)
results.54 <- create_results(1853, 1854)
results.55 <- create_results(1854, 1855)
results.56 <- create_results(1855, 1856)
results.57 <- create_results(1856, 1857)
results.51_57 <- create_results(1851, 1857)

# Although the model factors for each pair of years have been chosen to produce MCMC results
# that converge fairly reliably, individual runs sometimes fail to converge. The calls to
# the ei.MD.bayes function are wrapped in a loop that tests for convergence to ensure that
# results are useable.

tune_size <- 10000
# Define model factors
sample <- 1000
thin <- 10
burnin <- 1000
lambda1 <- 7
lambda2 <- 1

tune.52 <- tuneMD(
  cbind(Democrat_in_1852, Whig_in_1852, Free_Soil_in_1852, Abstaining_in_1852)
  ~ cbind(Democrat_in_1851, Whig_in_1851, Free_Soil_in_1851, Abstaining_in_1851),
  data = results.52, total = results.52$ELIG_1852, totaldraws = tune_size, ntunes = 10,
  lambda1 = lambda1, lambda2 = lambda2
)

h <- c(0, 1)
while (sum(h) != length(h)) {
  ei.52 <- ei.MD.bayes(
    cbind(Democrat_in_1852, Whig_in_1852, Free_Soil_in_1852, Abstaining_in_1852)
    ~ cbind(Democrat_in_1851, Whig_in_1851, Free_Soil_in_1851, Abstaining_in_1851),
    data = results.52, total = results.52$ELIG_1852,
    sample = sample, burnin = burnin, thin = thin,
    tune.list = tune.52, lambda1 = lambda1, lambda2 = lambda2
  )
  h <- heidel.diag(lambda.MD(ei.52, p52))[, 1]
}
if (sum(h) != length(h)) {
  print(paste("1852 doesn't converge"))
  print(h[h == 0])
} else {
  table.52 <- construct_contingency(results.52, betasMD(beta_simsMD(ei.52, p52)), 1851, 1852)
  n52 <- nrow(results.52)
  meriden <- results.52 %>% with(which(town == "Meriden"))
  meriden.52 <- construct_contingency(results.52, betasMD(beta_simsMD(ei.52, p52, meriden)), 1851, 1852, meriden)
  print(table.52)
}

sample <- 1000
thin <- 10000
burnin <- 1000
lambda1 <- 9
lambda2 <- 2

tune.53 <- tuneMD(
  cbind(Democrat_in_1853, Whig_in_1853, Free_Soil_in_1853, Abstaining_in_1853)
  ~ cbind(Democrat_in_1852, Whig_in_1852, Free_Soil_in_1852, Abstaining_in_1852),
  data = results.53, total = results.53$ELIG_1853, totaldraws = tune_size, ntunes = 10,
  lambda1 = lambda1, lambda2 = lambda2
)

h <- c(0, 1)
while (sum(h) != length(h)) {
  ei.53 <- ei.MD.bayes(
    cbind(Democrat_in_1853, Whig_in_1853, Free_Soil_in_1853, Abstaining_in_1853)
    ~ cbind(Democrat_in_1852, Whig_in_1852, Free_Soil_in_1852, Abstaining_in_1852),
    data = results.53, total = results.53$ELIG_1853,
    sample = sample, burnin = burnin, thin = thin,
    tune.list = tune.53, lambda1 = lambda1, lambda2 = lambda2
  )
  h <- heidel.diag(lambda.MD(ei.53, p53))[, 1]
}
if (sum(h) != length(h)) {
  print(paste("1853 doesn't converge"))
  print(h[h == 0])
} else {
  table.53 <- construct_contingency(results.53, betasMD(beta_simsMD(ei.53, p53)), 1852, 1853)
  n53 <- nrow(results.53)
  meriden <- results.53 %>% with(which(town == "Meriden"))
  meriden.53 <- construct_contingency(results.53, betasMD(beta_simsMD(ei.53, p53, meriden)), 1852, 1853, meriden)
  print(table.53)
}

sample <- 1000
thin <- 10000
burnin <- 10000
lambda1 <- 7
lambda2 <- 1

tune.54 <- tuneMD(
  cbind(Democrat_in_1854, Whig_in_1854, Free_Soil_in_1854, Temperance_in_1854, Abstaining_in_1854)
  ~ cbind(Democrat_in_1853, Whig_in_1853, Free_Soil_in_1853, Abstaining_in_1853),
  data = results.54, total = results.54$ELIG_1854, totaldraws = tune_size, ntunes = 10,
  lambda1 = lambda1, lambda2 = lambda2
)

h <- c(0, 1)
while (sum(h) != length(h)) {
  ei.54 <- ei.MD.bayes(
    cbind(Democrat_in_1854, Whig_in_1854, Free_Soil_in_1854, Temperance_in_1854, Abstaining_in_1854)
    ~ cbind(Democrat_in_1853, Whig_in_1853, Free_Soil_in_1853, Abstaining_in_1853),
    data = results.54, total = results.54$ELIG_1854,
    sample = sample, burnin = burnin, thin = thin,
    tune.list = tune.54, lambda1 = lambda1, lambda2 = lambda2
  )
  h <- heidel.diag(lambda.MD(ei.54, p54))[, 1]
}
if (sum(h) != length(h)) {
  print(paste("1854 doesn't converge"))
  print(h[h == 0])
} else {
  table.54 <- construct_contingency(results.54, betasMD(beta_simsMD(ei.54, p54)), 1853, 1854)
  n54 <- nrow(results.54)
  meriden <- results.54 %>% with(which(town == "Meriden"))
  meriden.54 <- construct_contingency(results.54, betasMD(beta_simsMD(ei.54, p54, meriden)), 1853, 1854, meriden)
  print(table.54)
}

sample <- 1000
thin <- 10000
burnin <- 20000
lambda1 <- 7
lambda2 <- 1

tune.55 <- tuneMD(
  cbind(Democrat_in_1855, Whig_in_1855, Know_Nothing_in_1855, Abstaining_in_1855)
  ~ cbind(Democrat_in_1854, Whig_in_1854, Free_Soil_in_1854, Temperance_in_1854, Abstaining_in_1854),
  data = results.55, total = results.55$ELIG_1855, totaldraws = tune_size, ntunes = 10,
  lambda1 = lambda1, lambda2 = lambda2
)

h <- c(0, 1)
while (sum(h) != length(h)) {
  ei.55 <- ei.MD.bayes(
    cbind(Democrat_in_1855, Whig_in_1855, Know_Nothing_in_1855, Abstaining_in_1855)
    ~ cbind(Democrat_in_1854, Whig_in_1854, Free_Soil_in_1854, Temperance_in_1854, Abstaining_in_1854),
    data = results.55, total = results.55$ELIG_1855,
    sample = sample, burnin = burnin, thin = thin,
    tune.list = tune.55, lambda1 = lambda1, lambda2 = lambda2
  )
  h <- heidel.diag(lambda.MD(ei.55, p55))[, 1]
}
if (sum(h) != length(h)) {
  print(paste("1855 doesn't converge"))
  print(h[h == 0])
} else {
  table.55 <- construct_contingency(results.55, betasMD(beta_simsMD(ei.55, p55)), 1854, 1855)
  n55 <- nrow(results.55)
  meriden <- results.55 %>% with(which(town == "Meriden"))
  meriden.55 <- construct_contingency(results.55, betasMD(beta_simsMD(ei.55, p55, meriden)), 1854, 1855, meriden)
  print(table.55)
}

tune.56 <- tuneMD(
  cbind(Democrat_in_1856, Whig_in_1856, Know_Nothing_in_1856, Republican_in_1856, Abstaining_in_1856)
  ~ cbind(Democrat_in_1855, Whig_in_1855, Know_Nothing_in_1855, Abstaining_in_1855),
  data = results.56, total = results.56$ELIG_1856, totaldraws = tune_size, ntunes = 10,
  lambda1 = lambda1, lambda2 = lambda2
)

h <- c(0, 1)
while (sum(h) != length(h)) {
  ei.56 <- ei.MD.bayes(
    cbind(Democrat_in_1856, Whig_in_1856, Know_Nothing_in_1856, Republican_in_1856, Abstaining_in_1856)
    ~ cbind(Democrat_in_1855, Whig_in_1855, Know_Nothing_in_1855, Abstaining_in_1855),
    data = results.56, total = results.56$ELIG_1856,
    sample = sample, burnin = burnin, thin = thin,
    tune.list = tune.56, lambda1 = lambda1, lambda2 = lambda2
  )
  h <- heidel.diag(lambda.MD(ei.56, p56))[, 1]
}
if (sum(h) != length(h)) {
  print(paste("1856 doesn't converge"))
  print(h[h == 0])
} else {
  table.56 <- construct_contingency(results.56, betasMD(beta_simsMD(ei.56, p56)), 1855, 1856)
  n56 <- nrow(results.56)
  meriden <- results.56 %>% with(which(town == "Meriden"))
  meriden.56 <- construct_contingency(results.56, betasMD(beta_simsMD(ei.56, p56, meriden)), 1855, 1856, meriden)
  print(table.56)
}

sample <- 1000
thin <- 100
burnin <- 20000
lambda1 <- 9
lambda2 <- 2

tune.57 <- tuneMD(
  cbind(Democrat_in_1857, Republican_in_1857, Abstaining_in_1857)
  ~ cbind(Democrat_in_1856, Whig_in_1856, Know_Nothing_in_1856, Republican_in_1856, Abstaining_in_1856),
  data = results.57, total = results.57$ELIG_1857, totaldraws = tune_size, ntunes = 10,
  lambda1 = lambda1, lambda2 = lambda2
)

h <- c(0, 1)
while (sum(h) != length(h)) {
  ei.57 <- ei.MD.bayes(
    cbind(Democrat_in_1857, Republican_in_1857, Abstaining_in_1857)
    ~ cbind(Democrat_in_1856, Whig_in_1856, Know_Nothing_in_1856, Republican_in_1856, Abstaining_in_1856),
    data = results.57, total = results.57$ELIG_1857,
    sample = sample, burnin = burnin, thin = thin,
    tune.list = tune.57, lambda1 = lambda1, lambda2 = lambda2
  )
  h <- heidel.diag(lambda.MD(ei.57, p57))[, 1]
}
if (sum(h) != length(h)) {
  print(paste("1857 doesn't converge"))
  print(h[h == 0])
} else {
  table.57 <- construct_contingency(results.57, betasMD(beta_simsMD(ei.57, p57)), 1856, 1857)
  n57 <- nrow(results.57)
  meriden <- results.57 %>% with(which(town == "Meriden"))
  meriden.57 <- construct_contingency(results.57, betasMD(beta_simsMD(ei.57, p57, meriden)), 1856, 1857, meriden)
  print(table.57)
}

sample <- 1000
thin <- 5000
burnin <- 10000
lambda1 <- 7
lambda2 <- 1

tune.51_57 <- tuneMD(
  cbind(Democrat_in_1857, Republican_in_1857, Abstaining_in_1857)
  ~ cbind(Democrat_in_1851, Whig_in_1851, Free_Soil_in_1851, Abstaining_in_1851),
  data = results.51_57, total = results.51_57$ELIG_1857, totaldraws = tune_size, ntunes = 10,
  lambda1 = lambda1, lambda2 = lambda2
)

h <- c(0, 1)
while (sum(h) != length(h)) {
  ei.51_57 <- ei.MD.bayes(
    cbind(Democrat_in_1857, Republican_in_1857, Abstaining_in_1857)
    ~ cbind(Democrat_in_1851, Whig_in_1851, Free_Soil_in_1851, Abstaining_in_1851),
    data = results.51_57, total = results.51_57$ELIG_1857,
    sample = sample, burnin = burnin, thin = thin,
    tune.list = tune.51_57, lambda1 = lambda1, lambda2 = lambda2
  )
  h <- heidel.diag(lambda.MD(ei.51_57, p57))[, 1]
}
if (sum(h) != length(h)) {
  print(paste("1851-57 doesn't converge"))
  print(h[h == 0])
} else {
  table.51_57 <- construct_contingency(results.51_57, betasMD(beta_simsMD(ei.51_57, p57)), 1851, 1857)
  n51_57 <- nrow(results.51_57)
  meriden <- results.51_57 %>% with(which(town == "Meriden"))
  meriden.51_57 <- construct_contingency(results.51_57, betasMD(beta_simsMD(ei.51_57, p57, meriden)), 1851, 1857, meriden)
  print(table.51_57)
}

save(
  results.52, results.53, results.54, results.55, results.56, results.51_57,
  ei.52, ei.53, ei.54, ei.55, ei.56, ei.51_57,
  file = "inferences.Rda"
)
