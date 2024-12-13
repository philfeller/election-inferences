# Test the whether code works by running a simple regression and generating
# plots and maps

library(coda)
library(eiPack)
library(sf)
library(stringr)
source("betas.R")
source("present.R")
source("prepare_results.R")

results.52 <- create_results(1851, 1852)

# Define model factors
tune_size <- 10000
sample <- 2000
thin <- 10
burnin <- 1000

# Set parameters for gamma-distribution based on observed alpha results
lambda2 <- 5 / 32.5
lambda1 <- 5 * lambda2

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

print(construct_contingency(results.52, betas.MD(beta.sims.MD(ei.52, p52)), 1851, 1852))
meriden <- results.52 %>% with(which(town == "Meriden"))
print(construct_contingency(results.52, betas.MD(beta.sims.MD(ei.52, p52, meriden)), 1851, 1852, meriden))

boxplot <- boxplotMD(beta.sims.MD(ei.52, p52), 1851, 1852)
ridge <- ridgelineMD(beta.sims.MD(ei.52, p52), 1851, 1852)
map <- create_map(1851, (results.1851 %>% filter(! str_detect(town, "^Weighted")))$Free_Soil)
