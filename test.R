# Test the whether code works by running a simple regression and generating
# plots and maps

library(coda)
library(eiPack)
library(sf)
source("betas.R")
source("present.R")
source("prepare_results.R")

results.52 <- create_results(1851,1852)

# Define model factors
sample <- 1000
thin <- 10
burnin <- 1000
lambda1 <- 7
lambda2 <- 1

tune.52 <- tuneMD(
  cbind(Democrat_in_1852, Whig_in_1852, Free_Soil_in_1852, Abstaining_in_1852)
  ~ cbind(Democrat_in_1851, Whig_in_1851, Free_Soil_in_1851, Abstaining_in_1851),
  data = results.52, total = results.52$ELIG_1852, totaldraws = 1000, ntunes = 10,
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
}
print(construct_contingency(results.52,betasMD(beta_simsMD(ei.52,p52)), 1851, 1852))
meriden <- results.52 %>% with(which(town=="Meriden"))
print(construct_contingency(results.52,betasMD(beta_simsMD(ei.52,p52,meriden)), 1851, 1852, meriden))

boxplot <- boxplotMD(beta_simsMD(ei.52,p52), 1851, 1852)
ridge <- ridgelineMD(beta_simsMD(ei.52,p52), 1851, 1852)
map <- create_map(1851,results.52$Free_Soil_vote_in_1852/results.52$total_1852)
