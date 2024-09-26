library(eiPack)
source("betas.R")
source("present.R")
source("variables.R")

# Load the saved voter-inference results and preform statistical analysis

load("inferences.Rda")

# Meriden is the only town that passes a rough statistical proxy for
# Tyler Anbinder's argument that the Know Nothings arose because of
# anti-Nebraska sentiment.

# Look for all towns that voted within 0.55 standard deviations of the
# 1853 and 1854 mean Democratic and Whig shares but voted more than
# 0.55 standard deviations above the 1854 Free Soil share and more
# than 0.55 standard deviations above the 1855 mean Know Nothing share
weights.53 <- results.54$ELIG_1853 / sum(results.54$ELIG_1853)
whig.mu.53 <- weighted.mean(results.54$Whig_in_1853, weights.53)
whig.sd.53 <- sqrt(sum(weights.53 * (results.54$Whig_in_1853 - whig.mu.53)^2))
dem.mu.53 <- weighted.mean(results.54$Democrat_in_1853, weights.53)
dem.sd.53 <- sqrt(sum(weights.53 * (results.54$Democrat_in_1853 - dem.mu.53)^2))
weights.54 <- results.55$ELIG_1854 / sum(results.55$ELIG_1854)
whig.mu.54 <- weighted.mean(results.55$Whig_in_1854, weights.54)
whig.sd.54 <- sqrt(sum(weights.54 * (results.55$Whig_in_1854 - whig.mu.54)^2))
dem.mu.54 <- weighted.mean(results.55$Democrat_in_1854, weights.54)
dem.sd.54 <- sqrt(sum(weights.54 * (results.55$Democrat_in_1854 - dem.mu.54)^2))
fs.mu.54 <- weighted.mean(results.55$Free_Soil_in_1854, weights.54)
fs.sd.54 <- sqrt(sum(weights.54 * (results.55$Free_Soil_in_1854 - fs.mu.54)^2))
weights.55 <- results.55$ELIG_1855 / sum(results.55$ELIG_1855)
kn.mu.55 <- weighted.mean(results.55$Know_Nothing_in_1855, weights.55)
kn.sd.55 <- sqrt(sum(weights.55 * (results.55$Know_Nothing_in_1855 - kn.mu.55)^2))
non.mu.55 <- weighted.mean(results.55$Abstaining_in_1855, weights.55)
non.sd.55 <- sqrt(sum(weights.55 * (results.55$Abstaining_in_1855 - non.mu.55)^2))

for (t in 1:nrow(results.55)) {
  if (abs((results.55$Whig_in_1854[t] - whig.mu.54) / whig.sd.54) < 0.55 &&
    abs((results.55$Democrat_in_1854[t] - dem.mu.54) / dem.sd.54) < 0.55 &&
    abs((results.54$Whig_in_1853[t] - whig.mu.53) / whig.sd.53) < 0.55 &&
    abs((results.54$Democrat_in_1853[t] - dem.mu.53) / dem.sd.53) < 0.55 &&
    (results.55$Know_Nothing_in_1855[t] - kn.mu.55) / kn.sd.55 > 0.55 &&
    (results.55$Free_Soil_in_1854[t] - fs.mu.54) / fs.sd.54 > 0.55) {
    betas <- betas.MD(beta.sims.MD(ei.55, p55, t))
    print(kable(construct_contingency(results.55, betas, 1854, 1855, t, 1),
      caption = paste(t, results.55$town[t], sep = ": ")
    ))
  }
}

# Test whether differences in Meriden party shifts are statistically significant.
# For example, in 1855 both Democrats and Whigs were more likely to shift to
# Know Nothing than they were statewide, even at a 99.99% confidence level.

meriden <- results.55 %>% with(which(town == "Meriden"))
meriden.sims.55 <- beta.sims.MD(ei.55, p55, meriden)
ct.sims.55 <- beta.sims.MD(ei.55, p55)
t.test(meriden.sims.55[1, 3, ], ct.sims.55[1, 3, ], conf.level = 0.9999)
t.test(meriden.sims.55[2, 3, ], ct.sims.55[2, 3, ], conf.level = 0.9999)
