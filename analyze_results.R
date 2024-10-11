library(eiPack)
library(ineq)

source("betas.R")
source("present.R")
source("variables.R")

# Load the saved voter-inference results and perform statistical analysis

load("results.Rda")
load("1855_inference.Rda")
load("ct_demographics.Rda")

# Meriden is the only town that passes a rough statistical proxy for
# Tyler Anbinder's argument that the Know Nothings arose because of
# anti-Nebraska sentiment.

# Look for all towns that voted within 0.35 standard deviations of the
# 1853 and 1854 mean Democratic and Whig shares but voted more than
# 0.6 standard deviations above the 1854 Free Soil share and more
# than 0.6 standard deviations above the 1855 mean Know Nothing share
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
  if (abs((results.55$Whig_in_1854[t] - whig.mu.54) / whig.sd.54) < 0.35 &&
    abs((results.55$Democrat_in_1854[t] - dem.mu.54) / dem.sd.54) < 0.35 &&
    abs((results.54$Whig_in_1853[t] - whig.mu.53) / whig.sd.53) < 0.35 &&
    abs((results.54$Democrat_in_1853[t] - dem.mu.53) / dem.sd.53) < 0.35 &&
    (results.55$Know_Nothing_in_1855[t] - kn.mu.55) / kn.sd.55 > 0.6 &&
    (results.55$Free_Soil_in_1854[t] - fs.mu.54) / fs.sd.54 > 0.6) {
    betas <- betas.MD(beta.sims.MD(ei.55, p55, t))
    print(kable(construct_contingency(results.55, betas, 1854, 1855, t, 1),
      caption = paste(t, results.55$town[t], sep = ": ")
    ))
  }
}

# GINI index is strongly correlated to the degree of urbanization, using the
# percentage of household-head non-farm occupations as a proxy for urbanization,
# is inversely correlated to the average age of households heads, and is
# correlated to average household wealth.

summary(lm(gini ~ pct_farm_1860 + age_1860 + wealth, factors))

# 1855 party support is best described by a town's longitude.

for (party in colnames(results.55)[grep("[aiocn][tgle]_in_1855", colnames(results.55))]) {
  form <- as.formula(paste(party, "LON", sep = "~"))
  print(summary(lm(form, results.55, weights = results.55$ELIG_1855)))
}

# Show GINI indices for Meriden subgroups

age_cats <- c("0 - 20", "20 - 30", "30 - 40", "40 - 50", "50 - 60", "Over 60")
gini_all <- data.frame(row.names = c("all"))
gini_by_birth <- data.frame(row.names = c("native", "immigrant"))
gini_by_job <- data.frame(row.names = c("farm", "nonfarm"))
for (age in age_cats) {
  g <- ct_1860_hh %>%
    filter(town == "Meriden") %>%
    filter(AGE_CAT == age) %>%
    summarise(gini = ineq(FAMILY_WEALTH, type = "Gini")) %>%
    select(gini)
  gini_all <- bind_cols(gini_all, g, .name_repair = "unique_quiet")
}
for (age in age_cats) {
  g <- ct_1860_hh %>%
    filter(town == "Meriden") %>%
    filter(AGE_CAT == age) %>%
    group_by(BIRTH) %>%
    summarise(gini = ineq(FAMILY_WEALTH, type = "Gini")) %>%
    select(gini)
  gini_by_birth <- bind_cols(gini_by_birth, g, .name_repair = "unique_quiet")
}
for (age in age_cats) {
  g <- ct_1860_hh %>%
    filter(town == "Meriden") %>%
    filter(AGE_CAT == age) %>%
    group_by(JOB) %>%
    summarise(gini = ineq(FAMILY_WEALTH, type = "Gini")) %>%
    select(gini)
  gini_by_job <- bind_cols(gini_by_job, g, .name_repair = "unique_quiet")
}
gini_grid <- bind_rows(gini_all, gini_by_birth, gini_by_job)
colnames(gini_grid) <- age_cats
kable(gini_grid, label = "1860 Meriden GINI Index by Age Range")
