library(eiPack)
library(ineq)
library(knitr)
library(tidyr)

source("betas.R")
source("present.R")
source("variables.R")

# Load the saved voter-inference results and perform statistical analysis

load("results.Rda")
load("1855_inference.Rda")
load("ct_demographics.Rda")

# Calculate weighted means and standard deviations for results

whig.mu.51 <- weighted.mean(results.1851$Whig, results.1851$weight)
whig.sd.51 <- sqrt(sum(results.1851$weight * (results.1851$Whig - whig.mu.51)^2))
dem.mu.51 <- weighted.mean(results.1851$Democrat , results.1851$weight)
dem.sd.51 <- sqrt(sum(results.1851$weight * (results.1851$Democrat - dem.mu.51)^2))
fs.mu.51 <- weighted.mean(results.1851$Free_Soil , results.1851$weight)
fs.sd.51 <- sqrt(sum(results.1851$weight * (results.1851$Free_Soil - fs.mu.51)^2))
non.mu.51 <- weighted.mean(results.1851$Abstaining, results.1851$weight)
non.sd.51 <- sqrt(sum(results.1851$weight * (results.1851$Abstaining - non.mu.51)^2))

whig.mu.52 <- weighted.mean(results.1852$Whig, results.1852$weight)
whig.sd.52 <- sqrt(sum(results.1852$weight * (results.1852$Whig - whig.mu.52)^2))
dem.mu.52 <- weighted.mean(results.1852$Democrat , results.1852$weight)
dem.sd.52 <- sqrt(sum(results.1852$weight * (results.1852$Democrat - dem.mu.52)^2))
fs.mu.52 <- weighted.mean(results.1852$Free_Soil , results.1852$weight)
fs.sd.52 <- sqrt(sum(results.1852$weight * (results.1852$Free_Soil - fs.mu.52)^2))
non.mu.52 <- weighted.mean(results.1852$Abstaining, results.1852$weight)
non.sd.52 <- sqrt(sum(results.1852$weight * (results.1852$Abstaining - non.mu.52)^2))

whig.mu.53 <- weighted.mean(results.1853$Whig, results.1853$weight)
whig.sd.53 <- sqrt(sum(results.1853$weight * (results.1853$Whig - whig.mu.53)^2))
dem.mu.53 <- weighted.mean(results.1853$Democrat , results.1853$weight)
dem.sd.53 <- sqrt(sum(results.1853$weight * (results.1853$Democrat - dem.mu.53)^2))
fs.mu.53 <- weighted.mean(results.1853$Free_Soil , results.1853$weight)
fs.sd.53 <- sqrt(sum(results.1853$weight * (results.1853$Free_Soil - fs.mu.53)^2))
non.mu.53 <- weighted.mean(results.1853$Abstaining, results.1853$weight)
non.sd.53 <- sqrt(sum(results.1853$weight * (results.1853$Abstaining - non.mu.53)^2))

whig.mu.54 <- weighted.mean(results.1854$Whig, results.1854$weight)
whig.sd.54 <- sqrt(sum(results.1854$weight * (results.1854$Whig - whig.mu.54)^2))
dem.mu.54 <- weighted.mean(results.1854$Democrat, results.1854$weight)
dem.sd.54 <- sqrt(sum(results.1854$weight * (results.1854$Democrat - dem.mu.54)^2))
fs.mu.54 <- weighted.mean(results.1854$Free_Soil, results.1854$weight)
fs.sd.54 <- sqrt(sum(results.1854$weight * (results.1854$Free_Soil - fs.mu.54)^2))
temp.mu.54 <- weighted.mean(results.1854$Temperence, results.1854$weight)
temp.sd.54 <- sqrt(sum(results.1854$weight * (results.1854$Temperence - temp.mu.54)^2))
non.mu.54 <- weighted.mean(results.1854$Abstaining, results.1854$weight)
non.sd.54 <- sqrt(sum(results.1854$weight * (results.1854$Abstaining - non.mu.54)^2))

whig.mu.55 <- weighted.mean(results.1855$Whig, results.1855$weight)
whig.sd.55 <- sqrt(sum(results.1855$weight * (results.1855$Whig - whig.mu.55)^2))
dem.mu.55 <- weighted.mean(results.1855$Democrat, results.1855$weight)
dem.sd.55 <- sqrt(sum(results.1855$weight * (results.1855$Democrat - dem.mu.55)^2))
kn.mu.55 <- weighted.mean(results.1855$Know_Nothing, results.1855$weight)
kn.sd.55 <- sqrt(sum(results.1855$weight * (results.1855$Know_Nothing - kn.mu.55)^2))
non.mu.55 <- weighted.mean(results.1855$Abstaining, results.1855$weight)
non.sd.55 <- sqrt(sum(results.1855$weight * (results.1855$Abstaining - non.mu.55)^2))

whig.mu.56 <- weighted.mean(results.1856$Whig, results.1856$weight)
whig.sd.56 <- sqrt(sum(results.1856$weight * (results.1856$Whig - whig.mu.56)^2))
dem.mu.56 <- weighted.mean(results.1856$Democrat, results.1856$weight)
dem.sd.56 <- sqrt(sum(results.1856$weight * (results.1856$Democrat - dem.mu.56)^2))
kn.mu.56 <- weighted.mean(results.1856$Know_Nothing, results.1856$weight)
kn.sd.56 <- sqrt(sum(results.1856$weight * (results.1856$Know_Nothing - kn.mu.56)^2))
rep.mu.56 <- weighted.mean(results.1856$Republican, results.1856$weight)
rep.sd.56 <- sqrt(sum(results.1856$weight * (results.1856$Republican - rep.mu.56)^2))
non.mu.56 <- weighted.mean(results.1856$Abstaining, results.1856$weight)
non.sd.56 <- sqrt(sum(results.1856$weight * (results.1856$Abstaining - non.mu.56)^2))

dem.mu.57 <- weighted.mean(results.1857$Democrat, results.1857$weight)
dem.sd.57 <- sqrt(sum(results.1857$weight * (results.1857$Democrat - dem.mu.57)^2))
rep.mu.57 <- weighted.mean(results.1857$Republican, results.1857$weight)
rep.sd.57 <- sqrt(sum(results.1857$weight * (results.1857$Republican - rep.mu.57)^2))
non.mu.57 <- weighted.mean(results.1857$Abstaining, results.1857$weight)
non.sd.57 <- sqrt(sum(results.1857$weight * (results.1857$Abstaining - non.mu.57)^2))

z_scores <- data.frame(matrix(ncol = 7, nrow = 7))
rownames(z_scores) <- (1851:1857)
colnames(z_scores) <- c("Democrat", "Whig", "Free Soil", "Temperence", "Know Nothing", "Republican", "Abstaining")
z_scores$Democrat[1] <- round((results.1851 %>% filter(town == "Meriden") %>% select(Democrat) - dem.mu.51) / dem.sd.51, 2)
z_scores$Democrat[2] <- round((results.1852 %>% filter(town == "Meriden") %>% select(Democrat) - dem.mu.52) / dem.sd.52, 2)
z_scores$Democrat[3] <- round((results.1853 %>% filter(town == "Meriden") %>% select(Democrat) - dem.mu.53) / dem.sd.53, 2)
z_scores$Democrat[4] <- round((results.1854 %>% filter(town == "Meriden") %>% select(Democrat) - dem.mu.54) / dem.sd.54, 2)
z_scores$Democrat[5] <- round((results.1855 %>% filter(town == "Meriden") %>% select(Democrat) - dem.mu.55) / dem.sd.55, 2)
z_scores$Democrat[6] <- round((results.1856 %>% filter(town == "Meriden") %>% select(Democrat) - dem.mu.56) / dem.sd.56, 2)
z_scores$Democrat[7] <- round((results.1857 %>% filter(town == "Meriden") %>% select(Democrat) - dem.mu.57) / dem.sd.57, 2)
z_scores$Abstaining[1] <- round((results.1851 %>% filter(town == "Meriden") %>% select(Abstaining) - non.mu.51) / non.sd.51, 2)
z_scores$Abstaining[2] <- round((results.1852 %>% filter(town == "Meriden") %>% select(Abstaining) - non.mu.52) / non.sd.52, 2)
z_scores$Abstaining[3] <- round((results.1853 %>% filter(town == "Meriden") %>% select(Abstaining) - non.mu.53) / non.sd.53, 2)
z_scores$Abstaining[4] <- round((results.1854 %>% filter(town == "Meriden") %>% select(Abstaining) - non.mu.54) / non.sd.54, 2)
z_scores$Abstaining[5] <- round((results.1855 %>% filter(town == "Meriden") %>% select(Abstaining) - non.mu.55) / non.sd.55, 2)
z_scores$Abstaining[6] <- round((results.1856 %>% filter(town == "Meriden") %>% select(Abstaining) - non.mu.56) / non.sd.56, 2)
z_scores$Abstaining[7] <- round((results.1857 %>% filter(town == "Meriden") %>% select(Abstaining) - non.mu.57) / non.sd.57, 2)
z_scores$Whig[1] <- round((results.1851 %>% filter(town == "Meriden") %>% select(Whig) - whig.mu.51) / whig.sd.51, 2)
z_scores$Whig[2] <- round((results.1852 %>% filter(town == "Meriden") %>% select(Whig) - whig.mu.52) / whig.sd.52, 2)
z_scores$Whig[3] <- round((results.1853 %>% filter(town == "Meriden") %>% select(Whig) - whig.mu.53) / whig.sd.53, 2)
z_scores$Whig[4] <- round((results.1854 %>% filter(town == "Meriden") %>% select(Whig) - whig.mu.54) / whig.sd.54, 2)
z_scores$Whig[5] <- round((results.1855 %>% filter(town == "Meriden") %>% select(Whig) - whig.mu.55) / whig.sd.55, 2)
z_scores$Whig[6] <- round((results.1856 %>% filter(town == "Meriden") %>% select(Whig) - whig.mu.56) / whig.sd.56, 2)
z_scores$"Free Soil"[1] <- round((results.1851 %>% filter(town == "Meriden") %>% select(Free_Soil) - fs.mu.51) / fs.sd.51, 2)
z_scores$"Free Soil"[2] <- round((results.1852 %>% filter(town == "Meriden") %>% select(Free_Soil) - fs.mu.52) / fs.sd.52, 2)
z_scores$"Free Soil"[3] <- round((results.1853 %>% filter(town == "Meriden") %>% select(Free_Soil) - fs.mu.53) / fs.sd.53, 2)
z_scores$"Free Soil"[4] <- round((results.1854 %>% filter(town == "Meriden") %>% select(Free_Soil) - fs.mu.54) / fs.sd.54, 2)
z_scores$Temperence[4] <- round((results.1854 %>% filter(town == "Meriden") %>% select(Temperence) - temp.mu.54) / temp.sd.54, 2)
z_scores$"Know Nothing"[5] <- round((results.1855 %>% filter(town == "Meriden") %>% select(Know_Nothing) - kn.mu.55) / kn.sd.55, 2)
z_scores$"Know Nothing"[6] <- round((results.1856 %>% filter(town == "Meriden") %>% select(Know_Nothing) - kn.mu.56) / kn.sd.56, 2)
z_scores$Republican[6] <- round((results.1856 %>% filter(town == "Meriden") %>% select(Republican) - rep.mu.56) / rep.sd.56, 2)
z_scores$Republican[7] <- round((results.1857 %>% filter(town == "Meriden") %>% select(Republican) - rep.mu.57) / rep.sd.57, 2)

print(kable(z_scores, caption = "Z-scores for Meriden results"))

# Meriden is the only town that passes a rough statistical proxy for
# Tyler Anbinder's argument that the Know Nothings arose because of
# anti-Nebraska sentiment.

# Look for all towns that voted within 0.35 standard deviations of the
# 1853 and 1854 mean Democratic and Whig shares but voted more than
# 0.6 standard deviations above the 1854 Free Soil share and more
# than 0.6 standard deviations above the 1855 mean Know Nothing share

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

# Meriden has a disproportionate number of native-born, young-adult males
# Stonington and New London, other Know Nothing hotbeds, show the same pattern
hist((ct_1850 %>% filter(BIRTH == "native" & SEX == 1))$AGE, main="Connecticut")
for (t in (distinct(ct_1850 %>% select(town))[[1]])) {
  native_male <- ct_1850 %>% filter(town == t & BIRTH == "native" & SEX == 1)
  if (nrow(native_male %>% filter(AGE >= 20 & AGE <= 25)) > nrow(native_male %>% filter(AGE >= 10 & AGE <= 15)) * 1.5) {
    hist(native_male$AGE, breaks = "Scott", main = t)
  }
}

# GINI index is strongly correlated to the degree of urbanization, using the
# percentage of household-head non-farm occupations as a proxy for urbanization,
# is inversely correlated to the average age of households heads, and is
# correlated to average household wealth.

summary(lm(gini ~ pct_farm_1860 + age_1860 + wealth, factors))

# The greatest uncertainty in the 1855 inference is in determining what became of
# 1854 Free Soil voters, followed by doubt about Temperance and Whig voters.
# Similar results are obtained when standard deviations are calculated for the
# models alphas.

combined_betas <- combine_betas(beta.sims.MD(ei.55, p55))

for (from_party in pull(distinct(combined_betas %>% select(from)))) {
  for (to_party in pull(distinct(combined_betas %>% select(to)))) {
    transition <- combined_betas %>%
      filter(from == from_party) %>%
      filter(to == to_party) %>%
      select(value)
    beta_sd <- sd(pull(transition))
    if (beta_sd > 2.5) {
      print(paste(from_party, " to ", to_party, ": ", as.character(beta_sd)))
    }
  }
}

test_cols <- paste(colnames(results.55 %>% select(LON, gini:pct_farm_1860)), collapse = " + ")

# A town's longitude is the factor that best explains support for the parties
# whose transitions are least certain.

test_parties <- c("Know_Nothing_in_1855", "Free_Soil_in_1854", "Temperance_in_1854", "Whig_in_1854")
for (party in test_parties) {
  form <- as.formula(paste(party, test_cols, sep = "~"))
  print(party)
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
kable(gini_grid, label = "1860 Meriden Wealth GINI Index by Age Range")

# Calculate real-property-only GINI indices for 1850 and 1860
gini_all <- data.frame(row.names = c("all"))
gini_by_birth <- data.frame(row.names = c("native", "immigrant"))
gini_by_job <- data.frame(row.names = c("farm", "nonfarm"))
for (age in age_cats) {
  g <- ct_1850_hh %>%
    filter(town == "Meriden") %>%
    filter(AGE_CAT == age) %>%
    summarise(gini = ineq(FAMILY_REALPROP, type = "Gini")) %>%
    select(gini)
  gini_all <- bind_cols(gini_all, g, .name_repair = "unique_quiet")
}
for (age in age_cats) {
  g <- ct_1850_hh %>%
    filter(town == "Meriden") %>%
    filter(AGE_CAT == age) %>%
    group_by(BIRTH) %>%
    summarise(gini = ineq(FAMILY_REALPROP, type = "Gini")) %>%
    select(gini)
  gini_by_birth <- bind_cols(gini_by_birth, g, .name_repair = "unique_quiet")
}
for (age in age_cats) {
  g <- ct_1850_hh %>%
    filter(town == "Meriden") %>%
    filter(AGE_CAT == age) %>%
    group_by(JOB) %>%
    summarise(gini = ineq(FAMILY_REALPROP, type = "Gini")) %>%
    select(gini)
  gini_by_job <- bind_cols(gini_by_job, g, .name_repair = "unique_quiet")
}
gini_grid <- bind_rows(gini_all, gini_by_birth, gini_by_job)
colnames(gini_grid) <- age_cats
kable(gini_grid, label = "1850 Meriden Real Property GINI Index by Age Range")

gini_all <- data.frame(row.names = c("all"))
gini_by_birth <- data.frame(row.names = c("native", "immigrant"))
gini_by_job <- data.frame(row.names = c("farm", "nonfarm"))
for (age in age_cats) {
  g <- ct_1860_hh %>%
    filter(town == "Meriden") %>%
    filter(AGE_CAT == age) %>%
    summarise(gini = ineq(FAMILY_REALPROP, type = "Gini")) %>%
    select(gini)
  gini_all <- bind_cols(gini_all, g, .name_repair = "unique_quiet")
}
for (age in age_cats) {
  g <- ct_1860_hh %>%
    filter(town == "Meriden") %>%
    filter(AGE_CAT == age) %>%
    group_by(BIRTH) %>%
    summarise(gini = ineq(FAMILY_REALPROP, type = "Gini")) %>%
    select(gini)
  gini_by_birth <- bind_cols(gini_by_birth, g, .name_repair = "unique_quiet")
}
for (age in age_cats) {
  g <- ct_1860_hh %>%
    filter(town == "Meriden") %>%
    filter(AGE_CAT == age) %>%
    group_by(JOB) %>%
    summarise(gini = ineq(FAMILY_REALPROP, type = "Gini")) %>%
    select(gini)
  gini_by_job <- bind_cols(gini_by_job, g, .name_repair = "unique_quiet")
}
gini_grid <- bind_rows(gini_all, gini_by_birth, gini_by_job)
colnames(gini_grid) <- age_cats
kable(gini_grid, label = "1860 Meriden Real Property GINI Index by Age Range")
