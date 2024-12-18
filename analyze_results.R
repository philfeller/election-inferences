library(dplyr)
library(eiPack)
library(haven)
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

tag_yr <- function(results, yr, type) {
  results %>%
    filter(town == type) %>%
    mutate(weight = NULL) %>%
    rename(year = town) %>%
    mutate(year = yr)
}

result_mu <- results.1851 %>%
  tag_yr(1851, "Weighted mean") %>%
  bind_rows(results.1852 %>% tag_yr(1852, "Weighted mean")) %>%
  bind_rows(results.1853 %>% tag_yr(1853, "Weighted mean")) %>%
  bind_rows(results.1854 %>% tag_yr(1854, "Weighted mean")) %>%
  bind_rows(results.1855 %>% tag_yr(1855, "Weighted mean")) %>%
  bind_rows(results.1856 %>% tag_yr(1856, "Weighted mean")) %>%
  bind_rows(results.1857 %>% tag_yr(1857, "Weighted mean"))

result_sd <- results.1851 %>%
  tag_yr(1851, "Weighted SD") %>%
  bind_rows(results.1852 %>% tag_yr(1852, "Weighted SD")) %>%
  bind_rows(results.1853 %>% tag_yr(1853, "Weighted SD")) %>%
  bind_rows(results.1854 %>% tag_yr(1854, "Weighted SD")) %>%
  bind_rows(results.1855 %>% tag_yr(1855, "Weighted SD")) %>%
  bind_rows(results.1856 %>% tag_yr(1856, "Weighted SD")) %>%
  bind_rows(results.1857 %>% tag_yr(1857, "Weighted SD"))

single_stat <- function(stats, yr, party) {
  stats %>%
    filter(year == yr) %>%
    select(party)
}

# Calculate the Z-score for a town's party share in a given year
z_score <- function(results, mus, sds, yr, party, town) {
  town_result <- results %>%
    filter(town == "Meriden") %>%
    select(party)
  mu <- single_stat(mus, yr, party)
  sd <- single_stat(sds, yr, party)
  round((town_result - mu) / sd, 2)
}

# Calculate the Z-scores for Meriden from 1851 - 1857
z_scores <- data.frame(matrix(ncol = 7, nrow = 7))
years <- 1851:1857
parties <- c("Democrat", "Whig", "Free_Soil", "Temperence", "Know_Nothing", "Republican", "Abstaining")
rownames(z_scores) <- years
colnames(z_scores) <- parties
results_files <- list(results.1851, results.1852, results.1853, results.1854, results.1855, results.1856, results.1857)
for (party in parties) {
  for (yr in years) {
    i <- yr - 1850
    j <- match(party, parties)
    yr_results <- results_files[[i]]
    yr_parties <- colnames(yr_results)
    if (party %in% yr_parties) {
      z_scores[i, j] <- z_score(yr_results, result_mu, result_sd, yr, party, "Meriden")
    }
  }
}

print(kable(z_scores, caption = "Z-scores for Meriden results"))

# Meriden is the town that best demonstrates a rough statistical proxy
# for Tyler Anbinder's argument that the Know Nothings arose because of
# anti-Nebraska sentiment.

# Look for all towns that voted within 0.35 standard deviations of the
# 1853 and 1854 mean Democratic and Whig shares but voted more than
# 0.6 standard deviations above the 1854 Free Soil share and more
# than 0.6 standard deviations above the 1855 mean Know Nothing share.

whig.mu.53 <- result_mu %>% single_stat(1853, "Whig")
whig.sd.53 <- result_sd %>% single_stat(1853, "Whig")
dem.mu.53 <- result_mu %>% single_stat(1853, "Democrat")
dem.sd.53 <- result_sd %>% single_stat(1853, "Democrat")
whig.mu.54 <- result_mu %>% single_stat(1854, "Whig")
whig.sd.54 <- result_sd %>% single_stat(1854, "Whig")
dem.mu.54 <- result_mu %>% single_stat(1854, "Democrat")
dem.sd.54 <- result_sd %>% single_stat(1854, "Democrat")
fs.mu.54 <- result_mu %>% single_stat(1854, "Free_Soil")
fs.sd.54 <- result_sd %>% single_stat(1854, "Free_Soil")
kn.mu.55 <- result_mu %>% single_stat(1855, "Know_Nothing")
kn.sd.55 <- result_sd %>% single_stat(1855, "Know_Nothing")

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

# Meriden has a disproportionate number of native-born, young-adult males.
# Stonington and New London, other Know Nothing hotbeds, show the same pattern,
# as do the cities of Hartford and New Haven.

# Define a function to identify maxima
find_peaks <- function(x, m = 10) {
  shape <- diff(sign(diff(x, na.pad = FALSE)))
  pks <- sapply(which(shape < 0), FUN = function(i) {
    z <- i - m + 1
    z <- ifelse(z > 0, z, 1)
    w <- i + m + 1
    w <- ifelse(w < length(x), w, length(x))
    if (all(x[c(z:i, (i + 2):w)] <= x[i + 1])) {
      return(i + 1)
    } else {
      return(numeric(0))
    }
  })
  pks <- unlist(pks)
  pks
}

ct_native_male_1850 <- ct_1850 %>%
  filter(BIRTH == "native" & SEX == 1)
ct_native_male_1860 <- ct_1860 %>%
  filter(BIRTH == "native" & SEX == 1)
ct_age_distribution <- ct_native_male_1850 %>%
  select(AGE) %>%
  arrange(AGE) %>%
  group_by(AGE) %>%
  summarize(n = n())
ct_native_male.ss <- smooth.spline(ct_age_distribution)
plot(ct_age_distribution)
lines(ct_native_male.ss)
title("Connecticut Native-born Male Age Distribution")
for (t in (distinct(ct_1850 %>% select(town))[[1]])) {
  native_male_1850 <- ct_native_male_1850 %>%
    filter(town == t) %>%
    select(AGE) %>%
    arrange(AGE) %>%
    group_by(AGE) %>%
    summarize(n = n())
  native_male_1850.ss <- smooth.spline(native_male_1850)
  native_male_1860 <- ct_native_male_1860 %>%
    filter(town == t) %>%
    select(AGE) %>%
    arrange(AGE) %>%
    group_by(AGE) %>%
    summarize(n = n())
  native_male_1860.ss <- smooth.spline(native_male_1860)
  # Find the maxima in the smoothed age-distribution data
  peaks <- find_peaks(native_male_1850.ss$y, m = 5)
  if (length(peaks) < 3) {
    for (peak in peaks) {
      # Plot age distributions for towns that have a well-defined peak
      if (native_male_1850.ss$y[peak] - native_male_1850.ss$y[10] > 10) {
        max_age <- max(max(native_male_1850$n), max(native_male_1850.ss$y), max(native_male_1860.ss$y))
        plot(native_male_1850, ylim = c(0, max_age))
        lines(native_male_1850.ss)
        lines(native_male_1860.ss, lty = "dashed")
        legend("topright", legend = c("1850 ages", "1850 smoothed ages", "1860 smoothed ages"), pch = c(1, NA, NA), lty = c(NA, 1, 2))
        title(paste(t, "Native-born Male Age Distribution"))
        break
      }
    }
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
# whose transitions are least certain. Among the covariates that explained
# individual parties, one that is of interest is that the percentage of 
# low-wealth, young-adult, native-born males in 1850 negatively correlated to
# 1854 Whig support.

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
