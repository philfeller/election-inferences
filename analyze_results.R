source("./global.R")
source("./betas.R", local = TRUE)
source("./present.R", local = TRUE)
source("./results_utils.R", local = TRUE)
source("./census_utils.R", local = TRUE)

# Load the saved voter-inference results and perform statistical analysis

load(file = "results.Rda")
load(file = "1855_inference.Rda")
load(file = "1855_covariate_inference.Rda")
load(file = "ct_demographics.Rda")
load(file = "yr_results.Rda")

# Calculate Meriden's rank among towns for Free Soil and Whig votes as a percentage
# of total votes cast; present in a table with vote tallies and percentages

table_1849_1851 <- data.frame(matrix(ncol = 9, nrow = 3))
years <- 1849:1851
stats <- c("Number of towns", "Whig share rank", "Free Soil share rank", "Whig votes", "Democratic votes", "Free soil votes", "Whig share", "Democratic share", "Free soil share")
rownames(table_1849_1851) <- years
colnames(table_1849_1851) <- stats
table_1849_1851[1, 1] <- nrow(vote_share.1849)
table_1849_1851[1, 2] <- which((vote_share.1849 %>% arrange(desc(Whig)))$town == "Meriden")
table_1849_1851[1, 3] <- which((vote_share.1849 %>% arrange(desc(Free_Soil)))$town == "Meriden")
table_1849_1851[2, 1] <- nrow(vote_share.1850)
table_1849_1851[2, 2] <- which((vote_share.1850 %>% arrange(desc(Whig)))$town == "Meriden")
table_1849_1851[2, 3] <- which((vote_share.1850 %>% arrange(desc(Free_Soil)))$town == "Meriden")
table_1849_1851[3, 1] <- nrow(vote_share.1851)
table_1849_1851[3, 2] <- which((vote_share.1851 %>% arrange(desc(Whig)))$town == "Meriden")
table_1849_1851[3, 3] <- which((vote_share.1851 %>% arrange(desc(Free_Soil)))$town == "Meriden")
table_1849_1851[1, 4:9] <- vote_share.1849 %>%
  dplyr::filter(town == "Meriden") %>%
  select(Whig_votes, Democrat_votes, Free_Soil_votes, Whig, Democrat, Free_Soil)
table_1849_1851[2, 4:9] <- vote_share.1850 %>%
  dplyr::filter(town == "Meriden") %>%
  select(Whig_votes, Democrat_votes, Free_Soil_votes, Whig, Democrat, Free_Soil)
table_1849_1851[3, 4:9] <- vote_share.1851 %>%
  dplyr::filter(town == "Meriden") %>%
  select(Whig_votes, Democrat_votes, Free_Soil_votes, Whig, Democrat, Free_Soil)
print(knitr::kable(table_1849_1851, caption = "Meriden election results and statewide rank, 1849-1851"))

# Calculate weighted means and standard deviations for results

tag_yr <- function(results, yr, type) {
  results %>%
    dplyr::filter(town == type) %>%
    select(-weight) %>%
    rename(year = town) %>%
    mutate(year = yr)
}

result_mu <- results.1849 %>%
  tag_yr(1849, "Weighted mean") %>%
  bind_rows(results.1850 %>% tag_yr(1850, "Weighted mean")) %>%
  bind_rows(results.1851 %>% tag_yr(1851, "Weighted mean")) %>%
  bind_rows(results.1852 %>% tag_yr(1852, "Weighted mean")) %>%
  bind_rows(results.1853 %>% tag_yr(1853, "Weighted mean")) %>%
  bind_rows(results.1854 %>% tag_yr(1854, "Weighted mean")) %>%
  bind_rows(results.1855 %>% tag_yr(1855, "Weighted mean")) %>%
  bind_rows(results.1856 %>% tag_yr(1856, "Weighted mean")) %>%
  bind_rows(results.1857 %>% tag_yr(1857, "Weighted mean"))

result_sd <- results.1849 %>%
  tag_yr(1849, "Weighted SD") %>%
  bind_rows(results.1850 %>% tag_yr(1850, "Weighted SD")) %>%
  bind_rows(results.1851 %>% tag_yr(1851, "Weighted SD")) %>%
  bind_rows(results.1852 %>% tag_yr(1852, "Weighted SD")) %>%
  bind_rows(results.1853 %>% tag_yr(1853, "Weighted SD")) %>%
  bind_rows(results.1854 %>% tag_yr(1854, "Weighted SD")) %>%
  bind_rows(results.1855 %>% tag_yr(1855, "Weighted SD")) %>%
  bind_rows(results.1856 %>% tag_yr(1856, "Weighted SD")) %>%
  bind_rows(results.1857 %>% tag_yr(1857, "Weighted SD"))

single_stat <- function(stats, yr, party) {
  stats %>%
    dplyr::filter(year == yr) %>%
    select(party)
}

# Calculate the Z-score for a town's party share in a given year
z_score <- function(results, mus, sds, yr, party, town) {
  town_result <- results %>%
    dplyr::filter(town == "Meriden") %>%
    select(party)
  mu <- single_stat(mus, yr, party)
  sd <- single_stat(sds, yr, party)
  round((town_result - mu) / sd, 2)
}

# Calculate the Z-scores for Meriden from 1849 - 1857
z_scores <- data.frame(matrix(ncol = 7, nrow = 9))
years <- 1849:1857
parties <- c("Democrat", "Whig", "Free_Soil", "Temperence", "Know_Nothing", "Republican", "Abstaining")
rownames(z_scores) <- years
colnames(z_scores) <- parties
results_dfs <- list(results.1849, results.1850, results.1851, results.1852, results.1853, results.1854, results.1855, results.1856, results.1857)
for (party in parties) {
  for (yr in years) {
    i <- yr - 1848
    j <- match(party, parties)
    yr_results <- results_dfs[[i]]
    yr_parties <- colnames(yr_results)
    if (party %in% yr_parties) {
      z_scores[i, j] <- z_score(yr_results, result_mu, result_sd, yr, party, "Meriden")
    }
  }
}

cat("\n\n         Z-scores for Meriden results")
print(z_scores)

# Meriden is the town that best demonstrates a rough statistical proxy
# for Tyler Anbinder's argument that the Know Nothings arose because of
# anti-Nebraska sentiment.

# Look for all towns that voted within 0.35 standard deviations of the
# 1853 and 1854 mean Democratic and Whig shares but voted more than
# 0.6 standard deviations above the 1854 Free Soil share and more
# than 0.6 standard deviations above the 1855 mean Know Nothing share.

cat("\n\n1855 transitions for towns with representative 1853 and 1854 Whig and Democratic votes but high 1854 Free Soil and 1855 Know Nothing votes")
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
    print(knitr::kable(construct_contingency(results.55, betas, 1854, 1855, t, 1),
      caption = paste(t, results.55$town[t], sep = ": ")
    ))
  }
}

# The greatest uncertainty in the 1855 inference is in determining what became of
# 1854 Free Soil voters, followed by doubt about Temperance and Whig voters.
# Similar results are obtained when standard deviations are calculated for the
# models alphas.

combined_betas <- combine_betas(beta.sims.MD(cov.ei.55, p55))
z_scores <- data.frame(beta_sd = numeric(), transition = character())
for (from_party in pull(distinct(combined_betas %>% select(from)))) {
  for (to_party in pull(distinct(combined_betas %>% select(to)))) {
    transition <- combined_betas %>%
      dplyr::filter(
        from == from_party,
        to == to_party
      ) %>%
      select(value)
    beta_sd <- sd(pull(transition))
    if (beta_sd > 2.5) {
      transition <- paste(from_party, " to ", to_party, sep = "")
      z_scores <- rbind(z_scores, data.frame(beta_sd = beta_sd, transition = transition))
    }
  }
}
cat("\n\n1855 transitions with high standard deviations in covariate model\n")
print(z_scores %>%
        arrange(desc(beta_sd)))

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
  dplyr::filter(BIRTH == "native" & SEX == 1)
ct_native_male_1860 <- ct_1860 %>%
  dplyr::filter(BIRTH == "native" & SEX == 1)
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
    dplyr::filter(town == t) %>%
    select(AGE) %>%
    arrange(AGE) %>%
    group_by(AGE) %>%
    summarize(n = n())
  native_male_1850.ss <- smooth.spline(native_male_1850)
  native_male_1860 <- ct_native_male_1860 %>%
    dplyr::filter(town == t) %>%
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

cat("\n\nRegression analysis of GINI index against demographic factors\n")
summary(lm(gini ~ pct_farm_1860 + age_1860 + wealth, factors))

# GINI index for native white males is strongly correlated to age, and age
# accounts for about 25% to the variance. Wealth is also correlated to age,
# and to GINI index, but it only accounts for about 5% of the variance in GINI

white_male_1860_hh <- ct_1860 %>%
  # Exclude servants and institutional housing
  dplyr::filter((GQ %in% c(1, 2, 5) && FAMUNIT == 1) || GQ == 4) %>%
  group_by(SERIAL * 100 + FAMUNIT) %>%
  summarise(
    FAMILY_REALPROP = sum(REALPROP),
    FAMILY_WEALTH = sum(WEALTH),
    AGE = first(AGE),
    AGE_CAT = first(AGE_CAT),
    BIRTH = first(BIRTH),
    SEX = first(SEX),
    RACE = first(RACE),
    JOB = first(JOB),
    CLASS = first(CLASS),
    town = first(town),
    combined = first(combined)
  ) %>%
  dplyr::filter(SEX == 1 && RACE == 1)

# Use the average age for each of the age categories when doing regression
avg_ages <- c()
age_cats <- c("20 - 29", "30 - 39", "40 - 49", "50 - 59", "60 and over")
for (age in age_cats) {
  avg_ages <- c(avg_ages, as.numeric(white_male_1860_hh %>%
    dplyr::filter(AGE_CAT == age) %>%
    dplyr::filter(BIRTH == "native") %>%
    summarise(avg_age = mean(AGE)) %>%
    select(avg_age)))
}

gini_by_town <- distinct(ct_1860_hh %>% select(town))
for (age in age_cats) {
  avg_age <- avg_ages[1]
  avg_ages <- avg_ages[2:length(avg_ages)]
  g <- white_male_1860_hh %>%
    dplyr::filter(AGE_CAT == age) %>%
    dplyr::filter(BIRTH == "native") %>%
    group_by(town) %>%
    summarise(
      gini = ineq::ineq(FAMILY_WEALTH, type = "Gini"),
      avg_wealth = mean(FAMILY_WEALTH)
    ) %>%
    select(town, gini, avg_wealth)
  colnames(g) <- c("town", avg_age, paste(avg_age, "wealth", sep = "_"))
  gini_by_town <- gini_by_town %>% left_join(g, by = join_by(town))
}

gini_age <- gini_by_town %>%
  select(!ends_with("wealth")) %>%
  pivot_longer(cols = !town, names_to = "age", values_to = "gini") %>%
  mutate(age = as.double(age))

wealth_age <- gini_by_town %>%
  select(town, ends_with("wealth")) %>%
  pivot_longer(cols = !town, names_to = "age", values_to = "wealth") %>%
  mutate(age = as.double(substring(age, 1, nchar(age) - 7)))

# Print out the proportion of variance in GINI index accounted for by age and wealth
cat("\n\nProportion of variance in GINI index accounted for by age and wealth\n")
lsr::etaSquared(lm(gini ~ wealth + age, gini_age %>% left_join(wealth_age, by = join_by(town, age))))

# Show GINI indices for Meriden subgroups

age_cats <- c("0 - 19", age_cats)
gini_all <- data.frame(row.names = c("all"))
gini_by_birth <- data.frame(row.names = c("native", "immigrant"))
gini_by_job <- data.frame(row.names = c("farm", "nonfarm"))
for (age in age_cats) {
  g <- ct_1860_hh %>%
    dplyr::filter(
      town == "Meriden",
      AGE_CAT == age
    ) %>%
    summarise(gini = ineq::ineq(FAMILY_WEALTH, type = "Gini")) %>%
    select(gini)
  gini_all <- bind_cols(gini_all, g, .name_repair = "unique_quiet")
}
for (age in age_cats) {
  g <- ct_1860_hh %>%
    dplyr::filter(
      town == "Meriden",
      AGE_CAT == age
    ) %>%
    group_by(BIRTH) %>%
    summarise(gini = ineq::ineq(FAMILY_WEALTH, type = "Gini")) %>%
    select(gini)
  gini_by_birth <- bind_cols(gini_by_birth, g, .name_repair = "unique_quiet")
}
for (age in age_cats) {
  g <- ct_1860_hh %>%
    dplyr::filter(
      town == "Meriden",
      AGE_CAT == age
    ) %>%
    group_by(JOB) %>%
    summarise(gini = ineq::ineq(FAMILY_WEALTH, type = "Gini")) %>%
    select(gini)
  gini_by_job <- bind_cols(gini_by_job, g, .name_repair = "unique_quiet")
}
gini_grid <- bind_rows(gini_all, gini_by_birth, gini_by_job)
colnames(gini_grid) <- age_cats
cat("\n\n1860 Meriden Wealth GINI Index by Age Range")
knitr::kable(gini_grid, label = "1860 Meriden Wealth GINI Index by Age Range")

# Calculate real-property-only GINI indices for 1850 and 1860
gini_all <- data.frame(row.names = c("all"))
gini_by_birth <- data.frame(row.names = c("native", "immigrant"))
gini_by_job <- data.frame(row.names = c("farm", "nonfarm"))
for (age in age_cats) {
  g <- ct_1850_hh %>%
    dplyr::filter(
      town == "Meriden",
      AGE_CAT == age
    ) %>%
    summarise(gini = ineq::ineq(FAMILY_REALPROP, type = "Gini")) %>%
    select(gini)
  gini_all <- bind_cols(gini_all, g, .name_repair = "unique_quiet")
}
for (age in age_cats) {
  g <- ct_1850_hh %>%
    dplyr::filter(
      town == "Meriden",
      AGE_CAT == age
    ) %>%
    group_by(BIRTH) %>%
    summarise(gini = ineq::ineq(FAMILY_REALPROP, type = "Gini")) %>%
    select(gini)
  gini_by_birth <- bind_cols(gini_by_birth, g, .name_repair = "unique_quiet")
}
for (age in age_cats) {
  g <- ct_1850_hh %>%
    dplyr::filter(
      town == "Meriden",
      AGE_CAT == age
    ) %>%
    group_by(JOB) %>%
    summarise(gini = ineq::ineq(FAMILY_REALPROP, type = "Gini")) %>%
    select(gini)
  gini_by_job <- bind_cols(gini_by_job, g, .name_repair = "unique_quiet")
}
gini_grid <- bind_rows(gini_all, gini_by_birth, gini_by_job)
colnames(gini_grid) <- age_cats
cat("\n\n1850 Meriden Real Property GINI Index by Age Range")
knitr::kable(gini_grid, label = "1850 Meriden Real Property GINI Index by Age Range")

gini_all <- data.frame(row.names = c("all"))
gini_by_birth <- data.frame(row.names = c("native", "immigrant"))
gini_by_job <- data.frame(row.names = c("farm", "nonfarm"))
for (age in age_cats) {
  g <- ct_1860_hh %>%
    dplyr::filter(
      town == "Meriden",
      AGE_CAT == age
    ) %>%
    summarise(gini = ineq::ineq(FAMILY_REALPROP, type = "Gini")) %>%
    select(gini)
  gini_all <- bind_cols(gini_all, g, .name_repair = "unique_quiet")
}
for (age in age_cats) {
  g <- ct_1860_hh %>%
    dplyr::filter(
      town == "Meriden",
      AGE_CAT == age
    ) %>%
    group_by(BIRTH) %>%
    summarise(gini = ineq::ineq(FAMILY_REALPROP, type = "Gini")) %>%
    select(gini)
  gini_by_birth <- bind_cols(gini_by_birth, g, .name_repair = "unique_quiet")
}
for (age in age_cats) {
  g <- ct_1860_hh %>%
    dplyr::filter(
      town == "Meriden",
      AGE_CAT == age
    ) %>%
    group_by(JOB) %>%
    summarise(gini = ineq::ineq(FAMILY_REALPROP, type = "Gini")) %>%
    select(gini)
  gini_by_job <- bind_cols(gini_by_job, g, .name_repair = "unique_quiet")
}
gini_grid <- bind_rows(gini_all, gini_by_birth, gini_by_job)
colnames(gini_grid) <- age_cats
cat("\n\n1860 Meriden Real Property GINI Index by Age Range")
knitr::kable(gini_grid, label = "1860 Meriden Real Property GINI Index by Age Range")

# Demonstrate the extent to which results for statewide offices are correlated
# This is to be expected because voters are given a single ballot with all offices
# listed on it, so their choices for different offices are not independent

governor_results <- read_results(results_file, alias_map, party_assignments, "Governor", 1849, 1857)
lt_governor_results <- read_results(results_file, alias_map, party_assignments, "Lieutenant Governor", 1849, 1857)
secretary_results <- read_results(results_file, alias_map, party_assignments, "Secretary of the State", 1849, 1857)
treasurer_results <- read_results(results_file, alias_map, party_assignments, "Treasurer", 1849, 1857)
probate_results <- read_probate_results(results_file, alias_map, party_assignments, 1849, 1857)

# Keep only those rows that have results for all five offices
governor_results <- semi_join(governor_results, probate_results, by = c("town", "yr")) %>%
  semi_join(lt_governor_results, by = c("town", "yr")) %>%
  semi_join(secretary_results, by = c("town", "yr")) %>%
  semi_join(treasurer_results, by = c("town", "yr"))
lt_governor_results <- semi_join(lt_governor_results, probate_results, by = c("town", "yr")) %>%
  semi_join(governor_results, by = c("town", "yr")) %>%
  semi_join(secretary_results, by = c("town", "yr")) %>%
  semi_join(treasurer_results, by = c("town", "yr"))
secretary_results <- semi_join(secretary_results, probate_results, by = c("town", "yr")) %>%
  semi_join(governor_results, by = c("town", "yr")) %>%
  semi_join(lt_governor_results, by = c("town", "yr")) %>%
  semi_join(treasurer_results, by = c("town", "yr"))
treasurer_results <- semi_join(treasurer_results, probate_results, by = c("town", "yr")) %>%
  semi_join(governor_results, by = c("town", "yr")) %>%
  semi_join(lt_governor_results, by = c("town", "yr")) %>%
  semi_join(secretary_results, by = c("town", "yr"))
probate_results <- semi_join(probate_results, governor_results, by = c("town", "yr")) %>%
  semi_join(lt_governor_results, by = c("town", "yr")) %>%
  semi_join(secretary_results, by = c("town", "yr")) %>%
  semi_join(treasurer_results, by = c("town", "yr"))

parties <- c("Democrat", "Whig", "Free_Soil", "Temperance", "Know_Nothing", "Republican")
for (party in parties) {
  corr_df <- corr_matrix_by_party(governor_results, lt_governor_results, secretary_results, treasurer_results, probate_results, party)
  print(knitr::kable(corr_df, caption = paste("Correlation of", party, "vote tallies for statewide offices, 1849-1857")))
}
