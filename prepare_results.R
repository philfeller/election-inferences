# Create and save election results for Connecticut towns from 1849 to 1857

source("./global.R")
source("./results_utils.R", local = TRUE)

raw_results <- read_results(results_file, alias_map, party_assignments, "Governor", 1849, 1857)

# Because estimates of eligible voters depend on both 1850 and 1860 census data,
# they must be calculated for town combinations that will be less granular than
# those used for ecological inference. Calculate for each of these combinations
# what percentage of eligible voters cast ballots in a given year. These can be
# used to estimate the number of eligible voters for each of the towns.
# In a few cases, estimates of eligible voters fall short of the total votes
# cast; set the percentage to 100% in these cases.
eligible_pct <- raw_results %>%
  select(yr, combined, total) %>%
  pivot_wider(names_prefix = "TOTAL_", names_from = yr, values_from = total, values_fn = sum, values_fill = 0) %>%
  left_join(ct_eligible, by = "combined") %>%
  mutate(
    ELIG_1849_PCT = sapply(TOTAL_1849 / ELIG_1850, cap),
    ELIG_1850_PCT = sapply(TOTAL_1850 / ELIG_1850, cap),
    ELIG_1851_PCT = sapply(TOTAL_1851 / ELIG_1851, cap),
    ELIG_1852_PCT = sapply(TOTAL_1852 / ELIG_1852, cap),
    ELIG_1853_PCT = sapply(TOTAL_1853 / ELIG_1853, cap),
    ELIG_1854_PCT = sapply(TOTAL_1854 / ELIG_1854, cap),
    ELIG_1855_PCT = sapply(TOTAL_1855 / ELIG_1855, cap),
    ELIG_1856_PCT = sapply(TOTAL_1856 / ELIG_1856, cap),
    ELIG_1857_PCT = sapply(TOTAL_1857 / ELIG_1857, cap)
  ) %>%
  select(combined, ends_with("_PCT"))

# Calculate vote shares for 1849 to 1851

vote_share.1849 <- yr_results(raw_results, 1849) %>%
  mutate(
    Free_Soil = Free_Soil_votes / total,
    Whig = Whig_votes / total,
    Democrat = Democrat_votes / total
  )

vote_share.1850 <- yr_results(raw_results, 1850) %>%
  mutate(
    Free_Soil = Free_Soil_votes / total,
    Whig = Whig_votes / total,
    Democrat = Democrat_votes / total
  )

vote_share.1851 <- yr_results(raw_results, 1851) %>%
  mutate(
    Free_Soil = Free_Soil_votes / total,
    Whig = Whig_votes / total,
    Democrat = Democrat_votes / total
  )

# Calculate results for individual years, in order to calculate a town's
# z-score for particular results
results.1849 <- yr_results(raw_results, 1849) %>%
  mutate(
    ELIG = round(total / ELIG_1850_PCT),
    weight = ELIG / sum(ELIG),
    Democrat = Democrat_votes / ELIG,
    Whig = Whig_votes / ELIG,
    Free_Soil = Free_Soil_votes / ELIG,
    Abstaining = remainder(Democrat + Whig + Free_Soil)
  ) %>%
  result_summary()

results.1850 <- yr_results(raw_results, 1850) %>%
  mutate(
    ELIG = round(total / ELIG_1850_PCT),
    weight = ELIG / sum(ELIG),
    Democrat = Democrat_votes / ELIG,
    Whig = Whig_votes / ELIG,
    Free_Soil = Free_Soil_votes / ELIG,
    Abstaining = remainder(Democrat + Whig + Free_Soil)
  ) %>%
  result_summary()

results.1851 <- yr_results(raw_results, 1851) %>%
  mutate(
    ELIG = round(total / ELIG_1851_PCT),
    weight = ELIG / sum(ELIG),
    Democrat = Democrat_votes / ELIG,
    Whig = Whig_votes / ELIG,
    Free_Soil = Free_Soil_votes / ELIG,
    Abstaining = remainder(Democrat + Whig + Free_Soil)
  ) %>%
  result_summary()

results.1852 <- yr_results(raw_results, 1852) %>%
  mutate(
    ELIG = round(total / ELIG_1852_PCT),
    weight = ELIG / sum(ELIG),
    Democrat = Democrat_votes / ELIG,
    Whig = Whig_votes / ELIG,
    Free_Soil = Free_Soil_votes / ELIG,
    Abstaining = remainder(Democrat + Whig + Free_Soil)
  ) %>%
  result_summary()

results.1853 <- yr_results(raw_results, 1853) %>%
  mutate(
    ELIG = round(total / ELIG_1853_PCT),
    weight = ELIG / sum(ELIG),
    Democrat = Democrat_votes / ELIG,
    Whig = Whig_votes / ELIG,
    Free_Soil = Free_Soil_votes / ELIG,
    Abstaining = remainder(Democrat + Whig + Free_Soil)
  ) %>%
  result_summary()

results.1854 <- yr_results(raw_results, 1854) %>%
  mutate(
    ELIG = round(total / ELIG_1854_PCT),
    weight = ELIG / sum(ELIG),
    Democrat = Democrat_votes / ELIG,
    Whig = Whig_votes / ELIG,
    Free_Soil = Free_Soil_votes / ELIG,
    Temperence = Temperance_votes / ELIG,
    Abstaining = remainder(Democrat + Whig + Free_Soil + Temperence)
  ) %>%
  result_summary()

results.1855 <- yr_results(raw_results, 1855) %>%
  mutate(
    ELIG = round(total / ELIG_1855_PCT),
    weight = ELIG / sum(ELIG),
    Democrat = Democrat_votes / ELIG,
    Whig = Whig_votes / ELIG,
    Know_Nothing = Know_Nothing_votes / ELIG,
    Abstaining = remainder(Democrat + Whig + Know_Nothing)
  ) %>%
  result_summary()

results.1856 <- yr_results(raw_results, 1856) %>%
  mutate(
    ELIG = round(total / ELIG_1856_PCT),
    weight = ELIG / sum(ELIG),
    Democrat = Democrat_votes / ELIG,
    Whig = Whig_votes / ELIG,
    Know_Nothing = Know_Nothing_votes / ELIG,
    Republican = Republican_votes / ELIG,
    Abstaining = remainder(Democrat + Whig + Know_Nothing + Republican)
  ) %>%
  result_summary()

results.1857 <- yr_results(raw_results, 1857) %>%
  mutate(
    ELIG = round(total / ELIG_1857_PCT),
    weight = ELIG / sum(ELIG),
    Democrat = Democrat_votes / ELIG,
    Republican = Republican_votes / ELIG,
    Abstaining = remainder(Democrat + Republican)
  ) %>%
  result_summary()

save(
  vote_share.1849, vote_share.1850, vote_share.1851,
  results.1849, results.1850, results.1851, results.1852, results.1853, results.1854, results.1855, results.1856, results.1857,
  file = "yr_results.Rda"
)
