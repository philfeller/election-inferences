# Define a function that will be used to create a tibble with gubernatorial
# results for a pair of years:
#   create_results - prepares a tibble with percentage election results

library(sf)

# Define constants
source("./variables.R")

# Estimate nonvoters and calculate potential covariates
source("./ct_demographics.R")

# Define functions to use when combining election results from multiple towns
combine_towns <- function(input_tibble, towns, combined_name) {
  return(input_tibble %>%
    filter(town %in% towns) %>%
    group_by(yr) %>%
    summarize(across(where(is.double), sum)) %>%
    add_column(town = combined_name, .after = "yr"))
}

combine_results <- function(input_tibble, beg_yr, end_yr) {
  # Combine results to account for new towns that were incorporated during time period
  # https://portal.ct.gov/SOTS/Register-Manual/Section-VII/Connecticut-Towns-in-the-Order-of-their-Establishment

  # When comparing pre-1856 results to 1856 or 1857, combine towns that split in 1856
  if (beg_yr < 1856 && end_yr >= 1856) {
    # Putnam, taken from Killingly, Thompson, and Pomfret, in is 1856 and 1857 results
    kp_combined <- combine_towns(input_tibble, killingly_et_al, "Killingly, Thompson, and Pomfret")
    # Bethel, taken from Danbury, is in 1856 and 1857 results
    d_combined <- combine_towns(input_tibble, danbury, "Danbury")
    # South Lyme, taken from Lyme, is in 1856 and 1857 results
    l_combined <- combine_towns(input_tibble, lyme, "Lyme")
  } else {
    kp_combined <- input_tibble %>% filter(town %in% killingly_et_al)
    d_combined <- input_tibble %>% filter(town %in% danbury)
    l_combined <- input_tibble %>% filter(town %in% lyme)
  }

  # When comparing pre-1855 results to 1855 or later, combine towns that split in 1855
  if (beg_yr < 1855 && end_yr >= 1855) {
    # West Hartford, taken from Hartford, is in 1855, 1856, and 1857 results
    h_combined <- combine_towns(input_tibble, hartford, "Hartford")
    # Windsor Locks taken from Windsor in 1855, 1856, and 1857 results
    w_combined <- combine_towns(input_tibble, windsor, "Windsor")
  } else {
    h_combined <- input_tibble %>% filter(town %in% hartford)
    w_combined <- input_tibble %>% filter(town %in% windsor)
  }

  # When comparing pre-1852 results to 1852 or later, combine town that split in 1851
  if (beg_yr < 1852 && end_yr >= 1852) {
    # Cromwell, taken from Middletown, is in results after 1851
    m_combined <- combine_towns(input_tibble, middletown, "Middletown")
  } else {
    m_combined <- input_tibble %>% filter(town %in% middletown)
  }

  # Division and renaming of parts of Saybrook affect results beginning in 1853
  os_combined <- combine_towns(input_tibble, saybrook, "Saybrook")
  # Bridgewater , taken from New Milford, is in 1857 results
  nm_combined <- combine_towns(input_tibble, new_milford, "New Milford")

  return(input_tibble %>%
    filter(!town %in% middletown) %>%
    add_row(m_combined) %>%
    filter(!town %in% saybrook) %>%
    add_row(os_combined) %>%
    filter(!town %in% hartford) %>%
    add_row(h_combined) %>%
    filter(!town %in% windsor) %>%
    add_row(w_combined) %>%
    filter(!town %in% danbury) %>%
    add_row(d_combined) %>%
    filter(!town %in% lyme) %>%
    add_row(l_combined) %>%
    filter(!town %in% killingly_et_al) %>%
    add_row(kp_combined) %>%
    filter(!town %in% new_milford) %>%
    add_row(nm_combined))
}

# Define function to assign party designation to candidates
assign_party <- function(input_tibble, other = FALSE) {
  # Nominees of major parties between 1852 and 1857
  major_candidates <- c(
    "Alexander H. Holley",
    "Charles Chapman",
    "Gideon Welles",
    "Henry Dutton",
    "Lafayette S. Foster",
    "Green Kendrick",
    "John A. Rockwell",
    "John Hooker",
    "Francis Gillette",
    "John Boyd",
    "Samuel Ingham",
    "Thomas H. Seymour",
    "William T. Minor"
  )

  # Alternate names used to vote for Seymour, Ingham, and Foster
  seymour <- c("T. H. Seymour", "T. Seymour", "Thomas Seymour", "Thomas S. Seymour")
  ingham <- c("Samuel D. Ingham", "Samuel Engraham")
  foster <- c("L. S. Foster", "Lafayette Foster")

  # Candidates whose names appeared on stray ballots when not their party's nominee
  stray_votes <- input_tibble %>%
    filter(yr == 1856 && candidate_name == "Henry Dutton") %>%
    filter(yr > 1853 && candidate_name == "Francis Gillette") %>%
    filter(yr == 1853 && candidate_name == "Samuel Ingham") %>%
    filter(yr == 1857 && candidate_name == "William T. Minor")

  minor_votes <- input_tibble %>%
    filter(!candidate_name %in% major_candidates) %>%
    add_row(stray_votes) %>%
    group_by(yr, town) %>%
    summarize(across(where(is.double), sum)) %>%
    add_column(candidate_party = "Other_votes", .after = "yr")

  major_votes <- input_tibble %>%
    mutate(candidate_name = ifelse(candidate_name %in% seymour, "Thomas H. Seymour", candidate_name)) %>%
    mutate(candidate_name = ifelse(candidate_name %in% foster, "Lafayette S. Foster", candidate_name)) %>%
    mutate(candidate_name = ifelse(candidate_name %in% ingham, "Samuel Ingham", candidate_name)) %>%
    filter(candidate_name %in% major_candidates) %>%
    # Remove candidates whose names appeared on stray ballots when their not party's nominee
    filter(!(yr == 1856 && candidate_name == "Henry Dutton")) %>%
    filter(!(yr > 1853 && candidate_name == "Francis Gillette")) %>%
    filter(!(yr == 1853 && candidate_name == "Samuel Ingham")) %>%
    filter(!(yr == 1857 && candidate_name == "William T. Minor")) %>%
    mutate(candidate_party = case_when(
      candidate_name == "Alexander H. Holley" ~ "Republican_votes",
      candidate_name == "Charles Chapman" ~ "Temperance_votes",
      candidate_name == "Gideon Welles" ~ "Republican_votes",
      candidate_name == "Lafayette S. Foster" ~ "Whig_votes",
      candidate_name == "Green Kendrick" ~ "Whig_votes",
      candidate_name == "Henry Dutton" ~ "Whig_votes",
      candidate_name == "John A. Rockwell" ~ "Whig_votes",
      candidate_name == "John Hooker" ~ "Free_Soil_votes",
      candidate_name == "Francis Gillette" ~ "Free_Soil_votes",
      candidate_name == "John Boyd" ~ "Free_Soil_votes",
      candidate_name == "Samuel Ingham" ~ "Democrat_votes",
      candidate_name == "Thomas H. Seymour" ~ "Democrat_votes",
      candidate_name == "William T. Minor" ~ "Know_Nothing_votes"
    )) %>%
    mutate(candidate_name = NULL)

  if (other) {
    return(major_votes %>%
      add_row(minor_votes))
  } else {
    return(major_votes)
  }
}

# Create a tibble with demographic factors for a particular year range,
# appropriately combining towns that separated during the period
create_factors <- function(tibble, beg_yr, end_yr) {
  towns <- c()
  if (beg_yr >= 1852) {
    towns <- c(towns, middletown)
  }
  if (beg_yr >= 1855) {
    towns <- c(towns, hartford, windsor)
  }
  if (beg_yr >= 1856) {
    towns <- c(towns, danbury, lyme)
  }
  if (end_yr <= 1855 || beg_yr >= 1856) {
    towns <- c(towns, killingly_et_al)
  }
  tibble %>%
    mutate(
      combined = ifelse(town %in% towns, town, combined),
      combined = ifelse(combined == "Old Lyme", "South Lyme", combined),
      gini = ifelse(town %in% towns, gini, comb_gini),
      wealth = ifelse(town %in% towns, wealth, comb_wealth),
      age_1860 = ifelse(town %in% towns, age_1860, comb_age_1860),
      town = NULL,
      comb_gini = NULL,
      comb_wealth = NULL,
      comb_age_1850 = NULL,
      comb_age_1860 = NULL
    ) %>%
    distinct(.keep_all = TRUE)
}

# Read and process the downloaded election results file
read_results <- function(file) {
  return(
    read_csv(file, show_col_types = FALSE) %>%
      select(
        election_date,
        candidate_name,
        candidate_party,
        granular_division_name,
        votes
      ) %>%
      rename(town = granular_division_name) %>%
      mutate(
        yr = as.integer(year(election_date)),
        election_date = NULL, .before = candidate_name
      ) %>%
      filter(!candidate_name == "Total Ballots Cast") %>%
      filter(!candidate_name == "Total Votes Cast") %>%
      filter(yr >= 1851) %>%
      assign_party() %>%
      # Change blank votes to zero
      replace(is.na(.), 0) %>%
      # Create tibble that combines all town results into one row
      pivot_wider(
        names_from = candidate_party,
        values_from = votes,
        values_fn = sum,
        values_fill = 0
      ) %>%
      # Create a column with total votes for town election
      mutate(
        total = rowSums(across(where(is.double))),
        combined = get_combined(town)
      )
  )
}

raw_results <- read_results(results_file)

# Helper function to cap results at 1
cap <- function(val) {
  if (val > 1) {
    return(1)
  } else {
    return(val)
  }
}

# Helper function to obtain the capped remainder
remainder <- function(arry) {
  if (sum(is.na(arry)) == 0) {
    return(1 - sapply(arry, cap))
  } else {
    return(arry)
  }
}

# Because estimates of eligible voters depend on both 1850 and 1860 census data,
# they must be calculated for town combinations that will be less granular than
# those used for ecological inference. Calculate for each of these combinations
# what percentage of eligible voters cast ballots in a given year. These can be
# used to estimate the number of eligible voters for each of the towns.
# In a few cases, estimates of eligible voters fall short of the total votes
# cast; set the percentage to 100% in these cases.
eligible_pct <- raw_results %>%
  select(yr, combined, total) %>%
  pivot_wider(names_prefix = "TOTAL_", names_from = yr, values_from = total, values_fn = sum) %>%
  left_join(ct_eligible, by = "combined") %>%
  mutate(
    ELIG_1851_PCT = sapply(TOTAL_1851 / ELIG_1851, cap),
    ELIG_1852_PCT = sapply(TOTAL_1852 / ELIG_1852, cap),
    ELIG_1853_PCT = sapply(TOTAL_1853 / ELIG_1853, cap),
    ELIG_1854_PCT = sapply(TOTAL_1854 / ELIG_1854, cap),
    ELIG_1855_PCT = sapply(TOTAL_1855 / ELIG_1855, cap),
    ELIG_1856_PCT = sapply(TOTAL_1856 / ELIG_1856, cap),
    ELIG_1857_PCT = sapply(TOTAL_1857 / ELIG_1857, cap)
  ) %>%
  select(combined, ends_with("_PCT"))

# Create a results tibble that is appropriate for the years to be compared;
# don't combine election results for towns unless necessary.
# This is particularly important in Windham County, where Killingly, Pomfret,
# and Thompson were split to form Putnam. Unless comparing to 1856 or 1857
# results, combining these towns results in a significant loss of detail.
create_results <- function(beg_yr, end_yr) {
  # Create separate tibbles for each year
  e51 <- raw_results %>%
    filter(yr == 1851) %>%
    combine_results(beg_yr, end_yr) %>%
    filter(yr >= beg_yr && yr <= end_yr) %>%
    mutate(
      yr = NULL,
      Temperance_votes = NULL,
      Know_Nothing_votes = NULL,
      Republican_votes = NULL
    ) %>%
    rename(
      Democrat_vote_in_1851 = Democrat_votes,
      Whig_vote_in_1851 = Whig_votes,
      Free_Soil_vote_in_1851 = Free_Soil_votes,
      total_1851 = total
    )
  e52 <- raw_results %>%
    filter(yr == 1852) %>%
    combine_results(beg_yr, end_yr) %>%
    filter(yr >= beg_yr && yr <= end_yr) %>%
    mutate(
      yr = NULL,
      Temperance_votes = NULL,
      Know_Nothing_votes = NULL,
      Republican_votes = NULL
    ) %>%
    rename(
      Democrat_vote_in_1852 = Democrat_votes,
      Whig_vote_in_1852 = Whig_votes,
      Free_Soil_vote_in_1852 = Free_Soil_votes,
      total_1852 = total
    )

  e53 <- raw_results %>%
    filter(yr == 1853) %>%
    combine_results(beg_yr, end_yr) %>%
    filter(yr >= beg_yr && yr <= end_yr) %>%
    mutate(
      yr = NULL,
      Temperance_votes = NULL,
      Know_Nothing_votes = NULL,
      Republican_votes = NULL
    ) %>%
    rename(
      Democrat_vote_in_1853 = Democrat_votes,
      Whig_vote_in_1853 = Whig_votes,
      Free_Soil_vote_in_1853 = Free_Soil_votes,
      total_1853 = total
    )
  e54 <- raw_results %>%
    filter(yr == 1854) %>%
    combine_results(beg_yr, end_yr) %>%
    filter(yr >= beg_yr && yr <= end_yr) %>%
    mutate(
      yr = NULL,
      Know_Nothing_votes = NULL,
      Republican_votes = NULL
    ) %>%
    rename(
      Democrat_vote_in_1854 = Democrat_votes,
      Whig_vote_in_1854 = Whig_votes,
      Free_Soil_vote_in_1854 = Free_Soil_votes,
      Temperance_vote_in_1854 = Temperance_votes,
      total_1854 = total
    )
  e55 <- raw_results %>%
    filter(yr == 1855) %>%
    combine_results(beg_yr, end_yr) %>%
    filter(yr >= beg_yr && yr <= end_yr) %>%
    mutate(
      yr = NULL,
      Temperance_votes = NULL,
      Free_Soil_votes = NULL,
      Republican_votes = NULL
    ) %>%
    rename(
      Democrat_vote_in_1855 = Democrat_votes,
      Whig_vote_in_1855 = Whig_votes,
      Know_Nothing_vote_in_1855 = Know_Nothing_votes,
      total_1855 = total
    )
  e56 <- raw_results %>%
    filter(yr == 1856) %>%
    combine_results(beg_yr, end_yr) %>%
    filter(yr >= beg_yr && yr <= end_yr) %>%
    mutate(
      yr = NULL,
      Temperance_votes = NULL,
      Free_Soil_votes = NULL
    ) %>%
    rename(
      Democrat_vote_in_1856 = Democrat_votes,
      Whig_vote_in_1856 = Whig_votes,
      Know_Nothing_vote_in_1856 = Know_Nothing_votes,
      Republican_vote_in_1856 = Republican_votes,
      total_1856 = total
    )
  e57 <- raw_results %>%
    filter(yr == 1857) %>%
    combine_results(beg_yr, end_yr) %>%
    filter(yr >= beg_yr && yr <= end_yr) %>%
    mutate(
      yr = NULL,
      Whig_votes = NULL,
      Temperance_votes = NULL,
      Free_Soil_votes = NULL,
      Know_Nothing_votes = NULL
    ) %>%
    rename(
      Democrat_vote_in_1857 = Democrat_votes,
      Republican_vote_in_1857 = Republican_votes,
      total_1857 = total
    )

  # Choose the appropriate shapefile from which to get longitude
  shapefile <- "./maps/1851_CT_towns.shp"
  if (beg_yr >= 1852) {
    shapefile <- "./maps/1852_CT_towns.shp"
  }
  if (beg_yr >= 1854) {
    shapefile <- "./maps/1854_CT_towns.shp"
  }
  if (beg_yr >= 1855) {
    shapefile <- "./maps/1855_CT_towns.shp"
  }
  if (end_yr >= 1856) {
    shapefile <- "./maps/1856_CT_towns.shp"
  }
  longitude <- as.data.frame(read_sf(shapefile)) %>%
    select(TOWN_NAME, LON) %>%
    rename(town = TOWN_NAME) %>%
    mutate(town = ifelse(town == "Putnam" & beg_yr < 1856, "Killingly, Thompson, and Pomfret", town))

  # Generate the demographic factors appropriate for the range of years
  demo_factors <- factors %>%
    create_factors(beg_yr, end_yr) %>%
    rename(town = combined)

  # Combine the year-by-year tibbles to have one record for each town
  full_results <- e51 %>%
    full_join(e52, by = "town") %>%
    full_join(e53, by = "town") %>%
    full_join(e54, by = "town") %>%
    full_join(e55, by = "town") %>%
    full_join(e56, by = "town") %>%
    full_join(e57, by = "town") %>%
    filter(!if_all(everything(), is.na)) %>%
    left_join(longitude, by = "town") %>%
    mutate(combined = get_combined(town)) %>%
    left_join(eligible_pct, by = "combined") %>%
    left_join(demo_factors, by = "town") %>%
    arrange(town) %>%
    mutate(
      ELIG_1851 = round(total_1851 / ELIG_1851_PCT),
      nonvote_in_1851 = ELIG_1851 - total_1851,
      Democrat_in_1851 = Democrat_vote_in_1851 / ELIG_1851,
      Whig_in_1851 = Whig_vote_in_1851 / ELIG_1851,
      Free_Soil_in_1851 = Free_Soil_vote_in_1851 / ELIG_1851,
      Abstaining_in_1851 = remainder(Democrat_in_1851 + Whig_in_1851 + Free_Soil_in_1851),
      ELIG_1852 = round(total_1852 / ELIG_1852_PCT),
      nonvote_in_1852 = ELIG_1852 - total_1852,
      Democrat_in_1852 = Democrat_vote_in_1852 / ELIG_1852,
      Whig_in_1852 = Whig_vote_in_1852 / ELIG_1852,
      Free_Soil_in_1852 = Free_Soil_vote_in_1852 / ELIG_1852,
      Abstaining_in_1852 = remainder(Democrat_in_1852 + Whig_in_1852 + Free_Soil_in_1852),
      ELIG_1853 = round(total_1853 / ELIG_1853_PCT),
      nonvote_in_1853 = ELIG_1853 - total_1853,
      Democrat_in_1853 = Democrat_vote_in_1853 / ELIG_1853,
      Whig_in_1853 = Whig_vote_in_1853 / ELIG_1853,
      Free_Soil_in_1853 = Free_Soil_vote_in_1853 / ELIG_1853,
      Abstaining_in_1853 = remainder(Democrat_in_1853 + Whig_in_1853 + Free_Soil_in_1853),
      ELIG_1854 = round(total_1854 / ELIG_1854_PCT),
      nonvote_in_1854 = ELIG_1854 - total_1854,
      Democrat_in_1854 = Democrat_vote_in_1854 / ELIG_1854,
      Whig_in_1854 = Whig_vote_in_1854 / ELIG_1854,
      Free_Soil_in_1854 = Free_Soil_vote_in_1854 / ELIG_1854,
      Temperance_in_1854 = Temperance_vote_in_1854 / ELIG_1854,
      Abstaining_in_1854 = remainder(Democrat_in_1854 + Whig_in_1854 + Free_Soil_in_1854 + Temperance_in_1854),
      ELIG_1855 = round(total_1855 / ELIG_1855_PCT),
      nonvote_in_1855 = ELIG_1855 - total_1855,
      Democrat_in_1855 = Democrat_vote_in_1855 / ELIG_1855,
      Whig_in_1855 = Whig_vote_in_1855 / ELIG_1855,
      Know_Nothing_in_1855 = Know_Nothing_vote_in_1855 / ELIG_1855,
      Abstaining_in_1855 = remainder(Democrat_in_1855 + Whig_in_1855 + Know_Nothing_in_1855),
      ELIG_1856 = round(total_1856 / ELIG_1856_PCT),
      nonvote_in_1856 = ELIG_1856 - total_1856,
      Democrat_in_1856 = Democrat_vote_in_1856 / ELIG_1856,
      Republican_in_1856 = Republican_vote_in_1856 / ELIG_1856,
      Know_Nothing_in_1856 = Know_Nothing_vote_in_1856 / ELIG_1856,
      Whig_in_1856 = Whig_vote_in_1856 / ELIG_1856,
      Abstaining_in_1856 = remainder(Democrat_in_1856 + Republican_in_1856 + Know_Nothing_in_1856 + Whig_in_1856),
      ELIG_1857 = round(total_1857 / ELIG_1857_PCT),
      nonvote_in_1857 = ELIG_1857 - total_1857,
      Democrat_in_1857 = Democrat_vote_in_1857 / ELIG_1857,
      Republican_in_1857 = Republican_vote_in_1857 / ELIG_1857,
      Abstaining_in_1857 = remainder(Democrat_in_1857 + Republican_in_1857)
    ) %>%
    select_if(function(x) !any(is.na(x)))

  return(full_results)
}

# Calculate results for individual years, in order to calculate a town's
# z-score for particular results
results.1851 <- raw_results %>%
  filter(yr == 1851) %>%
  mutate(
    yr = NULL,
    Temperance_votes = NULL,
    Know_Nothing_votes = NULL,
    Republican_votes = NULL
  ) %>%
  left_join(eligible_pct, by = "combined") %>%
  mutate(
    ELIG = round(total / ELIG_1851_PCT),
    weight = ELIG / sum(ELIG),
    Democrat = Democrat_votes / ELIG,
    Whig = Whig_votes / ELIG,
    Free_Soil = Free_Soil_votes / ELIG,
    Abstaining = remainder(Democrat + Whig + Free_Soil)) %>%
  select(town, weight, Democrat, Whig, Free_Soil, Abstaining)

results.1852 <- raw_results %>%
  filter(yr == 1852) %>%
  mutate(
    yr = NULL,
    Temperance_votes = NULL,
    Know_Nothing_votes = NULL,
    Republican_votes = NULL
  ) %>%
  left_join(eligible_pct, by = "combined") %>%
  mutate(
    ELIG = round(total / ELIG_1852_PCT),
    weight = ELIG / sum(ELIG),
    Democrat = Democrat_votes / ELIG,
    Whig = Whig_votes / ELIG,
    Free_Soil = Free_Soil_votes / ELIG,
    Abstaining = remainder(Democrat + Whig + Free_Soil)) %>%
  select(town, weight, Democrat, Whig, Free_Soil, Abstaining)

results.1853 <- raw_results %>%
  filter(yr == 1853) %>%
  mutate(
    yr = NULL,
    Temperance_votes = NULL,
    Know_Nothing_votes = NULL,
    Republican_votes = NULL
  ) %>%
  left_join(eligible_pct, by = "combined") %>%
  mutate(
    ELIG = round(total / ELIG_1853_PCT),
    weight = ELIG / sum(ELIG),
    Democrat = Democrat_votes / ELIG,
    Whig = Whig_votes / ELIG,
    Free_Soil = Free_Soil_votes / ELIG,
    Abstaining = remainder(Democrat + Whig + Free_Soil)) %>%
  select(town, weight, Democrat, Whig, Free_Soil, Abstaining)

results.1854 <- raw_results %>%
  filter(yr == 1854) %>%
  mutate(
    yr = NULL,
    Know_Nothing_votes = NULL,
    Republican_votes = NULL
  ) %>%
  left_join(eligible_pct, by = "combined") %>%
  mutate(
    ELIG = round(total / ELIG_1854_PCT),
    weight = ELIG / sum(ELIG),
    Democrat = Democrat_votes / ELIG,
    Whig = Whig_votes / ELIG,
    Free_Soil = Free_Soil_votes / ELIG,
    Temperence = Temperance_votes / ELIG,
    Abstaining = remainder(Democrat + Whig + Free_Soil + Temperence)) %>%
  select(town, weight, Democrat, Whig, Free_Soil, Temperence, Abstaining)

results.1855 <- raw_results %>%
  filter(yr == 1855) %>%
  mutate(
    yr = NULL,
    Free_Soil_votes = NULL,
    Temperence_votes = NULL,
    Republican_votes = NULL
  ) %>%
  left_join(eligible_pct, by = "combined") %>%
  mutate(
    ELIG = round(total / ELIG_1855_PCT),
    weight = ELIG / sum(ELIG),
    Democrat = Democrat_votes / ELIG,
    Whig = Whig_votes / ELIG,
    Know_Nothing = Know_Nothing_votes / ELIG,
    Abstaining = remainder(Democrat + Whig + Know_Nothing)) %>%
  select(town, weight, Democrat, Whig, Know_Nothing, Abstaining)

results.1856 <- raw_results %>%
  filter(yr == 1856) %>%
  mutate(
    yr = NULL,
    Free_Soil_votes = NULL,
    Temperence_votes = NULL
  ) %>%
  left_join(eligible_pct, by = "combined") %>%
  mutate(
    ELIG = round(total / ELIG_1856_PCT),
    weight = ELIG / sum(ELIG),
    Democrat = Democrat_votes / ELIG,
    Whig = Whig_votes / ELIG,
    Know_Nothing = Know_Nothing_votes / ELIG,
    Republican = Republican_votes / ELIG,
    Abstaining = remainder(Democrat + Whig + Know_Nothing + Republican)) %>%
  select(town, weight, Democrat, Whig, Know_Nothing, Republican, Abstaining)

results.1857 <- raw_results %>%
  filter(yr == 1857) %>%
  mutate(
    yr = NULL,
    Whig_votes = NULL,
    Free_Soil_votes = NULL,
    Temperence_votes = NULL,
    Know_Nothing_votes = NULL
  ) %>%
  left_join(eligible_pct, by = "combined") %>%
  mutate(
    ELIG = round(total / ELIG_1857_PCT),
    weight = ELIG / sum(ELIG),
    Democrat = Democrat_votes / ELIG,
    Republican = Republican_votes / ELIG,
    Abstaining = remainder(Democrat + Republican)) %>%
  select(town, weight, Democrat, Republican, Abstaining)