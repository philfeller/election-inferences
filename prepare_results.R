# Define a function that will be used to create a tibble with gubernatorial
# results for a pair of years:
#   create_results - prepares a tibble with percentage election results

library(magrittr)
library(sf)
library(spdep)

# Define constants
source("./variables.R")

# Estimate nonvoters and calculate potential covariates
source("./ct_demographics.R")

# Either the 1855 or 1857 map will be used as the base shapefile
map.1855 <- read_sf("./maps/1855_CT_towns.shp")
map.1857 <- read_sf("./maps/1857_CT_towns.shp")

# Define functions to use when combining election results from multiple towns
combine_towns <- function(input_tibble, towns, combined_name) {
  return(input_tibble %>%
    filter(town %in% towns) %>%
    group_by(yr) %>%
    summarize(across(where(is.double), sum)) %>%
    add_column(town = combined_name, .after = "yr"))
}

# Define function to dissolve towns in a spatial dataframe
dissolve_towns <- function(map, town_vector, new_name, name_col = "town") {
  geom <- st_union(map %>% filter(.data[[name_col]] %in% town_vector))
  comb_poly <- st_sf(setNames(data.frame(new_name), name_col), geometry = geom)
  rest <- map %>% filter(!(.data[[name_col]] %in% town_vector))
  bind_rows(rest, comb_poly)
}

# Define function to get list of town groupings based on year range
year_groupings <- function(beg_yr, end_yr) {
  groupings <- list()

  # When comparing pre-1856 results to 1856 or 1857, combine towns that split in 1856
  if (beg_yr < 1856 & end_yr >= 1856) {
    # Putnam splits out in 1856; combine back for earlier
    groupings <- c(groupings, list(list(towns = killingly_et_al, new = "Killingly, Thompson, and Pomfret")))
    # Bethel splits out in 1856; combine back for earlier
    groupings <- c(groupings, list(list(towns = danbury, new = "Danbury")))
    # South Lyme splits out in 1856; combine back for earlier
    groupings <- c(groupings, list(list(towns = lyme, new = "Lyme")))
  }

  # When comparing pre-1855 results to 1855 or later, combine towns that split in 1855
  if (beg_yr < 1855 & end_yr >= 1855) {
    # West Hartford splits out in 1855; combine back for earlier
    groupings <- c(groupings, list(list(towns = hartford, new = "Hartford")))
    # Windsor Locks splits out in 1855; combine back for earlier
    groupings <- c(groupings, list(list(towns = windsor, new = "Windsor")))
  }

  # When comparing pre-1852 results to 1852 or later, combine town that split in 1851
  if (beg_yr < 1852 & end_yr >= 1852) {
    # Cromwell splits out in 1852; combine back for earlier
    groupings <- c(groupings, list(list(towns = middletown, new = "Middletown")))
  }

  # Division and renaming of parts of Saybrook affect results beginning in 1853
  groupings <- c(groupings, list(list(towns = saybrook, new = "Saybrook")))

  # Bridgewater , taken from New Milford, is in 1857 results
  groupings <- c(groupings, list(list(towns = new_milford, new = "New Milford")))

  return(groupings)
}

combine_results <- function(results, shp, result_grp, map_grp, name_col = "town") {
  # results: the tibble/data.frame to combine
  # shp: the sf object to dissolve
  # result_grp: list of list(towns=..., new=...) for each combination
  # map_grp: list of list(towns=..., new=...) for each combination in the map
  # name_col: column with town names in both results and shp

  # Save towns not to be combined
  rest <- results %>%
    filter(!(town %in% unlist(lapply(result_grp, function(x) x$towns))))

  # Combine towns in results and dissolve polygons in shp for each grouping
  for (grp in result_grp) {
    # Combine results for town grouping
    combined <- combine_towns(results, grp$towns, grp$new)
    # Add to results
    rest <- bind_rows(rest, combined)
  }
  for (grp in map_grp) {
    # Dissolve in shapefile
    shp <- dissolve_towns(shp, grp$towns, grp$new, name_col = name_col)
  }

  # Make sure the column order and names are consistent
  results <- rest %>% arrange(town)
  shp <- shp %>% arrange(town)

  # Return both
  list(results = results, shp = shp)
}

# Define function to assign party designation to candidates
assign_party <- function(input_tibble, other = FALSE) {
  # Nominees of major parties between 1849 and 1857
  major_candidates <- c(
    "John M. Niles",
    "Joseph Trumbull",
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

  # Alternate names used to vote for Niles, Seymour, Ingham, and Foster
  niles <- c("J. M. Niles")
  seymour <- c("Jonas H. Senor", "T. H. Seymour", "T. Seymour", "Thomas Seymour", "Thomas S. Seymour", "Ths. H. Seymour", "Tom H. Seymour")
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
    mutate(candidate_name = ifelse(candidate_name %in% niles, "John M. Niles", candidate_name)) %>%
    mutate(candidate_name = ifelse(candidate_name %in% seymour, "Thomas H. Seymour", candidate_name)) %>%
    mutate(candidate_name = ifelse(candidate_name %in% foster, "Lafayette S. Foster", candidate_name)) %>%
    mutate(candidate_name = ifelse(candidate_name %in% ingham, "Samuel Ingham", candidate_name)) %>%
    filter(candidate_name %in% major_candidates) %>%
    # Remove candidates whose names appeared on stray ballots when their not party's nominee
    # Francis Gillette is coded as Free Soil in 1849, even though John M. Niles was the party
    # nominee; he was the Liberty party candidate in 1848, before the Free Soil party formed,
    # and he outpolled Niles 88-4 in East Haven
    filter(!(yr == 1856 && candidate_name == "Henry Dutton")) %>%
    filter(!(yr > 1853 && candidate_name == "Francis Gillette")) %>%
    filter(!(yr == 1853 && candidate_name == "Samuel Ingham")) %>%
    filter(!(yr == 1857 && candidate_name == "William T. Minor")) %>%
    mutate(candidate_party = case_when(
      candidate_name == "John M. Niles" ~ "Free_Soil_votes",
      candidate_name == "Joseph Trumbull" ~ "Whig_votes",
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
      filter(yr >= 1849) %>%
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

# Create a results tibble that is appropriate for the years to be compared;
# don't combine election results for towns unless necessary.
# This is particularly important in Windham County, where Killingly, Pomfret,
# and Thompson were split to form Putnam. Unless comparing to 1856 or 1857
# results, combining these towns results in a significant loss of detail.
create_results <- function(results, shp, beg_yr, end_yr) {
  # Create list of towns to be combined based on year range
  result_grp <- year_groupings(beg_yr, end_yr)
  # Create list of towns to be combined in shapefile based on year range
  # Because of Putnam's formation in 1856, the base map will be that of 1855
  # or 1857, depending on the end year
  map_grp <- year_groupings(beg_yr, max(1855, end_yr))
  # Create separate tibbles for each year
  e51 <- results %>%
    filter(yr == 1851) %>%
    combine_results(., shp, result_grp, map_grp) %>%
    extract2("results") %>%
    filter(yr >= beg_yr & yr <= end_yr) %>%
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
  e52 <- results %>%
    filter(yr == 1852) %>%
    combine_results(., shp, result_grp, map_grp) %>%
    extract2("results") %>%
    filter(yr >= beg_yr & yr <= end_yr) %>%
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

  e53 <- results %>%
    filter(yr == 1853) %>%
    combine_results(., shp, result_grp, map_grp) %>%
    extract2("results") %>%
    filter(yr >= beg_yr & yr <= end_yr) %>%
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
  e54 <- results %>%
    filter(yr == 1854) %>%
    combine_results(., shp, result_grp, map_grp) %>%
    extract2("results") %>%
    filter(yr >= beg_yr & yr <= end_yr) %>%
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
  e55 <- results %>%
    filter(yr == 1855) %>%
    combine_results(., shp, result_grp, map_grp) %>%
    extract2("results") %>%
    filter(yr >= beg_yr & yr <= end_yr) %>%
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
  e56 <- results %>%
    filter(yr == 1856) %>%
    combine_results(., shp, result_grp, map_grp) %>%
    extract2("results") %>%
    filter(yr >= beg_yr & yr <= end_yr) %>%
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
  e57 <- results %>%
    filter(yr == 1857) %>%
    combine_results(., shp, result_grp, map_grp) %>%
    extract2("results") %>%
    filter(yr >= beg_yr & yr <= end_yr) %>%
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

  # Choose the appropriate shapefile from which to get geographic data
  shp <- combine_results(results, shp, result_grp, map_grp) %>%
    extract2("shp") %>%
    arrange(town)

  # Create spatial weights for the towns
  nb <- poly2nb(shp)
  listw <- nb2listw(nb, style = "W")

  # Get geographic centroids for each town
  geo <- as.data.frame(st_coordinates(st_centroid(shp))) %>%
    rename(lon = X) %>%
    rename(lat = Y) %>%
    bind_cols(shp %>% st_set_geometry(NULL) %>% select(town))

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
    left_join(geo, by = "town") %>%
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

  # Calculate spatial lag variables for party strength in beginning year
  p_name <- paste0("p", substr(beg_yr, 3, 4))
  vars <- get(p_name) # e.g., c("Democrat_in_1854", ...)

  for (v in vars) {
    lag_name <- paste0("lag_", v)
    lag_var <- lag.listw(listw, full_results[[v]])
    full_results[[lag_name]] <- lag_var
  }

  return(list(results = full_results, shp = shp))
}

weighted.sd <- function(results) {
  returns <- as.matrix(results %>% filter(town != "Weighted mean") %>% select(Democrat:Abstaining))
  num_rows <- nrow(returns)
  mu <- unlist(rep((results %>% filter(town == "Weighted mean") %>% select(Democrat:Abstaining)), num_rows, byrow = TRUE))
  num_cols <- length(mu) / num_rows
  mu <- t(matrix(mu, nrow = num_cols))
  weights <- unlist(results %>% filter(town != "Weighted mean") %>% select(weight))
  sd <- as.data.frame(t(sqrt(colSums((returns - mu)^2 * weights))))
  cbind(data.frame(town = "Weighted SD", weight = NA), sd)
}

# Define functions to be used when generating yearly results
yr_results <- function(raw, filter_yr) {
  raw %>%
    filter(yr == filter_yr) %>%
    select_if(function(x) any(x > 0)) %>%
    mutate(yr = NULL) %>%
    left_join(eligible_pct, by = "combined")
}

result_summary <- function(results) {
  results_plus_mean <- results %>%
    select(!starts_with("ELIG") & !ends_with("votes") & !total & !combined) %>%
    bind_rows(summarise(., across(Democrat:Abstaining, ~ weighted.mean(.x, weight)))) %>%
    mutate(town = replace(town, is.na(town), "Weighted mean"))
  result_sd <- weighted.sd(results_plus_mean)
  results_plus_mean %>% bind_rows(result_sd)
}

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
