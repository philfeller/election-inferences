# Utility functions used when creating election results;
# functions are called by prepare_results.R and covariate_regressions.R
# and are written to apply specifically to a limited range of years
#   create_results - prepares a tibble with percentage election results
#                    for each pair of years, combining towns as necessary
#   read_results   - reads raw election results from CSV and processes them
#   read_probate_results - reads raw election results for Judge of Probate
#                           and assigns party affiliations based on Governor's
#   create_factors - prepares a tibble with demographic factors for each town
#   combine_results - combines election results and dissolves towns in shapefile
#                     as necessary based on year range

# Define constants
source("./global.R")

# Either the 1855 or 1857 map will be used as the base shapefile
map.1855 <- sf::read_sf("./maps/1855_CT_towns.shp")
map.1857 <- sf::read_sf("./maps/1857_CT_towns.shp")

# Define functions to use when combining election results from multiple towns
combine_towns <- function(input_tibble, towns, combined_name) {
  # input_tibble: tibble with election results
  # towns: vector of town names to be combined
  # combined_name: name for the combined towns

  return(input_tibble %>%
    dplyr::filter(town %in% towns) %>%
    group_by(yr) %>%
    summarize(across(where(is.double), sum)) %>%
    add_column(town = combined_name, .after = "yr"))
}

# Define function to dissolve towns in a spatial dataframe
dissolve_towns <- function(map, town_vector, new_name, name_col = "town") {
  # map: sf object with town polygons
  # town_vector: vector of town names to be combined
  # new_name: name for the combined towns
  # name_col: column with town names in the sf object

  geom <- sf::st_union(map %>% dplyr::filter(.data[[name_col]] %in% town_vector))
  comb_poly <- sf::st_sf(setNames(data.frame(new_name), name_col), geometry = geom)
  rest <- map %>% dplyr::filter(!(.data[[name_col]] %in% town_vector))
  bind_rows(rest, comb_poly)
}

# Define function to get list of town groupings based on year range
year_groupings <- function(beg_yr, end_yr) {
  # beg_yr: beginning year of election transition
  # end_yr: ending year of election transition

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
    dplyr::filter(!(town %in% unlist(lapply(result_grp, function(x) x$towns))))

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

# Function to replace aliases in a vector with their preferred names;
# used to correct mistranscriptions in election data.
replace_aliases <- function(vec, mapping) {
  # mapping is a named vector where names are preferred names
  # and values are vectors of aliases to be replaced
  for (preferred in names(mapping)) {
    vec[vec %in% mapping] <- preferred
  }
  return(vec)
}

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

# Exclude towns that appear with zero votes in elections returns for some offices
# but do not appear in others for the same year. Having these towns appear causes
# results for different offices to have different numbers of rows, preventing linear
# regressions from being run.
exclude_towns <- c("Plainville", "East Granby", "Newington", "Ansonia", "Beacon Falls")

# Read and process the downloaded election results file
read_results <- function(file, alias_map, party_assignments, office, beg_yr, end_yr) {
  # file: path to CSV file with election results
  # alias_map: named vector for replacing candidate name aliases
  # party_assignments: tibble with candidate party assignments by year
  # office: office for which election results are being processed
  # beg_yr: beginning year of election returns to be processed
  # end_yr: ending year of election returns to be processed

  return(
    read_csv(file, show_col_types = FALSE) %>%
      select(
        election_date,
        office_name,
        candidate_name,
        granular_division_name,
        votes
      ) %>%
      rename(town = granular_division_name) %>%
      mutate(
        yr = as.integer(lubridate::year(election_date)), .before = candidate_name
      ) %>%
      select(-election_date) %>%
      # Remove total rows and filter by year and office
      dplyr::filter(
        candidate_name != "Total Ballots Cast",
        candidate_name != "Total Votes Cast",
        votes != 0,
        !town %in% exclude_towns,
        yr >= beg_yr,
        yr <= end_yr,
        office_name == office
      ) %>%
      # Map candidate names for 1855 Treasurer race for towns
      # that recorded a large number of votes for "All Other Votes";
      # transcriptions of the election ledger show that these were votes for
      # candidates on the right-hand page.
      # https://electionhistory.ct.gov/eng/contests/get_source_documentation/26870/1
      mutate(
        candidate_name = case_when(
          candidate_name == "All Other Votes" & yr == 1855 & office_name == "Treasurer" & town %in% c("Greenwich", "Sharon") ~ "Daniel W. Camp",
          candidate_name == "All Other Votes" & yr == 1855 & office_name == "Treasurer" & town %in% c("Darien", "Huntington", "Newtown", "Redding", "Plainsfield", "Harwinton") ~ "Arthur B. Calef",
          candidate_name == "All Other Votes" & yr == 1855 & office_name == "Treasurer" & town %in% c("Colebrook", "Norfolk") ~ "Amos Townsend, Jr.",
          candidate_name == "All Other Votes" & yr == 1855 & office_name == "Treasurer" & town %in% c("New Hartford") ~ "Talcott Crosby",
          candidate_name == "All Other Votes" & yr == 1855 & office_name == "Treasurer" & town %in% c("Tollund", "Bolton") ~ "Nehemiah D. Sperry",
          TRUE ~ candidate_name
        )
      ) %>%
      # Clean candidate names and join with party assignments
      mutate(candidate_name = replace_aliases(candidate_name, alias_map)) %>%
      dplyr::filter(candidate_name %in% party_assignments$candidate_name) %>%
      left_join(party_assignments, by = c("candidate_name", "office_name", "yr")) %>%
      # Remove rows with candidates not assigned to a party
      dplyr::filter(!is.na(candidate_party)) %>%
      # Change blank votes to zero
      replace(is.na(.), 0) %>%
      # Remove candidate names after joining
      select(-candidate_name) %>%
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

# Create results for Judge of Probate, assuming that a candidate's party
# affiliation is the same as for the Governor's race in the same year; use the
# top three Governor candidates in each town/year to assign party affiliations.
# Using the top two or four yielded similar results.
read_probate_results <- function(file, alias_map, party_assignments, beg_yr, end_yr) {
  # file: path to CSV file with election results
  # alias_map: named vector for replacing candidate name aliases
  # party_assignments: tibble with candidate party assignments by year
  # governor_results: tibble with processed election results for Governor
  # beg_yr: beginning year of election returns to be processed
  # end_yr: ending year of election returns to be processed

  probate_results <- read_csv(file, show_col_types = FALSE) %>%
    select(
      election_date,
      office_name,
      candidate_name,
      granular_division_name,
      district_name,
      votes
    ) %>%
    rename(town = granular_division_name) %>%
    mutate(
      yr = as.integer(lubridate::year(election_date)), .before = candidate_name
    ) %>%
    select(-election_date) %>%
    # Remove total rows and filter by year and office
    dplyr::filter(
      candidate_name != "Total Ballots Cast",
      candidate_name != "Total Votes Cast",
      votes != 0,
      !town %in% exclude_towns,
      yr >= beg_yr,
      yr <= end_yr,
      office_name == "Judge of Probate"
    ) %>%
    # Change blank votes to zero
    replace(is.na(.), 0)

  governor_party <- read_csv(file, show_col_types = FALSE) %>%
    select(
      election_date,
      office_name,
      candidate_name,
      granular_division_name,
      votes
    ) %>%
    rename(town = granular_division_name) %>%
    mutate(
      yr = as.integer(lubridate::year(election_date)), .before = candidate_name
    ) %>%
    select(-election_date) %>%
    # Remove total rows and filter by year and office
    dplyr::filter(
      candidate_name != "Total Ballots Cast",
      candidate_name != "Total Votes Cast",
      votes != 0,
      !town %in% exclude_towns,
      yr >= beg_yr,
      yr <= end_yr,
      office_name == "Governor"
    ) %>%
    # Clean candidate names and join with party assignments
    mutate(candidate_name = replace_aliases(candidate_name, alias_map)) %>%
    dplyr::filter(candidate_name %in% party_assignments$candidate_name) %>%
    left_join(party_assignments, by = c("candidate_name", "office_name", "yr")) %>%
    group_by(yr, town) %>%
    arrange(desc(votes)) %>%
    mutate(gov_rank = row_number()) %>%
    dplyr::filter(gov_rank <= 3) %>%
    select(yr, town, gov_rank, candidate_party)

  return(probate_results %>%
    group_by(yr, town) %>%
    arrange(desc(votes)) %>%
    mutate(judge_rank = row_number()) %>%
    dplyr::filter(judge_rank <= 3) %>%
    # Join the top two Governor candidates in same town/year
    left_join(governor_party, by = c("yr", "town", "judge_rank" = "gov_rank")) %>%
    # Remove candidate names after joining
    select(-c(candidate_name, judge_rank)) %>%
    # Create tibble that combines all town results into one row
    pivot_wider(
      names_from = candidate_party,
      values_from = votes,
      values_fn = sum,
      values_fill = 0
    ) %>%
    # Combine results for towns with more than one probate district
    group_by(yr, town) %>%
    summarise_if(is.double, sum, na.rm = TRUE))
}

# Create a tibble with demographic factors for a particular year range,
# appropriately combining towns that separated during the period
create_factors <- function(tibble, beg_yr, end_yr) {
  # tibble: input tibble with demographic factors
  # beg_yr: beginning year of election transition
  # end_yr: ending year of election transition

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
      age_1860 = ifelse(town %in% towns, age_1860, comb_age_1860)
    ) %>%
    select(-c(town, starts_with("comb_"))) %>%
    distinct(.keep_all = TRUE)
}

prepare_election_year <- function(results, shp, result_grp, map_grp, yr, beg_yr, end_yr) {
  # results: tibble with raw election results
  # shp: sf object with town polygons
  # result_grp: list of list(towns=..., new=...) for each combination
  # map_grp: list of list(towns=..., new=...) for each combination in the map
  # yr: year for which to prepare results
  # beg_yr: beginning year of election transition
  # end_yr: ending year of election transition

  # Get the parties for the year
  p_var <- get(paste("p", yr, sep = ""))

  # Extract party names, removing the "_in_185x" suffix, and dropping "Abstaining"
  parties_to_keep <- str_remove_all(p_var, "_in_185[0-9]")[-1 * length(p_var)]

  party_vote_cols <- paste0(parties_to_keep, "_votes")
  year <- paste("18", yr, sep = "")

  results %>%
    dplyr::filter(yr == year) %>%
    combine_results(., shp, result_grp, map_grp) %>%
    magrittr::extract2("results") %>%
    dplyr::filter(yr >= beg_yr & yr <= end_yr) %>%
    select(any_of(c("town", party_vote_cols, "total")), -yr) %>%
    rename_with(
      ~ str_replace(.x, "_votes$", paste0("_vote_in_18", yr)),
      ends_with("_votes")
    ) %>%
    rename(!!paste0("total_18", yr) := total)
}

# Create a results tibble that is appropriate for the years to be compared;
# don't combine election results for towns unless necessary.
# This is particularly important in Windham County, where Killingly, Pomfret,
# and Thompson were split to form Putnam. Unless comparing to 1856 or 1857
# results, combining these towns results in a significant loss of detail.
create_results <- function(results, shp, beg_yr, end_yr, eligible_pct, factors) {
  # results: tibble with raw election results
  # shp: sf object with town polygons
  # beg_yr: beginning year of election transition
  # end_yr: ending year of election transition
  # eligible_pct: tibble with percentage of eligible voters who voted in each town combination
  # factors: tibble with demographic factors for each town

  # Create list of towns to be combined, based on year range
  result_grp <- year_groupings(beg_yr, end_yr)
  # Create list of towns to be combined in a shapefile, based on year range
  # Because of Putnam's formation in 1856, the base map will be that of 1855
  # or 1857, depending on the end year
  map_grp <- year_groupings(beg_yr, max(1855, end_yr))

  # Choose the appropriate shapefile from which to get geographic data
  shp <- combine_results(results, shp, result_grp, map_grp) %>%
    magrittr::extract2("shp") %>%
    arrange(town)

  # Create spatial weights for the towns
  nb <- spdep::poly2nb(shp)
  listw <- spdep::nb2listw(nb, style = "W")

  # Get geographic centroids for each town
  centroids <- sf::st_centroid(shp)
  geo <- as.data.frame(sf::st_coordinates(centroids)) %>%
    rename(lon = X) %>%
    rename(lat = Y) %>%
    bind_cols(shp %>% sf::st_set_geometry(NULL) %>% select(town))
  
  # Transform centroids to use same the CRS as the railroads shapefile
  centroids <- sf::st_transform(centroids, sf::st_crs(railroads))
  
  # Calculate distance to nearest railroad for each town centroid
  rr_dist <- as.data.frame(sf::st_distance(centroids, railroads)[,2:9] %>%
    apply(1, min) %>%
    as.numeric()
  )
  
  colnames(rr_dist) <- c("rr_dist")
  geo <- bind_cols(geo, rr_dist)

  # Generate the demographic factors appropriate for the range of years
  demo_factors <- factors %>%
    create_factors(beg_yr, end_yr) %>%
    rename(town = combined)

  # Create a tibble that combines the election data for each year in the period,
  # keeping only the data for the beginning and end years.
  yrs <- 51:57
  combined <- yrs %>%
    map(~ prepare_election_year(raw_results, shp, result_grp, map_grp, .x, beg_yr, end_yr)) %>%
    reduce(full_join, by = "town")

  # Create a tibble with percentage results, including estimated nonvoters.
  full_results <- combined %>%
    dplyr::filter(!if_all(everything(), is.na)) %>%
    left_join(geo, by = "town") %>%
    mutate(combined = get_combined(town)) %>%
    left_join(eligible_pct, by = "combined") %>%
    left_join(demo_factors, by = "town") %>%
    arrange(town) %>%
    # Create columns for estimated eligible voters and non-voting eligible voters
    mutate(across(
      starts_with("total_18"),
      ~ round(.x / get(paste0("ELIG_", str_sub(cur_column(), -4), "_PCT"))),
      .names = "ELIG_{str_extract(.col, '(\\\\d{4})')}"
    )) %>%
    mutate(across(
      starts_with("ELIG_18") & ! ends_with("_PCT"),
      ~ .x - get(paste0("total_", str_sub(cur_column(), -4))),
      .names = "nonvote_in_{str_extract(.col, '(\\\\d{4})')}"
    )) %>%
    mutate(
      Democrat_in_1851 = Democrat_vote_in_1851 / ELIG_1851,
      Whig_in_1851 = Whig_vote_in_1851 / ELIG_1851,
      Free_Soil_in_1851 = Free_Soil_vote_in_1851 / ELIG_1851,
      Abstaining_in_1851 = remainder(Democrat_in_1851 + Whig_in_1851 + Free_Soil_in_1851),
      Democrat_in_1852 = Democrat_vote_in_1852 / ELIG_1852,
      Whig_in_1852 = Whig_vote_in_1852 / ELIG_1852,
      Free_Soil_in_1852 = Free_Soil_vote_in_1852 / ELIG_1852,
      Abstaining_in_1852 = remainder(Democrat_in_1852 + Whig_in_1852 + Free_Soil_in_1852),
      Democrat_in_1853 = Democrat_vote_in_1853 / ELIG_1853,
      Whig_in_1853 = Whig_vote_in_1853 / ELIG_1853,
      Free_Soil_in_1853 = Free_Soil_vote_in_1853 / ELIG_1853,
      Abstaining_in_1853 = remainder(Democrat_in_1853 + Whig_in_1853 + Free_Soil_in_1853),
      Democrat_in_1854 = Democrat_vote_in_1854 / ELIG_1854,
      Whig_in_1854 = Whig_vote_in_1854 / ELIG_1854,
      Free_Soil_in_1854 = Free_Soil_vote_in_1854 / ELIG_1854,
      Temperance_in_1854 = Temperance_vote_in_1854 / ELIG_1854,
      Abstaining_in_1854 = remainder(Democrat_in_1854 + Whig_in_1854 + Free_Soil_in_1854 + Temperance_in_1854),
      Democrat_in_1855 = Democrat_vote_in_1855 / ELIG_1855,
      Whig_in_1855 = Whig_vote_in_1855 / ELIG_1855,
      Know_Nothing_in_1855 = Know_Nothing_vote_in_1855 / ELIG_1855,
      Abstaining_in_1855 = remainder(Democrat_in_1855 + Whig_in_1855 + Know_Nothing_in_1855),
      Democrat_in_1856 = Democrat_vote_in_1856 / ELIG_1856,
      Republican_in_1856 = Republican_vote_in_1856 / ELIG_1856,
      Know_Nothing_in_1856 = Know_Nothing_vote_in_1856 / ELIG_1856,
      Whig_in_1856 = Whig_vote_in_1856 / ELIG_1856,
      Abstaining_in_1856 = remainder(Democrat_in_1856 + Republican_in_1856 + Know_Nothing_in_1856 + Whig_in_1856),
      Democrat_in_1857 = Democrat_vote_in_1857 / ELIG_1857,
      Republican_in_1857 = Republican_vote_in_1857 / ELIG_1857,
      Abstaining_in_1857 = remainder(Democrat_in_1857 + Republican_in_1857)
    ) %>%
    select_if(function(x) !any(is.na(x)))

  # Calculate spatial lag variables for party strength in beginning and ending years;
  # spatial lag captures the extent to which results are affected by those in
  # nearby towns, weighted by adjacency.
  b_name <- paste0("p", substr(beg_yr, 3, 4))
  e_name <- paste0("p", substr(end_yr, 3, 4))
  vars <- get(b_name) # e.g., c("Democrat_in_1854", ...)
  vars <- c(vars, get(e_name))

  for (v in vars) {
    lag_name <- paste0("lag_", v)
    lag_var <- spdep::lag.listw(listw, full_results[[v]])
    full_results[[lag_name]] <- lag_var
  }

  return(list(results = full_results, shp = shp))
}

# Define function to calculate the weighted standard deviation
weighted.sd <- function(results) {
  # results: tibble with election results, including weighted mean row

  returns <- as.matrix(results %>% dplyr::filter(town != "Weighted mean") %>% select(Democrat:Abstaining))
  num_rows <- nrow(returns)
  mu <- unlist(rep((results %>% dplyr::filter(town == "Weighted mean") %>% select(Democrat:Abstaining)), num_rows, byrow = TRUE))
  num_cols <- length(mu) / num_rows
  mu <- t(matrix(mu, nrow = num_cols))
  weights <- unlist(results %>% dplyr::filter(town != "Weighted mean") %>% select(weight))
  sd <- as.data.frame(t(sqrt(colSums((returns - mu)^2 * weights))))
  cbind(data.frame(town = "Weighted SD", weight = NA), sd)
}

# Define functions to be used when generating yearly results
yr_results <- function(raw, filter_yr) {
  # raw: tibble with raw election results
  # filter_yr: year for which results are to be extracted

  raw %>%
    dplyr::filter(yr == filter_yr) %>%
    select_if(function(x) any(x > 0)) %>%
    select(-yr) %>%
    left_join(eligible_pct, by = "combined")
}

# Define function to generate summary results with weighted mean and SD
result_summary <- function(results) {
  # results: tibble with election results for a particular year

  results_plus_mean <- results %>%
    select(!starts_with("ELIG") & !ends_with("votes") & !total & !combined) %>%
    bind_rows(summarise(., across(Democrat:Abstaining, ~ weighted.mean(.x, weight)))) %>%
    mutate(town = replace(town, is.na(town), "Weighted mean"))
  result_sd <- weighted.sd(results_plus_mean)
  results_plus_mean %>% bind_rows(result_sd)
}
