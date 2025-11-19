# Utility functions used to process IPUMs census data
# - get_1850_town(SERIAL): Return the town name for a given SERIAL number from the 1850 census
# - get_1860_town(SERIAL): Return the town name for a given SERIAL number from the 1860 census
# - get_combined(town): Function to assign combined name, the name of two towns before splitting
# - estimate_voters(yr, native_1850, foreign_1850, native_change, poll_change):
#   Estimate the number of eligible voters between census years; assume a linear change for native-born


# Return the town name for a given SERIAL number from the 1850 census
get_1850_town <- function(SERIAL) {
  # SERIAL: IPUMS SERIAL number
  
  as.character(town_1850[SERIAL])
}

# Return the town name for a given SERIAL number from the 1860 census
get_1860_town <- function(SERIAL) {
  # SERIAL: IPUMS SERIAL number
  
  as.character(town_1860[SERIAL])
}

# Function to assign combined name, the name of two towns before splitting;
# combined name is used to match 1860 census towns to towns in the 1850 census
# and in the voting record and to calculate demographic data that depends only
# on one census
get_combined <- function(town) {
  # town: name of town
  
  case_when(
    town %in% canaan ~ "Canaan",
    town %in% litchfield ~ "Litchfield",
    town %in% windham ~ "Windham",
    town %in% granby ~ "Granby",
    town %in% new_milford ~ "New Milford",
    town %in% killingly_et_al ~ "Killingly, Thompson, and Pomfret",
    town %in% danbury ~ "Danbury",
    town %in% lyme ~ "Lyme",
    town %in% hartford ~ "Hartford",
    town %in% windsor ~ "Windsor",
    town %in% saybrook ~ "Saybrook",
    town %in% middletown ~ "Middletown",
    TRUE ~ town
  )
}

# Estimate the number of eligible voters between census years; assume a
# linear change for native-born, that immigration numbers follow the same
# trend as for the nation as a whole, and that adult male immigrants make
# up a consistent portion of all immigrants from year to year.
# Nationwide immigration by year is taken from the 2022 DHS Yearbook:
# https://www.dhs.gov/ohss/topics/immigration/yearbook/2022
estimate_voters <- function(yr, native_1850, foreign_1850, native_change, poll_change) {
  # yr: year for which to estimate eligible voters
  # native_1850: number of native-born, white, adult males in 1850 census
  # foreign_1850: number of foreign-born, white, adult males in 1850 census
  # native_change: estimated annual change in native-born, white, adult males
  # poll_change: estimated annual change in taxable polls (used for towns with
  #               incorrect birthplace transcriptions)
  
  # For towns with incorrect birthplace transcriptions estimate the total
  # number of eligible using the increase in taxable polls
  native_voters <- native_1850 + round((yr - 1850) * ifelse(poll_change == 0, native_change, poll_change))
  
  cum_1820_1845 <- 8385 + 9127 + 6911 + 6354 + 7912 + 10199 + 10837 + 18875 +
    27382 + 22520 + 23322 + 22633 + 60482 + 58640 + 65365 + 45374 + 76242 +
    79340 + 38914 + 68069 + 84066 + 80289 + 104565 + 52496 + 78615 + 114371
  imm_1846 <- 154416
  imm_1847 <- 234968
  imm_1848 <- 226527
  imm_1849 <- 297024
  imm_1850 <- 369980
  cum_1820_1850 <- cum_1820_1845 + imm_1846 + imm_1847 + imm_1848 + imm_1849 + imm_1850
  # Calculate the percentage of foreign-born people in the 1850 census who had
  # immigrated between 1820 and 1845. DHS immigration numbers are for the fiscal
  # year ending June 30, and there is a five-year lag before new immigrants meet
  # the Connecticut residency requirement, making those who immigrated before
  # June 30, 1845, eligible for the 1851 election.
  pct_1851 <- cum_1820_1845 / cum_1820_1850
  pct_1852 <- (cum_1820_1845 + imm_1846) / cum_1820_1850
  pct_1853 <- (cum_1820_1845 + imm_1846 + imm_1847) / cum_1820_1850
  pct_1854 <- (cum_1820_1845 + imm_1846 + imm_1847 + imm_1848) / cum_1820_1850
  pct_1855 <- (cum_1820_1845 + imm_1846 + imm_1847 + imm_1848 + imm_1849) / cum_1820_1850
  #  Assume that the literacy requirement enacted in 1855 effectively suppresses
  # significant addition of foreign-born voters for the 1856 and 1857 elections.
  foreign_voters <- case_when(
    yr == 1851 ~ round(foreign_1850 * pct_1851),
    yr == 1852 ~ round(foreign_1850 * pct_1852),
    yr == 1853 ~ round(foreign_1850 * pct_1853),
    yr == 1854 ~ round(foreign_1850 * pct_1854),
    yr >= 1855 ~ round(foreign_1850 * pct_1855)
  )
  
  total_voters <- native_voters + foreign_voters
  
  return(total_voters)
}