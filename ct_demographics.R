# Extract demographic data from IPUMS full-count census records for Connecticut.
# Variables and function defined by this script are used in other scripts:
#   get_combined - function that assigns a name to a town that changes shape
#   ct_eligible - tibble with estimated number of eligible voters by town and year
#   factors - tibble with calculated informtion about a town, such as gini

library(dplyr, warn.conflicts = FALSE)
library(readr)
library(tibble)
library(tidyr)
library(lubridate, warn.conflicts = FALSE)
library(ipumsr)
library(ineq)

source("./variables.R")

get_1850_town <- function(SERIAL) {
  as.character(town_1850[SERIAL])
}

get_1860_town <- function(SERIAL) {
  as.character(town_1860[SERIAL])
}

# Function to assign combined name, the name of two towns before splitting;
# combined name is used to match 1860 census towns to towns in the 1850 census
# and in the voting record and to calculate demographic data that depends only
# on one census
get_combined <- function(town) {
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

ddi1850 <- read_ipums_ddi(ipums_1850)
ddi1860 <- read_ipums_ddi(ipums_1860)
ct_1850 <- read_ipums_micro(ddi1850, verbose = FALSE) %>%
  select(HIK, SERIAL, GQ, FAMUNIT, RELATE, SEX, AGE, RACE, BPL,
         OCC1950, REALPROP, SCHOOL, LIT, PAUPER, CRIME) %>%
  mutate(
    family = SERIAL * 100 + FAMUNIT,
    BIRTH = ifelse(BPL < 100, "native", "immigrant"),
    JOB = ifelse(OCC1950 %in% c(810, 820, 830, 840, 100, 123), "farm", "nonfarm"),
    AGE_CAT = case_when(
      AGE < 20 ~ "0 - 19",
      AGE >= 20 & AGE < 30 ~ "20 - 29",
      AGE >= 30 & AGE < 40 ~ "30 - 39",
      AGE >= 40 & AGE < 50 ~ "40 - 49",
      AGE >= 50 & AGE < 60 ~ "50 - 59",
      AGE >= 60 ~ "60 and over"
    ),
    # The CLASS column categorizes people according to Doherty's 1977 model
    CLASS = case_when(
      SEX == 1 & AGE > 30 & REALPROP >= 10000 ~ "elite",
      SEX == 1 & AGE > 30 & REALPROP >= 750 & REALPROP < 10000 ~ "middle",
      SEX == 1 & AGE >= 16 & AGE <= 30 & REALPROP < 750 ~ "young",
      SEX == 1 & AGE > 30 & REALPROP < 750 ~ "casualties"
    ),
    town = get_1850_town(SERIAL),
    combined = get_combined(town)
  ) %>%
  filter(town != "NULL") %>%
  arrange(SERIAL, FAMUNIT, RELATE)

# Load file with data that is missing from IPUMS: first pages for Brooklyn and Hebron
load(file = "missing_1860_ipums_rows.Rda")

# The IPUMS data contains some transcription errors. Although GINI indices calculated
# from IPUMS data vary insignificantly from indices calculated with correct real and
# personal property values, update transcription errors for Meriden with corrected
# values.
ct_1860 <- read_ipums_micro(ddi1860, verbose = FALSE) %>%
  mutate(RELATE = ifelse(SERIAL == 294445, 1, RELATE)) %>%
  select(HIK, SERIAL, GQ, FAMUNIT, RELATE, SEX, AGE, RACE, BPL,
         OCC1950, REALPROP, SCHOOL, LIT, PAUPER, CRIME, PERSPROP) %>%
  bind_rows(missing_rows) %>%
  left_join(read_csv("corrected_meriden_wealth.csv", show_col_types = FALSE),
    by = c("SERIAL", "SEX", "AGE", "BPL", "REALPROP", "PERSPROP")
  ) %>%
  mutate(
    REALPROP = ifelse(!is.na(real) & is.numeric(real), real, REALPROP),
    PERSPROP = ifelse(!is.na(pers) & is.numeric(pers), pers, PERSPROP)
  ) %>%
  mutate(
    family = SERIAL * 100 + FAMUNIT,
    WEALTH = REALPROP + PERSPROP,
    BIRTH = ifelse(BPL < 100, "native", "immigrant"),
    JOB = ifelse(OCC1950 %in% c(810, 820, 830, 840, 100, 123), "farm", "nonfarm"),
    AGE_CAT = case_when(
      AGE < 20 ~ "0 - 19",
      AGE >= 20 & AGE < 30 ~ "20 - 29",
      AGE >= 30 & AGE < 40 ~ "30 - 39",
      AGE >= 40 & AGE < 50 ~ "40 - 49",
      AGE >= 50 & AGE < 60 ~ "50 - 59",
      AGE >= 60 ~ "60 and over"
    ),
    # The CLASS column categorizes people according to Doherty's 1977 model
    CLASS = case_when(
      SEX == 1 & AGE > 30 & REALPROP >= 10000 ~ "elite",
      SEX == 1 & AGE > 30 & REALPROP >= 750 & REALPROP < 10000 ~ "middle",
      SEX == 1 & AGE >= 16 & AGE <= 30 & REALPROP < 750 ~ "young",
      SEX == 1 & AGE > 30 & REALPROP < 750 ~ "casualties"
    ),
    town = get_1860_town(SERIAL),
    combined = get_combined(town)
  ) %>%
  filter(town != "NULL") %>%
  arrange(SERIAL, FAMUNIT, RELATE)

# Load 1860 religious-accommodation data, which is used to estimate degree of denominational affiliation.
# Data were hand-entered into a spreadsheet from FamilySearch Social Statistics census schedule images:
# https://www.familysearch.org/records/images/search-results?page=1&place=346&endDate=1860&startDate=1860&creator=Federal%20Census
religion_1860 <- read_csv(religion_file, show_col_types = FALSE) %>%
  select(Town, where(is.numeric)) %>%
  mutate(combined = get_combined(Town)) %>%
  rowwise() %>%
  mutate(total = sum(c_across(Congregational:Friends))) %>%
  group_by(combined) %>%
  summarise(
    cong = sum(Congregational),
    bap = sum(Baptist),
    meth = sum(Methodist),
    epis = sum(Episcopal),
    piet = sum(Baptist) + sum(Methodist) + sum(Christian) + sum(Disciples) + sum(Sandemanian) + sum(Free) + sum(Union) + sum(Adventist) + sum(Spiritualist) + sum(Friends),
    total = sum(total)
  ) %>%
  mutate(
    pct_cong = cong / total,
    pct_bap = bap / total,
    pct_meth = meth / total,
    pct_epis = epis / total,
    pct_piet = piet / total
  ) %>%
  select(combined, starts_with("pct_"))

# Construct tibbles with data about household wealth and demographics
ct_1850_hh <- ct_1850 %>%
  # Exclude servants and institutional housing
  filter((GQ %in% c(1, 2, 5) && FAMUNIT == 1) || GQ == 4) %>%
  mutate(family = SERIAL * 100 + FAMUNIT) %>%
  group_by(family) %>%
  summarise(
    FAMILY_REALPROP = sum(REALPROP),
    FAMILY_HEAD_AGE = first(AGE),
    AGE_CAT = first(AGE_CAT),
    BIRTH = first(BIRTH),
    JOB = first(JOB),
    CLASS = first(CLASS),
    town = first(town),
    combined = first(combined)
  )

ct_1860_hh <- ct_1860 %>%
  # Exclude servants and institutional housing
  filter((GQ %in% c(1, 2, 5) && FAMUNIT == 1) || GQ == 4) %>%
  mutate(family = SERIAL * 100 + FAMUNIT) %>%
  group_by(family) %>%
  summarise(
    FAMILY_REALPROP = sum(REALPROP),
    FAMILY_WEALTH = sum(WEALTH),
    FAMILY_HEAD_AGE = first(AGE),
    AGE_CAT = first(AGE_CAT),
    BIRTH = first(BIRTH),
    JOB = first(JOB),
    CLASS = first(CLASS),
    town = first(town),
    combined = first(combined)
  )

white_male_1860_hh <- ct_1860 %>%
  # Exclude servants and institutional housing
  filter((GQ %in% c(1, 2, 5) && FAMUNIT == 1) || GQ == 4) %>%
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
  filter(SEX == 1 && RACE == 1)

ct_pop <- ct_1850 %>%
  group_by(combined) %>%
  summarize(POP_1850 = n()) %>%
  inner_join(ct_1860 %>%
    group_by(combined) %>%
    summarize(POP_1860 = n()), by = "combined")

# Count number of native-born, white, adult males to estimate eligible voters;
# the number of naturalized foreign-born citizens who meet the residency
# requirement will be estimated in the estimate_voters function.
ct_eligible <- ct_1850 %>%
  filter(SEX == 1) %>%
  filter(BPL < 100) %>%
  filter(RACE == 1) %>%
  filter(AGE > 20) %>%
  group_by(combined) %>%
  summarize(ELIG_1850 = n()) %>%
  inner_join(ct_1850 %>%
    filter(SEX == 1) %>%
    filter(BPL >= 100) %>%
    filter(RACE == 1) %>%
    filter(AGE > 20) %>%
    group_by(combined) %>%
    summarize(FOREIGN_1850 = n()), by = "combined") %>%
  inner_join(ct_1860 %>%
    filter(SEX == 1) %>%
    filter(BPL < 100) %>%
    filter(RACE == 1) %>%
    filter(AGE > 20) %>%
    group_by(combined) %>%
    summarize(ELIG_1860 = n()), by = "combined") %>%
  mutate(POP_CHANGE = (ELIG_1860 - ELIG_1850) / 10) %>%
  mutate(POLL_CHANGE = case_when(
    combined == "Avon" ~ avon_poll_change,
    combined == "Burlington" ~ burlington_poll_change,
    combined == "Farmington" ~ farmington_poll_change,
    combined == "New Britain" ~ new_britain_poll_change,
    !combined %in% bad_birthplace ~ 0
  )) %>%
  mutate(ELIG_1851 = estimate_voters(1851, ELIG_1850, FOREIGN_1850, POP_CHANGE, POLL_CHANGE)) %>%
  mutate(ELIG_1852 = estimate_voters(1852, ELIG_1850, FOREIGN_1850, POP_CHANGE, POLL_CHANGE)) %>%
  mutate(ELIG_1853 = estimate_voters(1853, ELIG_1850, FOREIGN_1850, POP_CHANGE, POLL_CHANGE)) %>%
  mutate(ELIG_1854 = estimate_voters(1854, ELIG_1850, FOREIGN_1850, POP_CHANGE, POLL_CHANGE)) %>%
  mutate(ELIG_1855 = estimate_voters(1855, ELIG_1850, FOREIGN_1850, POP_CHANGE, POLL_CHANGE)) %>%
  mutate(ELIG_1856 = estimate_voters(1856, ELIG_1850, FOREIGN_1850, POP_CHANGE, POLL_CHANGE)) %>%
  mutate(ELIG_1857 = estimate_voters(1857, ELIG_1850, FOREIGN_1850, POP_CHANGE, POLL_CHANGE))

# Find information about people who moved out of state between 1850 and 1860.
# The 1860 records need to have an HIK variable in both censuses.
# The intent is to capture the age of the person who initiated the move, either
# as a family member or a independent adult
ct_hik_1850 <- ct_1850 %>%
  filter(HIK != "") %>%
  left_join(ct_1850_hh %>% select(family, FAMILY_HEAD_AGE), by = c("family"))

ct_hik_1860 <- ct_1860 %>%
  filter(HIK != "")

intrastate_moved_1860 <- ct_hik_1850 %>%
  inner_join(ct_hik_1860, by = c("HIK"), suffix = c("_1850", "_1860")) %>%
  filter(combined_1850 != combined_1860) %>%
  mutate(migrate_age = ifelse(RELATE_1860 == 1, AGE_1850, FAMILY_HEAD_AGE))

factors <- ungroup(ct_1860_hh %>%
  group_by(town, combined) %>%
  summarise(gini = ineq(FAMILY_WEALTH, type = "Gini")) %>%
  left_join(ct_1860_hh %>%
    group_by(combined) %>%
    summarise(comb_gini = ineq(FAMILY_WEALTH, type = "Gini")), by = c("combined")) %>%
  left_join(ct_1860_hh %>%
    group_by(combined) %>%
    summarise(real_1860_gini = ineq(FAMILY_REALPROP, type = "Gini")), by = c("combined")) %>%
  left_join(ct_1850_hh %>%
    group_by(combined) %>%
    summarise(real_1850_gini = ineq(FAMILY_REALPROP, type = "Gini")), by = c("combined")) %>%
  left_join(ct_1860 %>%
    group_by(town) %>%
    summarise(
      wealth = sum(WEALTH),
      age_1860 = mean(AGE),
      pop = n()
    ), by = c("town")) %>%
  left_join(ct_1860 %>%
    group_by(combined) %>%
    summarise(
      comb_wealth = sum(WEALTH),
      comb_age_1860 = mean(AGE),
      comb_pop = n()
    ), by = c("combined")) %>%
  left_join(ct_1850_hh %>%
    group_by(combined) %>%
    summarize(num_1850_hh = n()), by = "combined") %>%
  left_join(ct_1850_hh %>%
    filter(JOB == "farm") %>%
    group_by(combined) %>%
    summarize(farm_1850_hh = n()), by = "combined") %>%
  left_join(ct_1850 %>%
    group_by(combined) %>%
    summarise(
      age_1850 = mean(AGE)
    ), by = c("combined")) %>%
  left_join(ct_1850 %>%
    filter(BPL == 414) %>%
    group_by(combined) %>%
    summarize(irish_1850 = n()), by = "combined") %>%
  left_join(ct_1850 %>%
    filter(AGE >= 20 & AGE <= 30 & SEX == 1 & BIRTH == "native" & REALPROP <= 2000) %>%
    group_by(combined) %>%
    summarize(ya_male_1850 = n()), by = "combined") %>%
  left_join(ct_1860_hh %>%
    group_by(combined) %>%
    summarize(num_1860_hh = n()), by = "combined") %>%
  left_join(ct_1860_hh %>%
    filter(JOB == "farm") %>%
    group_by(combined) %>%
    summarize(farm_1860_hh = n()), by = "combined") %>%
  left_join(ct_pop, by = "combined") %>%
  left_join(religion_1860, by = "combined") %>%
  mutate(
    wealth = wealth / pop,
    comb_wealth = comb_wealth / comb_pop,
    pct_irish_1850 = irish_1850 / POP_1850,
    pct_ya_male_1850 = ya_male_1850 / POP_1850,
    pct_farm_1850 = farm_1850_hh / num_1850_hh,
    pct_farm_1860 = farm_1860_hh / num_1860_hh
  ) %>%
  filter(combined != "UNKNOWN")) %>%
  select(
    town, combined, ends_with("gini"), ends_with("wealth"),
    ends_with("age_1850"), ends_with("age_1860"),
    starts_with("pct"), ends_with("change")
  )

save(factors, ct_1850, ct_1860, ct_1860_hh, ct_1850_hh, white_male_1860_hh, file = "ct_demographics.Rda")
