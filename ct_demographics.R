# Extract demographic data from IPUMS full-count census records for Connecticut.
# Variables and function defined by this script are used in other scripts:
#   get_combined - function that assigns a name to a town that changes shape
#   ct_eligible - tibble with estimated number of eligible voters by town and year
#   factors - tibble with calculated informtion about a town, such as gini

source("./global.R")
source("./census_utils.R", local = TRUE)

# Load IPUMS data files; suppress messages about the API key

ddi1850 <- ipumsr::read_ipums_ddi(ipums_1850)
ddi1860 <- ipumsr::read_ipums_ddi(ipums_1860)

# Create a tibble with individual 1850 census data
ct_1850 <- ipumsr::read_ipums_micro(ddi1850, verbose = FALSE) %>%
  select(
    HIK, SERIAL, GQ, FAMUNIT, RELATE, SEX, AGE, RACE, BPL,
    OCC1950, REALPROP, SCHOOL, LIT, PAUPER, CRIME
  ) %>%
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
  dplyr::filter(town != "NULL") %>%
  arrange(SERIAL, FAMUNIT, RELATE)

# Create an intermediate tibble with 1850 individual data summarized by combined,
# a town name that is consistent for 1850 and 1860, taking boundary changes into
# account.
ct_1850_combined <- ct_1850 %>%
  group_by(combined) %>%
  summarise(
    POP_1850 = n(),
    age_1850 = mean(AGE),
    irish_1850 = sum(BPL == 414),
    ya_male_1850 = sum(AGE >= 20 & AGE <= 30 & SEX == 1 & BIRTH == "native" & REALPROP <= 2000)
  )

# Load file with data that is missing from IPUMS: first pages for Brooklyn and Hebron
load(file = "missing_1860_ipums_rows.Rda")

# The IPUMS data contains some transcription errors. Although GINI indices calculated
# from IPUMS data vary insignificantly from indices calculated with correct real and
# personal property values, update transcription errors for Meriden with corrected
# values.
ct_1860 <- ipumsr::read_ipums_micro(ddi1860, verbose = FALSE) %>%
  mutate(RELATE = ifelse(SERIAL == 294445, 1, RELATE)) %>%
  select(
    HIK, SERIAL, GQ, FAMUNIT, RELATE, SEX, AGE, RACE, BPL,
    OCC1950, REALPROP, SCHOOL, LIT, PAUPER, CRIME, PERSPROP
  ) %>%
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
  dplyr::filter(town != "NULL") %>%
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

# Construct tibble with 1850 household wealth and demographics data
ct_1850_hh <- ct_1850 %>%
  # Exclude servants and institutional housing
  dplyr::filter((GQ %in% c(1, 2, 5) & FAMUNIT == 1) | GQ == 4) %>%
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

# Create an intermediate tibble with 1850 household data summarized by combined,
# a town name that is consistent for 1850 and 1860, taking boundary changes into
# account.
ct_1850_hh_combined <- ct_1850_hh %>%
  group_by(combined) %>%
  summarise(
    num_1850_hh = n(),
    farm_1850_hh = sum(JOB == "farm"),
    real_1850_gini = ineq::ineq(FAMILY_REALPROP, type = "Gini")
  )

# Create an intermediate tibble with summarized individual and household data
ct_1850_summary <- left_join(ct_1850_combined, ct_1850_hh_combined, by = "combined")

# Construct tibble with 1860 household wealth and demographics data
ct_1860_hh <- ct_1860 %>%
  # Exclude servants and institutional housing
  dplyr::filter((GQ %in% c(1, 2, 5) & FAMUNIT == 1) | GQ == 4) %>%
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

# Create an intermediate tibble with 1860 household data summarized by combined,
# a town name that is consistent for 1850 and 1860, taking boundary changes into
# account.
ct_1860_hh_combined <- ct_1860_hh %>%
  group_by(combined) %>%
  summarise(
    num_1860_hh = n(),
    farm_1860_hh = sum(JOB == "farm"),
    real_1860_gini = ineq::ineq(FAMILY_REALPROP, type = "Gini"),
    comb_gini = ineq::ineq(FAMILY_WEALTH, type = "Gini"),
    POP_1860 = n()
  )

# Count number of native-born, white, adult males to estimate eligible voters;
# the number of naturalized foreign-born citizens who meet the residency
# requirement will be estimated in the estimate_voters function.
ct_eligible <- ct_1850 %>%
  dplyr::filter(
    SEX == 1,
    BPL < 100,
    RACE == 1,
    AGE > 20
  ) %>%
  group_by(combined) %>%
  summarize(ELIG_1850 = n()) %>%
  inner_join(ct_1850 %>%
    dplyr::filter(
      SEX == 1,
      BPL >= 100,
      RACE == 1,
      AGE > 20
    ) %>%
    group_by(combined) %>%
    summarize(FOREIGN_1850 = n()), by = "combined") %>%
  inner_join(ct_1860 %>%
    dplyr::filter(
      SEX == 1,
      BPL < 100,
      RACE == 1,
      AGE > 20
    ) %>%
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

# IPUMS has a field that shows individuals who have been algorithmically matched
# across censuses. Because matching is partial, this is not currently being
# used (compute_migration is set to FALSE in global.R), and these tibbles aren't
# created in order to speed execution.
if (compute_migration == TRUE) {
  # Find information about people who moved out of state between 1850 and 1860.
  # The 1860 records need to have an HIK variable in both censuses.
  # The intent is to capture the age of the person who initiated the move, either
  # as a family member or a independent adult
  ct_hik_1850 <- ct_1850 %>%
    dplyr::filter(HIK != "") %>%
    left_join(ct_1850_hh %>% select(family, FAMILY_HEAD_AGE), by = c("family"))

  ct_hik_1860 <- ct_1860 %>%
    dplyr::filter(HIK != "")

  intrastate_moved_1860 <- ct_hik_1850 %>%
    inner_join(ct_hik_1860, by = c("HIK"), suffix = c("_1850", "_1860")) %>%
    dplyr::filter(combined_1850 != combined_1860) %>%
    mutate(migrate_age = ifelse(RELATE_1860 == 1, AGE_1850, FAMILY_HEAD_AGE))

  net_intrastate_migration <- intrastate_moved_1860 %>%
    group_by(combined_1850) %>%
    summarise(number = n()) %>%
    select(combined_1850, number) %>%
    full_join(intrastate_moved_1860 %>%
      group_by(combined_1860) %>%
      summarise(number = n()) %>%
      select(combined_1860, number), join_by(combined_1850 == combined_1860)) %>%
    mutate(net_migration = number.y - number.x) %>%
    rename(combined = combined_1850) %>%
    select(combined, net_migration)
}

# Create a tibble with economic and demographic data by town, to be used as
# inference covariates and in other analysis.
factors <- ungroup(ct_1860_hh %>%
  group_by(town, combined) %>%
  summarise(gini = ineq::ineq(FAMILY_WEALTH, type = "Gini")) %>%
  left_join(ct_1860_hh_combined, by = "combined") %>%
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
  left_join(ct_1850_summary, by = "combined") %>%
  left_join(religion_1860, by = "combined") %>%
  mutate(
    wealth = wealth / pop,
    comb_wealth = comb_wealth / comb_pop,
    pct_irish_1850 = irish_1850 / POP_1850,
    pct_ya_male_1850 = ya_male_1850 / POP_1850,
    pct_farm_1850 = farm_1850_hh / num_1850_hh,
    pct_farm_1860 = farm_1860_hh / num_1860_hh
  ) %>%
  dplyr::filter(combined != "UNKNOWN")) %>%
  select(
    town, combined, ends_with("gini"), ends_with("wealth"),
    ends_with("age_1850"), ends_with("age_1860"),
    starts_with("pct"), ends_with("change")
  )

save(factors, ct_1850, ct_1860, ct_1860_hh, ct_1850_hh, file = "ct_demographics.Rda")
