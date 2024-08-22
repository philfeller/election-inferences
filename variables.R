library(ipumsr)
# Variables used in other scripts are consolidated here

# Codes used in IPUMS census records
connecticut_icp <- 01
# Fairfield County
fairfield_ipcsr <- 0010
# Hartford County
hartford_ipcsr <- 0030
# Litchfield County
litchfield_ipcsr <- 0050
# Middlesex County
middlesex_ipcsr <- 0070
# New London County
new_london_ipcsr <- 0110
# Tollund County
tollund_ipcsr <- 0130
# Windham County
windham_ipcsr <- 0150

# Load lookup data used to match an IPUMS serial number to the town in
# which the household identified by that number resided. These data were
# created by comparing IPUMS data to scanned manuscript census images,
# noting that IPUMS records are in the same order.
load(file = "lookup_1850.Rda")

town_1850 <- list()
for (i in 1:nrow(lookup_1850)) {
  town <- lookup_1850[i, ]$town
  beg <- lookup_1850[i, ]$beg_ser
  end <- lookup_1850[i, ]$end_ser
  for (j in beg:end) {
    town_1850[j] <- town
  }
}

load(file = "lookup_1860.Rda")
town_1860 <- list()
for (i in 1:nrow(lookup_1860)) {
  town <- lookup_1860[i, ]$town
  beg <- lookup_1860[i, ]$beg_ser
  end <- lookup_1860[i, ]$end_ser
  for (j in beg:end) {
    town_1860[j] <- town
  }
}

# Lists, by county, of the towns that appear in the election results
litchfield_county <- c(
  "Barkhamsted", "Bethlehem", "Bridgewater", "Canaan", "Colebrook",
  "Cornwall", "Goshen", "Harwinton", "Kent", "Litchfield", "New Hartford",
  "New Milford", "Norfolk", "Plymouth", "Roxbury", "Salisbury", "Sharon",
  "Torrington", "Warren", "Washington", "Watertown", "Winchester",
  "Woodbury"
)
hartford_county <- c(
  "Avon", "Berlin", "Bloomfield", "Bristol", "Burlington", "Canton",
  "East Hartford", "East Windsor", "Enfield", "Farmington", "Glastonbury",
  "Granby", "Hartford", "Hartland", "Manchester", "Marlborough",
  "New Britain", "New Canaan", "Rocky Hill", "Simsbury", "South Windsor",
  "Southington", "Suffield", "West Hartford", "Wethersfield", "Windsor",
  "Windsor Locks"
)
tolland_county <- c(
  "Andover", "Bolton", "Columbia", "Coventry", "Ellington", "Hebron",
  "Mansfield", "Somers", "Stafford", "Tolland", "Union", "Vernon",
  "Willington"
)
windham_county <- c(
  "Ashford", "Brooklyn", "Canterbury", "Chaplin", "Eastford", "Hampton",
  "Killingly", "Plainfield", "Pomfret", "Putnam", "Sterling", "Thompson",
  "Windham", "Woodstock"
)
new_london_county <- c(
  "Bozrah", "Colchester", "East Lyme", "Franklin", "Griswold", "Groton",
  "Lebanon", "Ledyard", "Lisbon", "Lyme", "Montville", "New London",
  "North Stonington", "Norwich", "Preston", "Salem", "South Lyme",
  "Stonington", "Voluntown", "Waterford"
)
middlesex_county <- c(
  "Chatham", "Chester", "Clinton", "Cromwell", "Durham", "East Haddam",
  "Essex", "Haddam", "Killingworth", "Middletown", "Old Saybrook",
  "Portland", "Saybrook", "Westbrook"
)
new_haven_county <- c(
  "Bethany", "Branford", "Cheshire", "Derby", "East Haven", "Guilford",
  "Hamden", "Madison", "Meriden", "Middlebury", "Milford", "Naugatuck",
  "New Haven", "North Branford", "North Haven", "Orange", "Oxford",
  "Prospect", "Seymour", "Southbury", "Wallingford", "Waterbury",
  "Wolcott", "Woodbridge"
)
fairfield_county <- c(
  "Bethel", "Bridgeport", "Brookfield", "Danbury", "Darien", "Easton",
  "Fairfield", "Greenwich", "Huntington", "Monroe", "New Fairfield",
  "Newtown", "Norwalk", "Redding", "Ridgefield", "Sherman", "Stamford",
  "Stratford", "Trumbull", "Weston", "Westport", "Wilton"
)

# Lists of towns to be combine to account for new towns that were created during period
# https://portal.ct.gov/SOTS/Register-Manual/Section-VII/Connecticut-Towns-in-the-Order-of-their-Establishment

# Morris, taken from Litchfield, is in 1860 results and census
litchfield <- c("Morris", "Litchfield")
# East Granby, taken from Windsor Locks and Granby, is in 1859 results, 1860 census.
# Comparing 1855 and 1869 maps, a large majority of its population came from Granby,
# and Windsor Locks retained a large majority of its population. For purposes of
# computing demographic factors, it's safe to ignore Windsor Locks.
granby <- c("East Granby", "Granby")
# Scotland, taken from Windham, is in 1858 results, 1860 census
windham <- c("Scotland", "Windham")
# Bridgewater , taken from New Milford, is in 1857 results
new_milford <- c("Bridgewater", "New Milford")
# Putnam, taken from Killingly, Thompson, and Pomfret, in is 1856 and 1857 results
killingly_et_al <- c("Putnam", "Killingly", "Thompson", "Pomfret")
# Bethel, taken from Danbury, is in 1856 and 1857 results
danbury <- c("Bethel", "Danbury")
# South Lyme, taken from Lyme, is in 1856 and 1857 results;
# renamed Old Lyme in 1857
lyme <- c("South Lyme", "Lyme", "Old Lyme")
# West Hartford, taken from Hartford, is in 1855, 1856, and 1857 results
hartford <- c("West Hartford", "Hartford")
# Windsor Locks, taken from Windsor, is in 1855, 1856, and 1857 results
windsor <- c("Windsor Locks", "Windsor")
# Division and renaming of parts of Saybrook affect results beginning in 1853
saybrook <- c("Saybrook", "Essex", "Old Saybrook")
# Cromwell, taken from Middleton, is in 1852 results
middletown <- c("Cromwell", "Middletown")

# Locations of data files

# Full-count census data taken using ./download_ipums.R
# It can also be requested from https://usa.ipums.org/usa-action/variables/group
# and downloaded from https://usa.ipums.org/usa-action/data_requests/download

api_key <- Sys.getenv("API_KEY")
set_ipums_api_key(api_key)

# Last two extracts will be those requested by download_ipums.R
# The names of the files with the downloaded data include the extract number
ipums_extracts <- get_extract_history("usa", how_many = 2)
for (ipums_extract in ipums_extracts) {
  if (grepl("1850", ipums_extract$description)) last_1850 <- ipums_extract$number
  if (grepl("1860", ipums_extract$description)) last_1860 <- ipums_extract$number
}

ipums_1850 <- paste("./usa_", formatC(last_1850, width=5, flag="0"), ".xml", sep="")
ipums_1860 <- paste("./usa_", formatC(last_1860, width=5, flag="0"), ".xml", sep="")
# ipums_1850 <- "./usa_00023.xml"
# ipums_1860 <- "./usa_00024.xml"

# Election data taken from https://electionhistory.ct.gov/eng/contests/search/year_from:1849/year_to:1857/office_id:4/stage:et-id-3
results_file <- "./electionhistory_ct_gov_eng_contests_search_year_from_1849_year_to_1857_office_id_4_show_granularity_dt_id_1.csv"

# Parties contesting in each election year
p51 <- c("Democrat_in_1851", "Whig_in_1851", "Free_Soil_in_1851", "Abstaining_in_1851")
p52 <- c("Democrat_in_1852", "Whig_in_1852", "Free_Soil_in_1852", "Abstaining_in_1852")
p53 <- c("Democrat_in_1853", "Whig_in_1853", "Free_Soil_in_1853", "Abstaining_in_1853")
p54 <- c("Democrat_in_1854", "Whig_in_1854", "Free_Soil_in_1854", "Temperance_in_1854", "Abstaining_in_1854")
p55 <- c("Democrat_in_1855", "Whig_in_1855", "Know_Nothing_in_1855", "Abstaining_in_1855")
p56 <- c("Democrat_in_1856", "Whig_in_1856", "Know_Nothing_in_1856", "Republican_in_1856", "Abstaining_in_1856")
p57 <- c("Democrat_in_1857", "Republican_in_1857", "Abstaining_in_1857")

all_parties <- c("Democrat", "Whig", "Know Nothing", "Republican", "Free Soil", "Temperance", "Abstaining")

# Define consistent colors for each party;
# colors are from https://colorbrewer2.org/#type=qualitative&scheme=Dark2&n=7
party_colors <- c(
  "Democrat" = "#e7298a",
  "Whig" = "#1b9e77",
  "Free Soil" = "#a6761d",
  "Temperance" = "#7570b3",
  "Know Nothing" = "#d95f02",
  "Abstaining" = "#e6ab02",
  "Republican" = "#66a61e"
)

# Right-hand side of formula when evaluating residuals for covariates;
# y will be updated with appropriate dependent variable
formula_template <- as.formula("y ~ LON + gini + wealth + age_1850 + pct_farm_1850 + pct_irish_1850 + age_1860 + pct_farm_1860 + pct_irish_1860")
