# Variables used in other scripts are consolidated here

library(ipumsr)

# Suppress warnings; the ei.MD.bayes function produces uninformative warnings because the
# row and column marginals are proportions.
# Suppress messages about the default group that deplyr.summarise chooses.
options(
  warn = -1,
  dplyr.summarise.inform = FALSE
)

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
  "New Britain", "Rocky Hill", "Simsbury", "South Windsor",
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
  "Fairfield", "Greenwich", "Huntington", "Monroe", "New Canaan", "New Fairfield",
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
# North Canaan, taken from Canaan, is in 1860 census
canaan <- c("Canaan", "North Canaan")
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

# 1860 IPUMS data has bad birthplace transcriptions for Avon, Burlington,
# Farmington, and New Britain, and 1860 census data can't be used to estimate
# eligible voters for those towns. Taxable polls have been used as a proxy for
# eligible voters, and the rate of change in polls for these four towns will be
# used to estimate the annual rate of change of change for eligible voters.
# Numbers for 1852 and 1861 are taken from town statistics in the 1854 and 1862
# Connecticut Registers:
# https://babel.hathitrust.org/cgi/pt?id=hvd.li2mn1&view=1up&seq=105
# https://babel.hathitrust.org/cgi/pt?id=nyp.33433081898748&seq=127

avon_poll_change <- (201 - 197) / 9
burlington_poll_change <- (284 - 208) / 9
farmington_poll_change <- (645 - 534) / 9
new_britain_poll_change <- (1277 - 628) / 9
bad_birthplace <- c("Avon", "Burlington", "Farmington", "New Britain")

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
  if (grepl("1850 CT", ipums_extract$description)) last_1850 <- ipums_extract$number
  if (grepl("1860 CT", ipums_extract$description)) last_1860_ct <- ipums_extract$number
}

ipums_data_path <- "./data/licensed-data/ipums/"
if (exists("last_1850")) {
  ipums_1850 <- paste(ipums_data_path, "usa_", formatC(last_1850, width = 5, flag = "0"), ".xml", sep = "")
}
if (exists("last_1860_ct")) {
  ipums_1860 <- paste(ipums_data_path, "usa_", formatC(last_1860_ct, width = 5, flag = "0"), ".xml", sep = "")
}

# Election data taken from https://electionhistory.ct.gov/eng/contests/search/year_from:1849/year_to:1857
results_file <- "./data/electionhistory_ct_gov_eng_contests_search_year_from_1849_year_to_1857_stage_et_id_3_show_granularity_dt_id_1.csv.gz"

# Religious accomodation data was taken from images of the 1860 Census Social Statistics schedule for each town.
religion_file <- "./data/1860_CT_religious_accomodation.csv"

# Parties contesting in each election year
p51 <- c("Democrat_in_1851", "Whig_in_1851", "Free_Soil_in_1851", "Abstaining_in_1851")
p52 <- c("Democrat_in_1852", "Whig_in_1852", "Free_Soil_in_1852", "Abstaining_in_1852")
p53 <- c("Democrat_in_1853", "Whig_in_1853", "Free_Soil_in_1853", "Abstaining_in_1853")
p54 <- c("Democrat_in_1854", "Whig_in_1854", "Free_Soil_in_1854", "Temperance_in_1854", "Abstaining_in_1854")
p55 <- c("Democrat_in_1855", "Whig_in_1855", "Know_Nothing_in_1855", "Abstaining_in_1855")
p56 <- c("Democrat_in_1856", "Whig_in_1856", "Know_Nothing_in_1856", "Republican_in_1856", "Abstaining_in_1856")
p57 <- c("Democrat_in_1857", "Republican_in_1857", "Abstaining_in_1857")

all_parties <- c("Democrat", "Whig", "Know Nothing", "Republican", "Free Soil", "Temperance", "Abstaining")

# Map of candidate name variants to consistent names
alias_map <- list(
  "Lafayette S. Foster" = c("L. S. Foster", "Lafayette Foster"),
  "Samuel Ingham" = c("Samuel D. Ingham", "Samuel Engraham"),
  "John M. Niles" = c("J. M. Niles"),
  "Thomas H. Seymour" = c("Jonas H. Senor", "T. H. Seymour", "T. Seymour", "Thomas Seymour", "Thomas S. Seymour", "Ths. H. Seymour", "Tom H. Seymour"),
  "Alexander W Holley" = c("A H Holly", "Alexander H. Holley"),
  "Charles H. Pond" = c("C H Pond"),
  "John T. Wait" = c("J. T. Wait", "John T Waite"),
  "Josiah M. Carter" = c("Joseph M Carter", "Joshua M Carter", "Josiah B Carter"),
  "William Alexander" = c("William Alexander Jr."),
  "Hiram Weed" = c("Miriam Weed"),
  "J. Hammond Trumbull" = c("Hammond Trumbull", "Jonathan H. Trumbull"),
  "John Boyd" = c("John D. Boyd"),
  "Joseph W. White" = c("Joseph White"),
  "Nehemiah D. Sperry" = c("N. D. Sperry", "Nehemiah Sperry", "Jeremiah D. Sperry"),
  "Orville H. Platt" = c("O. H. Platt", "Orville Platt"),
  "Amos Townsend, Jr." = c("Amos Townsend"),
  "Arthur B. Calef" = c("A.B. Calef", "Arthur B. Caleff", "Arthur B. Califf"),
  "Daniel W. Camp" = c("Daniel Camp"),
  "Frederick P. Coe" = c("Frederick A. Poe"),
  "George Reed" = c("George Read"),
  "Jesse G. Baldwin" = c("Jesse Baldwin", "Jesse D. Baldwin", "Jesse W. Baldwin", "J. G. Baldwin"),
  "Joseph B. Gay" = c("Joseph Gay"),
  "Talcott Crosby" = c("Talcott L. Crosby"),
  "Thomas Clark" = c("Th. Clark")
)

# Data frame assigning parties to candidates for each office and year
nominated_party_assignments <- data.frame(
  candidate_name = c(
    # Governor
    rep("Thomas H. Seymour", 5),
    rep("Samuel Ingham", 4),
    rep("Lafayette S. Foster", 3),
    "Green Kendrick",
    rep("Henry Dutton", 3),
    "John A. Rockwell",
    "John M. Niles",
    rep("John Boyd", 2),
    rep("Francis Gillette", 2),
    "John Hooker",
    "Charles Chapman",
    rep("William T. Minor", 2),
    "Gideon Welles",
    "Alexander H. Holley",
    # Lt. Governor
    rep("Charles H. Pond", 5),
    rep("John T. Wait", 4),
    "Thomas Backus",
    rep("Green Kendrick", 2),
    "Roger H. Mills",
    "Nathaniel I Kellogg",
    rep("Alexander W Holley", 2),
    "Josiah M. Carter",
    rep("John Boyd", 2),
    "Julius Clark",
    rep("William Field", 2),
    "William Alexander",
    "Erastus Lester",
    "William Field",
    "Albert Day",
    "Henry B. Harrison",
    "Alfred A. Burnham",
    # Secretary of State
    rep("Hiram Weed", 2),
    rep("John P. C. Mather", 3),
    rep("Roger Averill", 4),
    rep("Roger H. Mills", 3),
    rep("J. Hammond Trumbull", 2),
    "Thomas Robinson",
    rep("Oliver H. Perry", 2),
    "Edmund Perkins",
    "Joseph W. White",
    "George Read",
    "Walter Webb",
    "John S. Baldwin",
    rep("David Lyman", 2),
    "Edward S. Moseley",
    rep("Nehemiah D. Sperry", 2),
    "John Boyd",
    "Orville H. Platt",
    # Treasurer
    rep("Henry D. Smith", 3),
    rep("Edwin Stearns", 2),
    rep("Talcott Crosby", 4),
    "Stephen Taylor",
    rep("Thomas Clark", 3),
    rep("Daniel W. Camp", 3),
    "Stephen D. Pardee",
    rep("George Reed", 2),
    rep("Jesse G. Baldwin", 2),
    "Walter Webb",
    rep("Amos Townsend, Jr.", 2),
    "Thomas Clark",
    "Arthur B. Calef",
    "Frederick P. Coe",
    "Joseph B. Gay",
    "Frederick S. Wildman"
  ),
  office_name = c(
    rep("Governor", 28),
    rep("Lieutenant Governor", 28),
    rep("Secretary of the State", 29),
    rep("Treasurer", 29)
    # rep("Comptroller", 8)
  ),
  yr = c(
    # Governor
    1849:1853,
    1854:1857,
    1849:1851,
    1852,
    1853:1855,
    1856,
    1849,
    1850:1851,
    1852:1853,
    1854,
    1854,
    1855:1856,
    1856,
    1857,
    # Lt. Governor
    1849:1853,
    1854:1857,
    1849,
    1850:1851,
    1852,
    1853,
    1854:1855,
    1856,
    1849, 1853,
    1850,
    1851:1852,
    1854,
    1854,
    1855,
    1856,
    1856,
    1857,
    # Secretary of State
    1849:1850,
    1851:1853,
    1854:1857,
    1849:1851,
    1852, 1856,
    1853,
    1854:1855,
    1849,
    1850,
    1851,
    1852,
    1853,
    1854:1855,
    1854,
    1855:1856,
    1856,
    1857,
    # Treasurer
    1849:1851,
    1852:1853,
    1854:1857,
    1849,
    1850:1852,
    1853:1855,
    1856,
    1849:1850,
    1851:1852,
    1853,
    1854:1855,
    1854,
    1855:1856,
    1856:1857
  ),
  candidate_party = c(
    # Governor
    rep("Democrat_votes", 9),
    rep("Whig_votes", 8),
    rep("Free_Soil_votes", 6),
    "Temperance_votes",
    rep("Know_Nothing_votes", 2),
    rep("Republican_votes", 2),
    # Lt. Governor
    rep("Democrat_votes", 9),
    rep("Whig_votes", 8),
    rep("Free_Soil_votes", 6),
    "Temperance_votes",
    rep("Know_Nothing_votes", 2),
    rep("Republican_votes", 2),
    # Secretary of State
    rep("Democrat_votes", 9),
    rep("Whig_votes", 8),
    rep("Free_Soil_votes", 7),
    "Temperance_votes",
    rep("Know_Nothing_votes", 2),
    rep("Republican_votes", 2),
    # Treasurer
    rep("Democrat_votes", 9),
    rep("Whig_votes", 8),
    rep("Free_Soil_votes", 7),
    "Temperance_votes",
    rep("Know_Nothing_votes", 2),
    rep("Republican_votes", 2)
  ),
  stringsAsFactors = FALSE
)
nominated_party_assignments$yr <- as.integer(nominated_party_assignments$yr)

# Data frame assigning parties to candidates in towns that voted in large numbers
# for someone other than the party nominee; party assignments were determined
# by looking at other races in which the candidate ran, by known affiliations, or
# by looking at vote patterns in towns where the candidate did well.
imputed_party_assignments <- data.frame(
  candidate_name = c(
    # Secretary of the State
    "Thomas Black",
    "Erwin Swames",
    "C. C. Cornstock",
    "Daniel W. Camp",
    "Walter Webb",
    rep("Talcott Crosby", 2),
    "Thomas Clark",
    "Daniel H. Camp",
    "Amos Townsend, Jr.",
    "Arthur B. Calef",
    "Joseph B. Gay",
    "Stephen D. Pardre",
    "Narnmore Trumbull",
    # Treasurer
    rep("J. Hammond Trumbull", 2),
    "John J.C. Mather",
    "Walter Booth",
    "Charles Robinson",
    "Jesse G. Baldwin",
    "Roger Averill",
    "Edward S. Mosely",
    "Nehemiah D. Sperry",
    "John Boyd"
  ),
  office_name = c(
    rep("Secretary of the State", 14),
    rep("Treasurer", 10)
  ),
  yr = c(
    # Secretary of the State
    1852,
    rep(1853, 4),
    1854, 1855,
    rep(1854, 3),
    1855,
    rep(1856, 3),
    # Treasurer
    1852, 1856,
    rep(1853, 4),
    rep(1854, 2),
    1855,
    1856
  ),
  candidate_party = c(
    # Secretary of the State
    "Whig_votes",
    "Democrat_votes",
    "Free_Soil_votes",
    "Whig_votes",
    "Free_Soil_votes",
    rep("Democrat_votes", 2),
    "Temperance_votes",
    "Whig_votes",
    "Free_Soil_votes",
    "Know_Nothing_votes",
    "Republican_votes",
    rep("Whig_votes", 2),
    # Treasurer
    rep("Whig_votes", 2),
    "Democrat_votes",
    "Free_Soil_votes",
    "Whig_votes",
    "Free_Soil_votes",
    "Democrat_votes",
    "Temperance_votes",
    "Know_Nothing_votes",
    "Republican_votes"
  ),
  stringsAsFactors = FALSE
)
imputed_party_assignments$yr <- as.integer(imputed_party_assignments$yr)

party_assignments <- rbind(
  nominated_party_assignments,
  imputed_party_assignments
)

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
