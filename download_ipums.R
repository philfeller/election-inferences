library(ipumsr)

# You must have an IPUMS account and API key to be able to work with data via the API
# https://developer.ipums.org/docs/v2/get-started/

# API key is retreived from an HCP secrets vault and defined as an environment variable
api_key <- Sys.getenv("API_KEY")
set_ipums_api_key(api_key)

# Select data only for Connecticut and limit the columns to only those that will be used
variables_1850 <- list(
  var_spec("STATEICP", case_selections = c("01")),
  "HIK", "SERIAL", "GQ", "FAMUNIT", "RELATE", "SEX", "AGE", "RACE", "BPL", "OCC1950", "REALPROP"
)

# The 1860 census added personal property value
variables_1860 <- append(variables_1850, "PERSPROP")

# Define that data will be extracted from the 1850 full-count US census
extract_1850 <- submit_extract(
  define_extract_usa(
    samples = "us1850c",
    description = "CT records from 1850",
    variables = variables_1850
  )
)

ddi_1850 <- wait_for_extract(extract_1850) %>%
  download_extract()

# Define that data will be extracted from the 1860 full-count US census
extract_1860 <- submit_extract(
  define_extract_usa(
    samples = "us1860c",
    description = "CT records from 1860",
    variables = variables_1860
  )
)

ddi_1860 <- wait_for_extract(extract_1860) %>%
  download_extract()
