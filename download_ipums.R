library(ipumsr)
library(magrittr)

# You must have an IPUMS account and API key to be able to work with data via the API
# https://developer.ipums.org/docs/v2/get-started/

# API key is retreived from GitHub secrets and defined as an environment variable
# The API key is set in variables.R
source("./variables.R")

# Specify the variables to be downloaded in the data extracts
ipums_vars <- list("HIK", "SERIAL", "GQ", "FAMUNIT", "RELATE", "SEX", "AGE", "RACE", "BPL", "OCC1950", "REALPROP",
                "SCHOOL", "LIT", "PAUPER", "CRIME")

# Select data only for Connecticut and limit the columns to only those that will be used
variables_1850 <- list(var_spec("STATEICP", case_selections = c("01")),
                       "HIK", "SERIAL", "GQ", "FAMUNIT", "RELATE", "SEX", "AGE", "RACE", "BPL",
                       "OCC1950", "REALPROP", "SCHOOL", "LIT", "PAUPER", "CRIME")

# The 1860 census added personal property value
variables_1860 <- append(variables_1850, "PERSPROP")

# Define that data will be extracted from the 1850 full-count US census
extract_1850 <- submit_extract(
  define_extract_usa(
    samples = "us1850c",
    description = "1850 CT records",
    variables = variables_1850
  )
)

ddi_1850 <- wait_for_extract(extract_1850) %>%
  download_extract(download_dir = ipums_data_path)

# Define that data will be extracted from the 1860 full-count US census
extract_1860 <- submit_extract(
  define_extract_usa(
    samples = "us1860c",
    description = "1860 CT records",
    variables = variables_1860
  )
)

ddi_1860 <- wait_for_extract(extract_1860) %>%
  download_extract(download_dir = ipums_data_path)

