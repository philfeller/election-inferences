# 1850s CT election inferences

## Getting started

This project has been tested using the Docker image that is defined in Dockerfile. That same image could also be used to create a container in which to execute the code for the elections project.
- Download and install [Docker Desktop](https://www.docker.com/products/docker-desktop/)
- docker pull ghcr.io/philfeller/election-inferences

If you do not use the Docker container, download and install [RStudio](https://posit.co/products/open-source/rstudio/), You will also need to install the R packages that are installed in the Docker image. The renv.lock file has all of the necessary packages and their dependencies, and it can be used to install them in another R instance. See the [renv package documentation](https://rstudio.github.io/renv/articles/renv.html) for more information on how to use this file.

## Installation

Download the project's source code, either by cloning this repo or by downloading a file that contains the source code. Note that the *.Rda files are stored using LFS; after executing the 'git clone' command, execute 'git lfs install && git lfs fetch --all && git lfs checkout'.

If you choose to use the Docker image:
- Create a container from the image.
- Use 'docker cp' to copy the source code to the container.

If you choose to use RStudio, set the directory with the source code as the working directory.

## Using the programs

A typical execution is done in three steps:
1. Create and download an extract from IPUMS US census data.
```
export API_KEY=your API key
Rscript download_ipums.R
```
2. Calculate inferences:
```
Rscript voter_regressions.R
```
3. Analyze the data generated by the previous step:
```
Rscript analyze_results.R
```
## Files included in project

### R
- lookup_1850.Rda - Defines the range of IPUMS serial numbers for a particular town in 1850
- lookup_1860.Rda - Defines the range of IPUMS serial numbers for a particular town in 1860
- missing_1860_ipums_rows.Rda - Data that is missing from IPUMS digitized transcription
- download_ipums.R - Download census data from IPUMS; requires an API key and that the key is defined in environment variable API_KEY.
- variables.R - Defines variables that are used elsewhere, including the paths for the files that were loadloaded from IPUMS.
- ct_demographics.R - Creates Tidyverse tibbles from the IPUMS data, summarizing by town and calculating information, such as GINI index.
- prepare_results.R - Defines functions that take election results for a pair of year and prepare them for analysis.
- betas.R - Defines functions to extract beta values from the MCMC data that is created by the eiPack inference.
- present.R - Defines functions for presenting betas in various formats.
- prepare_results.R - This is the script that calculate election inferences. It sources the scripts that need to have been run first.
- inferences.Rda - Output from a typical execution of prepare_results.R and the data that was used in the statistical analysis for the paper.
- analyze_results.R - Loads the inferences.Rda data file that was created by prepare_results.R and performs the statistical analysis that is used in the paper.
- test.R - Contains one simple inference and a few graphs; it is used by the CI/CD pipeline to test that all needed files are included and working as expected.

### IPUMS data
Census population data come from [IPUMS](https://usa.ipums.org/usa/). These files can be used if you do not want to set up an IPUMS API keys and download your own data. If you use them, you will need to modify variables.R.
- usa_00023.xml
- usa_00023.dat.gz
- usa_00024.xml
- usa_00024.dat.gz

### Other census data
The Social Statistics schedule includes the number of religious accommodations (seats) by demonination; and schedule images are available on [FamilySearch](https://www.familysearch.org/records/images/search-results?page=1&place=346&endDate=1860&startDate=1860&creator=Federal%20Census). These data were hand-entered into a CSV file.
- 1860_CT_religious_accomodation.csv

### Elections data
These data were downloaded from [the Connecticut Secetary of State](https://electionhistory.ct.gov/eng/contests/search/year_from:1849/year_to:1857/office_id:4/stage:et-id-3).
- electionhistory_ct_gov_eng_contests_search_year_from_1849_year_to_1857_office_id_4_show_granularity_dt_id_1.csv

### ESRI shapefiles
These files were created from a [Connecticut GIS file](https://ct-deep-gis-open-data-website-ctdeep.hub.arcgis.com/maps/82672ae5f3764021b9a4804f524f928b/about). Where towns were created from parts of other towns, I've drawn boundaries that conform as closely as possible to historic maps or town histories.
- 1851_CT_towns.shp
- 1851_CT_towns.shx
- 1851_CT_towns.dbf
- 1852_CT_towns.shp
- 1852_CT_towns.shx
- 1852_CT_towns.dbf
- 1854_CT_towns.shp
- 1854_CT_towns.shx
- 1854_CT_towns.dbf
- 1855_CT_towns.shp
- 1855_CT_towns.shx
- 1855_CT_towns.dbf
- 1856_CT_towns.shp
- 1856_CT_towns.shx
- 1856_CT_towns.dbf
- 1857_CT_towns.shp
- 1857_CT_towns.shx
- 1857_CT_towns.dbf

### Docker files
- docker/Dockerfile - Docker file used to build the ghcr.io/philfeller/election-inferences image
- docker/renv.lock - File that defines the R packages to be installed in the image
- runner/Dockerfile - Docker file used to build a self-hosted runner image
- runner/entrypoint.sh - BASH script used within the image to register and start a self-hosted runner

### GitHub files
- .github/workflows/deploy-image.yml - Builds and deploys the main Docker image
- .github/workflows/deploy-runner-image.yml - Builds and deploys the Docker image used for runners
- .github/workflows/test-code.yml - Tests the integrity of the code on a self-hosted runner using the ghcr.io/philfeller/election-inferences image.
- .gitattributes - Defines the files that are stored using LFS
