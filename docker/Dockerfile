# Start with the tidyverse image from the Rocker Project
# https://rocker-project.org/images/

FROM rocker/r-ver:4.1.2

# Install geospatial libraries needed by the 'sf' package

RUN apt-get update && apt-get install -y libudunits2-dev \
    libgdal-dev \
    libgeos-dev \
    libproj-dev && \
    apt-get clean && apt-get autoremove

# Use renv version 1.0.7
ENV RENV_VERSION v1.0.7

# Create a directory named after our project directory
WORKDIR /election-inferences

# Copy the lockfile over to the Docker image
COPY renv.lock renv.lock

# Install renv and all R packages in renv.lock
RUN Rscript -e "install.packages('remotes')" \
            -e "remotes::install_github('rstudio/renv@${RENV_VERSION}')" \
            -e "renv::restore()"

# Default to bash terminal when running docker image
CMD ["bash"]
