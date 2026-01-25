# Start with the tidyverse image from the Rocker Project
# https://rocker-project.org/images/

FROM rocker/r-ver:4.5.2

# Install geospatial libraries needed by the 'sf' package
# and libraries needed by systemfonts and textshaping

RUN apt-get update && apt-get install -y libudunits2-dev \
    cmake \
    libabsl-dev \
    libudunits2-dev \
    libgdal-dev \
    libgeos-dev \
    libproj-dev \
    proj-data \
    proj-bin \
    libfontconfig1-dev \
    libfreetype6-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev && \
    apt-get clean && apt-get autoremove && \
    rm -rf /var/lib/apt/lists/*

# Use renv version 1.1.6
ENV RENV_VERSION=v1.1.6

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
