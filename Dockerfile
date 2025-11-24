# Use rocker/shiny as base - includes R and Shiny Server
FROM rocker/shiny:4.3.2

# Install system dependencies
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    && rm -rf /var/lib/apt/lists/*

# Install required R packages
RUN R -e "install.packages(c('shiny', 'deSolve', 'ggplot2', 'tidyr', 'dplyr', 'scales', 'shinyBS', 'bslib', 'httpuv'), repos='https://cloud.r-project.org/')"

# Copy app files to Shiny server directory
COPY app.R /srv/shiny-server/

# Make port 3838 available
EXPOSE 3838

# Run Shiny Server
CMD ["/usr/bin/shiny-server"]