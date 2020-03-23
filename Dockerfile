# Use https://hub.docker.com/r/rocker/r-ver as base image
FROM rocker/r-ver:3.3.2

# install R packages
COPY install.R /tmp/install.R
RUN R -f /tmp/install.R
