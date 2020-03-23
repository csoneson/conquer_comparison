# Use https://hub.docker.com/r/rocker/tidyverse as base image
FROM rocker/tidyverse:3.3.2

# install R packages
COPY install.R /tmp/install.R
RUN R -f /tmp/install.R

# Set working directory for pipeline
WORKDIR /home/conquer

# Launch bash shell in WORKDIR when running container 
CMD ["bash"]
