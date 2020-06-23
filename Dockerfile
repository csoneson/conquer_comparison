# Use https://hub.docker.com/r/rocker/tidyverse as base image
FROM rocker/tidyverse:3.3.2

# install dependencies
RUN apt-get update && apt-get -y --no-install-recommends install \
    xauth \
    libxt-dev \
    libcairo2-dev \
    libnlopt-dev \
    # Clean up
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*
    
# install R packages
COPY install.R /tmp/install.R
RUN R -f /tmp/install.R

# install Seurat 1.4.0
RUN cd /tmp \
    && wget https://github.com/satijalab/seurat/archive/v1.4.0.tar.gz \
    && R CMD INSTALL v1.4.0.tar.gz \
    && rm -rf /tmp/v1.4.0.tar.gz

# Set working directory for pipeline
WORKDIR /home/conquer

# Launch bash shell in WORKDIR when running container 
CMD ["bash"]
