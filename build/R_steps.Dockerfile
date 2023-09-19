# FROM rocker/r-ubuntu
# FROM bioconductor/bioconductor_docker:latest
FROM bioconductor/bioconductor_docker:3.17-R-4.3.0

# ALREADY INSTALLED
# run the following command and try loading the libraries to see
# docker run -it --rm bioconductor/bioconductor_docker:3.17-R-4.3.0 R
# stringr, tools, yaml

# NEED TO INSTALL 
# Use install2.r which is a utility tool pre-installed on the bioconductor image 
# see https://rocker-project.org/use/extending.html#install2.r 
# dplyr, readxl

RUN install2.r --error --deps TRUE dplyr
RUN install2.r --error --deps TRUE readxl
RUN install2.r --error --deps TRUE BiocManager

# Need to install using BiocManager bc its on Bioconductor
# maybe can pass --repo option to install2.r 
# library(MuData)
# library(MultiAssayExperiment)
# library(SummarizedExperiment)
# library(rhdf5)

# RUN install2.r --error --deps TRUE rhdf5
# RUN install2.r --error --deps TRUE MuData
# RUN install2.r --error --deps TRUE MultiAssayExperiment
# RUN install2.r --error --deps TRUE SummarizedExperiment

RUN R -e 'BiocManager::install(c("MultiAssayExperiment", "SummarizedExperiment", "MuData", "rhdf5"))'
RUN rm -rf /tmp/*

