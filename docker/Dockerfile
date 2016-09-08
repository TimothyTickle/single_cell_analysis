FROM ubuntu:xenial
MAINTAINER ttickle@broadinstitute.org

######################
## Environment
######################

## Constants
### SOFTWARE versions
ENV R_VERSION 3.3.1-1xenial0

### locations
ENV BIN /usr/local/bin
ENV DATA /usr/local/data
ENV SRC /usr/local/src

######################
## Dependencies and Tools
######################
##############
## Helper tools
RUN apt-get update && \
    apt-get install -y unzip && \
    apt-get install -y wget && \
    apt-get install -y git

##############
## System tools
## devtools: libssl-dev, libcurl4-openssl-dev, libxml2-dev
## Cairo (scde): libxml2-dev 
## rstudio: libjpeg62, libgstreamer0.10-0, libgstreamer-plugins-base0.10-0
RUN apt-get install -y libssl-dev && \
    apt-get install -y libcurl4-openssl-dev && \
    apt-get install -y libxml2-dev && \
    apt-get install -y libxt-dev && \
    apt-get install -y libgstreamer0.10-0 && \
    apt-get install -y libjpeg62 && \
    apt-get install -y ibgstreamer-plugins-base0.10-0

## For SCDE
RUN apt-get install -y libcairo2-dev

##############
## R for scone needs R >=3.3
RUN echo "deb http://cran.rstudio.com/bin/linux/ubuntu xenial/" | tee -a /etc/apt/sources.list && \
    gpg --keyserver keyserver.ubuntu.com --recv-key E084DAB9 && \
    gpg -a --export E084DAB9 | apt-key add - && \
    apt-get update && \
    apt-get install -y r-base=${R_VERSION} && \
    apt-get install -y r-base-dev=${R_VERSION}

##############
## devtools
## scone
## seurat
## RaceID
## FastProject
RUN echo "source(\"https://bioconductor.org/biocLite.R\")" > ${SRC}/install_pkgs.R && \
    echo "biocLite(\"devtools\")" >> ${SRC}/install_pkgs.R && \
    echo "biocLite()" >> ${SRC}/install_pkgs.R && \
    echo "biocLite(\"YosefLab/scone\", dependencies=TRUE)" >> ${SRC}/install_pkgs.R && \
    echo "devtools::install_github(\"satijalab/seurat\")" >> ${SRC}/install_pkgs.R && \
    echo "biocLite(c(\"RcppArmadillo\", \"Rook\", \"rjson\", \"RMTstat\", \"extRemes\", \"pcaMethods\", \"Cairo\"), dependencies=TRUE)" >> ${SRC}/install_pkgs.R && \
    echo "biocLite(c(\"pheatmap\", \"amap\"), dependencies=TRUE)" >> ${SRC}/install_pkgs.R && \
    echo "biocLite()" >> ${SRC}/install_pkgs.R && \
    Rscript ${SRC}/install_pkgs.R

WORKDIR ${SRC}
RUN wget https://github.com/hms-dbmi/scde/archive/1.99.1.tar.gz && \
    R CMD INSTALL 1.99.1.tar.gz
    
## FastProject
# Needs python
WORKDIR ${BIN}
RUN wget https://rawgit.com/YosefLab/FastProject/master/FP_Linux_Install.sh && \
    yes | bash FP_Linux_Install.sh

## RaceID
WORKDIR ${BIN}
RUN git clone https://github.com/dgrun/RaceID

## RStudio
RUN wget https://download1.rstudio.org/rstudio-0.99.903-amd64.deb && \
    dpkg -i rstudio-0.99.903-amd64.deb
