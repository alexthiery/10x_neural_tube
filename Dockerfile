FROM rocker/tidyverse:3.6.3

LABEL authors="alex.thiery@crick.ac.uk" \
      description="Docker image containing all requirements to run 10x analysis"

ARG WHEN


# Install cellranger
RUN cd /tmp/ && \
	wget -O cellranger-3.0.2.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-3.0.2.tar.gz?Expires=1594958998&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cHM6Ly9jZi4xMHhnZW5vbWljcy5jb20vcmVsZWFzZXMvY2VsbC1leHAvY2VsbHJhbmdlci0zLjAuMi50YXIuZ3oiLCJDb25kaXRpb24iOnsiRGF0ZUxlc3NUaGFuIjp7IkFXUzpFcG9jaFRpbWUiOjE1OTQ5NTg5OTh9fX1dfQ__&Signature=g5grwCXWbt-grrRm6Szx1VUqiwffXGUvrJo4FcFuA-~ydhVw0Zl~bbrWX6~ds-0d5AlvDWxjm0AL2DvmhdzmlaQDT2tzkFkwMsrwkuKyOfdZ~9WBGTZl6FXSQD-TJgHmGP9SuzrMvUL38Ryh3NOmAGuBh3a-fhK66xQ7Kpyk~d5HyIuGdi8RRKfVL48WP-EMrsJzhkDonSd352slctD0jiMxJqR-eBrqsjanX1Xs46sPzm5ZxUWqJPoqeYV4WBOlgm1RNqlPtF9RddOVPaUjmALoZR3VrFR4s2ZLSbC~JELxxbueC3oK6uHAtMFMdU4do1ha01vvu5yEY5x~CWxBFA__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA" && \	
	mv cellranger-3.0.2.tar.gz /opt/ && \
	cd /opt/ && \
	tar -xzvf cellranger-3.0.2.tar.gz && \
	rm -f cellranger-3.0.2.tar.gz

# Set path
ENV PATH /opt/cellranger-3.0.2:$PATH


# Install apt packages
RUN apt-get update \
 && apt-get install -y --no-install-recommends \
 git=1:2.20.1-2+deb10u3 \
 apt-utils=1.8.2 \
 unzip=6.0-23+deb10u1 \
 procps=2:3.3.15-2 \
 build-essential=12.6 \
 zlib1g-dev=1:1.2.11.dfsg-1 \
 python3-pandas \
 libxt-dev \
 libmagick++-dev

RUN R -e "devtools::install_version('Seurat', version = '3.1.5', dependencies= NA)"
RUN R -e "devtools::install_version('future', version = '1.17.0', dependencies= T)"
RUN R -e "devtools::install_version('cowplot', version = '1.0.0', dependencies= T)"
RUN R -e "devtools::install_version('clustree', version = '0.4.2', dependencies= NA)"
RUN R -e "devtools::install_version('gridExtra', version = '2.3', dependencies=T)"
RUN R -e "devtools::install_version('getopt', version = '1.20.3', dependencies=T)"
RUN R -e "devtools::install_version('pheatmap', version = '1.0.12', dependencies=T)"
RUN R -e "BiocManager::install('limma')"