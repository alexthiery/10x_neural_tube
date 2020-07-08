FROM rocker/tidyverse:3.6.3

LABEL authors="alex.thiery@crick.ac.uk" \
      description="Docker image containing all requirements to run 10x analysis"

ARG WHEN


# Install cellranger
RUN cd /tmp/ && \
	wget -O cellranger-3.0.2.tar.gz "http://cf.10xgenomics.com/releases/cell-exp/cellranger-3.0.2.tar.gz?Expires=1592638197&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cDovL2NmLjEweGdlbm9taWNzLmNvbS9yZWxlYXNlcy9jZWxsLWV4cC9jZWxscmFuZ2VyLTMuMC4yLnRhci5neiIsIkNvbmRpdGlvbiI6eyJEYXRlTGVzc1RoYW4iOnsiQVdTOkVwb2NoVGltZSI6MTU5MjYzODE5N319fV19&Signature=KZ5nfy4MGXjmkBoLVHdGaEy8qycVmDS6m6jBNTjsoQ3R3wZkaGIVGbtXOznjPQW~DD4IrcNtLvmkM9pZi53UBgMZq3rjtSY80Rk~H8e49Hr1PrZrmEurySz1L1T4BG~Hjv8xfsM78eFUWNQQ70Yol4pRuQz6ggv8k2ixjcBNld-mmBpQLCgdgq6wSkVjRCHUqCVRj~pLr6DAvXIdD69~4uO2hk4OcHr1AnyN~RTRQ8O6cAlK0h3hhxQ0Cig5sQu0ZJPO-LFe6Ju-S90UuAx5w3-scN~ZGXdTIknKcw9L2zi1NmCzNTnjeOUKYAojMoIscBqL3Uz6J4ku7h9vAj~oQw__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA" && \	
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

