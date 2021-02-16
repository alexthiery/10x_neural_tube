FROM rocker/tidyverse:3.6.3

LABEL authors="alex.thiery@crick.ac.uk" \
      description="Docker image containing all requirements to run 10x analysis"


# Install cellranger
RUN cd /tmp/ && \
	wget -O cellranger-3.0.2.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-3.0.2.tar.gz?Expires=1613526989&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cHM6Ly9jZi4xMHhnZW5vbWljcy5jb20vcmVsZWFzZXMvY2VsbC1leHAvY2VsbHJhbmdlci0zLjAuMi50YXIuZ3oiLCJDb25kaXRpb24iOnsiRGF0ZUxlc3NUaGFuIjp7IkFXUzpFcG9jaFRpbWUiOjE2MTM1MjY5ODl9fX1dfQ__&Signature=b0iFXuBx48UbKDawH2e8AuhIWhSJHaeF~zSFjupoQ--K92PD5XQ6w6gJQHavTlnRLAOA4UCUdVwiCZLA4PWIroFFmm7DeBDlrt9i-1XdCci9VCa2S0NEMypjUZNfyyI4SHIRJFd5CETP2cqUFqrp1U1gWB2czzmuJXrY0v4ef~wSQq8sFdJr3YHKlJrf4J2U2JP9hbdT-Hgwi9p9bJkbcY5wLgYe~FXg-aWXN5h91lcIc5c4Q0oBVJPZQlEVercJDnbmlRO-yz24g7koOsMaozfLYvZywxBO~QobK1LJpwPpnoW9PZXJ~bjKv-knmzpF~dH35XN7R9VKBSv2FvGeKQ__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA" && \	
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
 apt-utils \
 unzip=6.0-23+deb10u1 \
 procps=2:3.3.15-2 \
 build-essential=12.6 \
 zlib1g-dev=1:1.2.11.dfsg-1 \
 libxt-dev \
 libmagick++-dev \
 python-pip

RUN pip install pandas

RUN R -e "install.packages('ggbeeswarm')" && \
	R -e "BiocManager::install('monocle')" && \
	R -e "devtools::install_version('Seurat', version = '3.1.5', dependencies= NA)" && \
	R -e "devtools::install_version('future', version = '1.17.0', dependencies= T)" && \
	R -e "devtools::install_version('cowplot', version = '1.0.0', dependencies= T)" && \
	R -e "devtools::install_version('clustree', version = '0.4.2', dependencies= NA)" && \
	R -e "devtools::install_version('gridExtra', version = '2.3', dependencies=T)" && \
	R -e "devtools::install_version('getopt', version = '1.20.3', dependencies=T)" && \
	R -e "devtools::install_version('pheatmap', version = '1.0.12', dependencies=T)" && \
	R -e "BiocManager::install('limma')" && \
	R -e "devtools::install_github('juliendelile/Antler', dependencies = TRUE)"