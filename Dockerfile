FROM rocker/tidyverse:3.6.3

LABEL authors="alex.thiery@crick.ac.uk" \
      description="Docker image containing all requirements to run 10x analysis"

ARG WHEN


# Install cellranger
RUN cd /tmp/ && \
	wget -O cellranger-3.0.2.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-3.0.2.tar.gz?Expires=1595277739&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cHM6Ly9jZi4xMHhnZW5vbWljcy5jb20vcmVsZWFzZXMvY2VsbC1leHAvY2VsbHJhbmdlci0zLjAuMi50YXIuZ3oiLCJDb25kaXRpb24iOnsiRGF0ZUxlc3NUaGFuIjp7IkFXUzpFcG9jaFRpbWUiOjE1OTUyNzc3Mzl9fX1dfQ__&Signature=LMn2xLIs5qxoNlchuWrpiBMG7rKHj9VQGDG6mL0-8Wj2TF2~vXw8PoOd8TCJQ9E8IELmM~B-8OoeCKx1LTsH--KLEfvndsxX~WXOsZDsDO3KKZMHRKBx8Be2yqmedaVtk4BOYkPpjRKmMOXgXWt0xrOJpYy6yWvUgdeUW~t5EBP~P22hJ-FPFcbmG~t4aqNLzVmqVC6D1RQaaifR2LQA1NTdVvMiG~t-itLC3Bg4RMECweT6ascHNpZu~eXpQtqlnxthsw-cC7MNaPawDZR9y27E1acC1aEVyq01gCLOFy8aagm9XyY8EDpvdizGgKDTp1qGiTxb2rbhOgkH9PEmnA__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA" && \	
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
 libmagick++-dev \
 python-pip

RUN pip install pandas

RUN R -e "devtools::install_version('Seurat', version = '3.1.5', dependencies= NA)"
RUN R -e "devtools::install_version('future', version = '1.17.0', dependencies= T)"
RUN R -e "devtools::install_version('cowplot', version = '1.0.0', dependencies= T)"
RUN R -e "devtools::install_version('clustree', version = '0.4.2', dependencies= NA)"
RUN R -e "devtools::install_version('gridExtra', version = '2.3', dependencies=T)"
RUN R -e "devtools::install_version('getopt', version = '1.20.3', dependencies=T)"
RUN R -e "devtools::install_version('pheatmap', version = '1.0.12', dependencies=T)"
RUN R -e "BiocManager::install('limma')"