# Chick 10x Neural Development

<p align="center">

<img src="./suppl_files/TissueCollection.png" width="100%">

</p>

This repository provides the code to run the 10x single cell RNA-eq analysis.


- [Data availability](#data-availability)
- [Analysis pre-requisites](#analysis-pre-requisites)
- [Docker](#docker)
- [Nextflow](#nextflow)
- [Interactive downstream analysis](#interactive-downstream-analysis)
- [Downstream analysis pipeline](#downstream-analysis-pipeline)

#
## Data availability
#
This repository contains the required code to run the entire alignment and downstream analysis pipeline. For those who would like to re-run the analysis, the following files should be downloaded:

- to run the analysis including the alignment, the raw fastq sequencing files can be found [here]().
- to run the downstream analysis from the UMI counts generated from 10x Genomics Cell Ranger are embedded can be found within the repository [here]("./alignmentOut/cellrangerCounts_renamed").

#
## Analysis pre-requisites
#
The pipeline is run using Nextflow and Docker to ensure reproducibility. The repository can be called directly from GitHub, so to re-run the analysis you just need to install [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html#installation) and [Docker](https://docs.docker.com/get-docker/).

If you are wanting to run the downstream analysis interactively outside of Nextflow, you still need to download Docker and you will also need to download this repository.

#
## Docker
#
Docker allows us to run our analysis through a container with all required libraries and dependencies. This ensures cross-platform reproducibility of our results.

The docker image used can be found [here](https://hub.docker.com/r/alexthiery/10x_neural_tube)

We have also included Rstudio within our docker image to allow you to run the downstream analysis interactively if you wish - for details on how to do this click [here](#interactive-downstream-analysis).

#
## Nextflow
#
You can easily re-run our entire pipeline in Nextflow using the following steps:
1) Install Nextflow and Docker
2) Download chick genome ([galgal6](http://ftp.ensembl.org/pub/release-97/fasta/gallus_gallus/dna/Gallus_gallus.GRCg6a.dna.toplevel.fa.gz))
3) Download annotation gtf ([galgal6](http://ftp.ensembl.org/pub/release-97/gtf/gallus_gallus/Gallus_gallus.GRCg6a.97.gtf.gz))
2) Download raw reads from [here]()
3) Make a sampleInfo.csv file containing the sample names and corresponding paths using this [template](sampleInfo.csv)
4) Run Nextflow using the following command

``` sh
nextflow run alexthiery/10x_neural_tube \
-profile docker \
--metadata <path to sampleInfo.csv> \
--gtf <path to gtf> \
--fa <path to genome>
```

This pipeline is configured to be ran on a cluster with 224GB memory and 32CPUs by default. The `-profile` flag can be used to set  either 'docker' or 'singularity', depending on the container application installed on your system. These settings can be adjusted  by replacing the `-profile` flag with a custom config file as below.

``` sh
nextflow run alexthiery/10x_neural_tube \
-c <path to custom.config file> \
--metadata <path to sampleInfo.csv> \
--gtf <path to gtf> \
--fa <path to genome>
```

 For a template of a custom.config file, see [crick.config](conf/crick.config). Further information on Nextflow config files can be found [here](https://www.nextflow.io/docs/latest/config.html#configuration-file).



#
## Interactive downstream analysis
#

If do not want to re-run the alignment, but would like to run the downstream analysis from the count files, you can run RStudio from within the docker container.

To do this, follow these steps:

1) clone this repository to your local computer
2) start a terminal session and download the docker image - `docker pull alexthiery/10x_neural_tube:v1.0`
3) within terminal launch a docker container interactively - `docker run --rm -e PASSWORD=test -p 8787:8787 -v <PATH_TO_LOCAL_REPOSITORY>:/home/rstudio alexthiery/10x_neural_tube:v1.0`
4) go to `localhost:8787` in the browser to open RStudio

#
## Downstream analysis pipeline
#

This analysis is ran using Seurat v3.1. For more information see https://satijalab.org/seurat/


Set environmental variables, load data and make Seurat object
In order to be able to run the script from either Rstudio, local terminal, or cluster terminal, I add a switch which looks for command line arguments. This then sets the directory paths accordingly.

``` R
#!/usr/bin/env Rscript

library('getopt')

# set arguments for Rscript
spec = matrix(c(
  'runtype', 'l', 2, "character",
  'cores'   , 'c', 2, "integer",
  'customFuncs', 'm', 2, "character",
  'networkGenes', 'd', 2, "character"
), byrow=TRUE, ncol=4)
opt = getopt(spec)

# set default location
if(length(commandArgs(trailingOnly = TRUE)) == 0){
  cat('No command line arguments provided, user defaults paths are set for running interactively in Rstudio on docker\n')
  opt$runtype = "user"
} else {
  if(is.null(opt$runtype)){
    stop("--runtype must be either 'user', 'nextflow' or 'docker'")
  }
  if(tolower(opt$runtype) != "docker" & tolower(opt$runtype) != "user" & tolower(opt$runtype) != "nextflow"){
    stop("--runtype must be either 'user', 'nextflow' or 'docker'")
  }
  if(tolower(opt$runtype) == "nextflow"){
    if(is.null(opt$customFuncs)){
      stop("--customFuncs path must be specified")
    }
    if(is.null(opt$networkGenes)){
      stop("--networkGenes path must be specified")
    }
  }
}

####################################################################
# user paths need to be defined here in order to run interactively #
####################################################################
if (opt$runtype == "user"){
  sapply(list.files('./bin/R/custom_functions/', full.names = T), source)
  
  plot.path = "./results/R_results/plots/"
  rds.path = "./results/R_results/RDS.files/"
  dir.create(plot.path, recursive = T)
  dir.create(rds.path, recursive = T)
  
  ##################################
  # set path where data is located #
  ##################################
  data_path = "./alignmentOut/cellrangerCounts_renamed"
  
  # read all files from dir
  files <- list.files(data_path, recursive = T, full.names = T)
  # remove file suffix
  file.path <- dirname(files)[!duplicated(dirname(files))]
  # make dataframe with tissue matching directory
  tissue = c("hh4", "hh6", "ss4", "ss8")
  matches <- sapply(tissue, function(x) file.path[grep(pattern = x, x = file.path)])
  sample.paths <- data.frame(tissue = names(matches), path = matches, row.names = NULL)
  
  # Read in favourite genes
  network_genes <- list.files("./bin/network_genes/", full.names = T)
  hh4_genes <- read.table(network_genes[grepl("HH4", network_genes)], stringsAsFactors = F)[,1]
  hh6_genes <- read.table(network_genes[grepl("HH6", network_genes)], stringsAsFactors = F)[,1]

} else if (opt$runtype == "nextflow"){
  cat('pipeling running through nextflow\n')
  
  sapply(list.files(opt$customFuncs, full.names = T), source)
  
  plot.path = "plots/"
  dir.create(plot.path, recursive = T)
  rds.path = "RDS.files/"
  dir.create(rds.path, recursive = T)
  
  # read all files from folder and keep only those from chr_edit
  files <- list.files("./", recursive = T, full.names = T)
  # remove file suffix
  file.path <- dirname(files)[!duplicated(dirname(files))]
  # make dataframe with tissue matching directory
  tissue = c("hh4", "hh6", "ss4", "ss8")
  matches <- sapply(tissue, function(x) file.path[grep(pattern = x, x = file.path)])
  sample.paths <- data.frame(tissue = names(matches), path = matches, row.names = NULL)
  
  # Read in favourite genes
  network_genes <- list.files(opt$networkGenes, full.names = T)
  hh4_genes <- read.table(network_genes[grepl("HH4", network_genes)], stringsAsFactors = F)[,1]
  hh6_genes <- read.table(network_genes[grepl("HH6", network_genes)], stringsAsFactors = F)[,1]

} else if (opt$runtype == "docker"){
  cat('R script running through docker\n')
  
  sapply(list.files('/home/bin/R/custom_functions/', full.names = T), source)
  
  plot.path = "/home/results/R_results/plots/"
  rds.path = "/home/results/R_results/RDS.files/"
  dir.create(plot.path, recursive = T)
  dir.create(rds.path, recursive = T)
  
  data_path = "/home/alignmentOut/cellrangerCounts_renamed"
  # read all files from dir
  files <- list.files(data_path, recursive = T, full.names = T)
  # remove file suffix
  file.path <- dirname(files)[!duplicated(dirname(files))]
  # make dataframe with tissue matching directory
  tissue = c("hh4", "hh6", "ss4", "ss8")
  matches <- sapply(tissue, function(x) file.path[grep(pattern = x, x = file.path)])
  sample.paths <- data.frame(tissue = names(matches), path = matches, row.names = NULL)
  
  # Read in favourite genes
  network_genes <- list.files("/home/bin/network_genes/", full.names = T)
  hh4_genes <- read.table(network_genes[grepl("HH4", network_genes)], stringsAsFactors = F)[,1]
  hh6_genes <- read.table(network_genes[grepl("HH6", network_genes)], stringsAsFactors = F)[,1]
}

# set number of cores to use for parallelisation
if(is.null(opt$cores)){ncores = 4}else{ncores= opt$cores}
cat(paste0("script ran with ", ncores, " cores\n"))
```

Load packages

``` R
reticulate::use_python('/usr/bin/python3.7')
library(Seurat)

library(future)
library(dplyr)
library(cowplot)
library(clustree)
library(gridExtra)
library(grid)
library(pheatmap)
library(RColorBrewer)
```
