# Chick 10x Neural Development

<p align="center">

<img src="./suppl_files/TissueCollection.png" width="100%">

</p>

This repository provides the code to run the 10x single cell RNA-eq analysis


- [Data availability](#data-availability)
- [Analysis pre-requisites](#analysis-pre-requisites)
- [Docker](#docker)
- [Nextflow](#nextflow)
- [Interactive downstream analysis](#interactive-downstream-analysis)
- [Downstream analysis pipeline](#downstream-analysis-pipeline)

#
## Data availability
#
This repository contains the required code to run the entire alignment and downstream analysis pipeline. The pipeline is run using [Nextflow](#nextflow) and [Docker](#docker) to ensure reproducibility. For those who would like to re-run the analysis, the following files should be downloaded:

- to run the analysis including the alignment, the raw fastq sequencing files can be found [here]()
- to run the downstream analysis from the UMI counts generated from 10x Genomics Cell Ranger are embedded can be found within the repository [here]("./alignmentOut/cellrangerCounts_renamed")

#
## Analysis pre-requisites
#
The repository can be called directly from GitHub, so to re-run the analysis you just need to install [Nextflow]('https://www.nextflow.io/docs/latest/getstarted.html#installation') and [Docker]('https://www.docker.com/').


#
## Docker
#

#
## Nextflow
#

#
## Interactive downstream analysis
#

#
## Downstream analysis pipeline
#
