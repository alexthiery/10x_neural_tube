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
2) Download chick genome ([galgal6](ftp://ftp.ensembl.org/pub/release-97/fasta/gallus_gallus/dna/Gallus_gallus.GRCg6a.dna.toplevel.fa.gz))
3) Download annotation gtf ([galgal6](ftp://ftp.ensembl.org/pub/release-97/gtf/gallus_gallus/Gallus_gallus.GRCg6a.97.gtf.gz))
2) Download raw reads from [here]()
3) Make a csv file containing the sample names and corresponding paths using this [template](sampleInfo.csv)
4) Run Nextflow using the following command

``` sh
nextflow run alexthiery/10x_neural_tube -r Dev \
-hub github \
-profile singularity \
--metadata /camp/home/thierya/scratch/10x_neural_tube/sampleInfo.csv \
--gtf /camp/home/thierya/working/genomes/galgal6/Gallus_gallus.GRCg6a.97.gtf \
--fa /camp/home/thierya/working/genomes/galgal6/Gallus_gallus.GRCg6a.dna.toplevel.fa
```
#
## Interactive downstream analysis
#

#
## Downstream analysis pipeline
#
