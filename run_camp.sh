#!/bin/bash

#SBATCH --job-name=10x_neural_tube
#SBATCH --mail-type=ALL,ARRAY_TASKS
#SBATCH --mail-user=alex.thiery@crick.ac.uk

ml purge
ml Anaconda2
ml R/3.6.0-foss-2019a
source activate 10x_neural_tube
cd /camp/home/thierya/working/analysis/10x_neural_tube/repo
# Rscript /camp/home/thierya/working/analysis/10x_neural_tube/repo/scripts/1_seurat_full.R --location CAMP --cores 16
Rscript /camp/home/thierya/working/analysis/10x_neural_tube/repo/scripts/2_neural_subset.R --location CAMP --cores 16