#!/bin/bash
#SBATCH -t 2:00:00
#SBATCH --mem=64gb
#SBATCH --cpus-per-task=16
#SBATCH -N 1
#SBATCH --partition=cpu
#SBATCH --job-name=10x_neural_tube
#SBATCH --mail-type=ALL,ARRAY_TASKS
#SBATCH --mail-user=alex.thiery@crick.ac.uk

ml purge
ml Anaconda2
ml R/3.6.0-foss-2019a
source activate 10x_neural_tube
cd /camp/home/thierya/working/analysis/10x_neural_tube/repo
Rscript /camp/home/thierya/working/analysis/10x_neural_tube/repo/scripts/1_seurat_full.R --location CAMP --cores 16