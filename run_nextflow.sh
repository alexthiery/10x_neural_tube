#!/bin/sh
#SBATCH -t 2:00:00
#SBATCH --mem=64gb
#SBATCH --cpus-per-task=16
#SBATCH -N 1
#SBATCH --partition=cpu
#SBATCH --job-name=10x_neural_tube
#SBATCH --mail-type=ALL,ARRAY_TASKS
#SBATCH --mail-user=alex.thiery@crick.ac.uk
export TERM=xterm

## LOAD REQUIRED MODULES
ml purge
ml Nextflow/20.01.0
ml Singularity/3.4.2
ml Graphviz

## UPDATE PIPLINE
nextflow pull alexthiery/10x_neural_tube -r keep_sex_genes -hub github

## RUN PIPELINE

nextflow run alexthiery/10x_neural_tube -r keep_sex_genes \
-hub github \
-profile crick \
--sampleDir /camp/home/thierya/working/analysis/10x_neural_tube/cellranger_ouput \
--email alex.thiery@crick.ac.uk