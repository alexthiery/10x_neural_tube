#!/bin/bash
#SBATCH --job-name=10x
#SBATCH -t 72:00:00
#SBATCH --mail-type=ALL,ARRAY_TASKS
#SBATCH --mail-user=alex.thiery@crick.ac.uk

export TERM=xterm

## LOAD REQUIRED MODULES
ml purge
ml Nextflow/20.01.0
ml Singularity/3.4.2
ml Graphviz

nextflow pull alexthiery/10x_neural_tube -r R -hub github

NXF_VER=20.01.0 nextflow run alexthiery/10x_neural_tube -r R \
-hub github \
-profile crick \
--counts /camp/home/thierya/scratch/10x_neural_tube/storedat/alignmentOut/cellrangerCounts_renamed