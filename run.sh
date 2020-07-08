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

nextflow pull alexthiery/10x_neural_tube -r master -hub github

nextflow run alexthiery/10x_neural_tube -r master \
-hub github \
-profile crick \
--metadata /camp/home/thierya/scratch/10x_neural_tube/sampleInfo.csv \
--gtf /camp/home/thierya/working/genomes/galgal6/Gallus_gallus.GRCg6a.97.gtf \
--fa /camp/home/thierya/working/genomes/galgal6/Gallus_gallus.GRCg6a.dna.toplevel.fa


# NXF_VER=20.01.0 nextflow run alexthiery/10x_neural_tube -r Dev \
# -hub github \
# -profile custom \
# --metadata /camp/home/thierya/scratch/10x_neural_tube/sampleInfo.csv \
# --gtf ~/dev/genomes/galgal6/Gallus_gallus.GRCg6a.97.gtf \
# --fa ~/dev/genomes/galgal6/Gallus_gallus.GRCg6a.dna.toplevel.fa 