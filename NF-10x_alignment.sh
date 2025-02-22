#!/bin/bash
#SBATCH --job-name=10x
#SBATCH -t 72:00:00
#SBATCH --mail-type=ALL,ARRAY_TASKS
#SBATCH --mail-user=alex.thiery@crick.ac.uk

## LOAD REQUIRED MODULES
ml purge
ml Nextflow/20.07.1
ml Singularity/3.4.2
ml Graphviz

export TERM=xterm
export NXF_VER=20.07.1
export NXF_SINGULARITY_CACHEDIR=/camp/home/thierya/working/NF_singularity

nextflow run NF-10x_alignment/main.nf \
-c conf/crick.config \
--samplesheet NF-10x_alignment/crick_samplesheet.csv \
--output output/NF-10x_alignment \
--gtf /camp/home/thierya/working/genomes/galgal5/Gallus_gallus.Gallus_gallus-5.0.94.gtf \
--fasta /camp/home/thierya/working/genomes/galgal5/Gallus_gallus.Gallus_gallus-5.0.dna.toplevel.fa \
-resume