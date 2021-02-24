#!/bin/bash
#SBATCH --job-name=10x
#SBATCH -t 72:00:00
#SBATCH --mail-type=ALL,ARRAY_TASKS
#SBATCH --mail-user=alex.thiery@crick.ac.uk

export TERM=xterm

## LOAD REQUIRED MODULES
ml purge
ml Nextflow/20.07.1
ml Singularity/3.4.2
ml Graphviz

nextflow pull alexthiery/10x_neural_tube -hub github

NXF_VER=20.07.1

nextflow run NF-10x_alignment/main.nf \
-profile crick \
--samplesheet ./NF-10x_alignment/crick_samplesheet.csv \
--gtf /camp/home/thierya/working/genomes/galgal5/Gallus_gallus.Gallus_gallus-5.0.94_modified.gtf \
--fasta /camp/home/thierya/working/genomes/galgal5/Gallus_gallus.Gallus_gallus-5.0.dna.toplevel.fa \
-resume


nextflow run NF-10x_alignment/main.nf \
-profile local \
--samplesheet ./NF-10x_alignment/crick_samplesheet.csv \
--gtf ~/dev/genomes/galgal5/Gallus_gallus.Gallus_gallus-5.0.94_modified.gtf \
--fasta ~/dev/genomes/galgal5/Gallus_gallus.Gallus_gallus-5.0.dna.toplevel.fa \
-resume