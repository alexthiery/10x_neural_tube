#!/bin/sh

export TERM=xterm

## LOAD REQUIRED MODULES
ml purge
ml Nextflow/20.01.0
ml Singularity/3.4.2
ml Graphviz


nextflow pull alexthiery/10x_neural_tube -r nextflow -hub github

nextflow run alexthiery/10x_neural_tube -r nextflow \
-hub github \
-profile docker \
--gtf /camp/home/thierya/working/genomes/galgal6/Gallus_gallus.GRCg6a.97.gtf \
--fa /camp/home/thierya/working/genomes/galgal6/Gallus_gallus.GRCg6a.dna.toplevel.fa \
--threads 16 \
--mem 64 \
-with-report report.html \
-N alex.thiery@crick.ac.uk



nextflow run main.nf \
-with-docker alexthiery/10x_neural_tube:dev \
--gtf ~/dev/genomes/test.gtf \
--fa ~/dev/genomes/galgal6/Gallus_gallus.GRCg6a.dna.toplevel.fa \
--cpus 4 \
--ram 8



