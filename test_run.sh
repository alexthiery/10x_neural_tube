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
--cpus 16 \
--ram 64
#--testfile /camp/home/thierya/working/raw_data/10x_2019/scRNA_Alex_raw/SC19040/190524_K00102_0343_AH7MYHBBXY/fastq/THI300A1 \



# nextflow run main.nf \
# -with-docker alexthiery/10x_neural_tube:dev \
# --gtf ~/dev/genomes/galgal6/Gallus_gallus.GRCg6a.97.gtf \
# --fa ~/dev/genomes/galgal6/Gallus_gallus.GRCg6a.dna.toplevel.fa \
# --cpus 2 \
# --ram 4 \
# -N alex.thiery@crick.ac.uk


