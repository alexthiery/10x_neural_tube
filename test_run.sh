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
-profile crick \
--metadata /camp/home/thierya/scratch/10x_neural_tube/sampleInfo.csv
--gtf /camp/home/thierya/working/genomes/galgal6/Gallus_gallus.GRCg6a.97.gtf \
--fa /camp/home/thierya/working/genomes/galgal6/Gallus_gallus.GRCg6a.dna.toplevel.fa \
--cpus 16 \
--ram 64 \
-with-report report.html \
-N alex.thiery@crick.ac.uk

# nextflow run alexthiery/10x_neural_tube -r nextflow \
# -hub github \
# -profile docker \
# --folder1 /Users/alex/dev/repos/10x_neural_tube/nf-testData/cellrangerCount \
# --folder2 /Users/alex/dev/repos/10x_neural_tube/nf-testData/cellrangerCount2 \
# --gtf /Users/alex/dev/genomes/galgal6/Gallus_gallus.GRCg6a.97.gtf \
# --fa /Users/alex/dev/genomes/galgal6//Gallus_gallus.GRCg6a.dna.toplevel.fa \
# -with-report report.html \
# -N alex.thiery@crick.ac.uk

