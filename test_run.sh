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

nextflow pull alexthiery/10x_neural_tube -r nextflow -hub github

nextflow run alexthiery/10x_neural_tube -r nextflow \
-hub github \
-profile crick \
--metadata /camp/home/thierya/scratch/10x_neural_tube/sampleInfo.csv \
--gtf /camp/home/thierya/working/genomes/galgal6/Gallus_gallus.GRCg6a.97.gtf \
--fa /camp/home/thierya/working/genomes/galgal6/Gallus_gallus.GRCg6a.dna.toplevel.fa \
-resume \
-N alex.thiery@crick.ac.uk


# nextflow run alexthiery/10x_neural_tube -r nextflow \
# -hub github \
# -profile crick \
# --metadata /camp/home/thierya/scratch/10x_neural_tube/sampleInfo.csv \
# --ref /camp/svc/scratch/luscomben/thierya/10x_neural_tube/null/reference_genome \
# -N alex.thiery@crick.ac.uk

# nextflow run test.nf \
# -with-docker alexthiery/10x_neural_tube:dev \
# --metadata /Users/alex/dev/repos/10x_neural_tube/sampleInfo.csv \
# --ref /Users/alex/dev/repos/10x_neural_tube/nf-testData/ \
# --cpus 4 \
# --ram 16