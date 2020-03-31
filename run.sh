#!/bin/bash

ml purge
ml Anaconda2
ml R/3.6.0-foss-2019a

sbatch -N 1 -c 8 --mem=100G -t 12:00:00 --wrap="Rscript /camp/home/thierya/working/alexthiery/analysis/10x_scRNAseq_2019/repo/scripts/1_seurat_full.R CAMP"