#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process filterGTF {
    label 'low_memory'

    input:
        path(gtf)

    output:
        path("filtered_genome.gtf")

    """
    #!/bin/bash

    # this step filters out genes based on the gene biotypes listed in attributes.
    cellranger mkgtf ${gtf} filtered_genome.gtf \
    --attribute=gene_biotype:protein_coding \
    --attribute=gene_biotype:lncRNA \
    --attribute=gene_biotype:pseudogene
    """
}