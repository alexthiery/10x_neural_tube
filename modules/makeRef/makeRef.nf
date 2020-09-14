#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process makeRef {
    label 'mid_memory'

    container 'alexthiery/10x_neural_tube:v1.0'

    input:
        path(filt_genome)
        path(fasta)

    output:
        path("reference_genome")

    """
    #!/bin/bash

    # make reference
    cellranger mkref --genome=reference_genome \
    --genes=${filt_genome} \
    --fasta=${fasta} \
    --nthreads=${task.cpus}

    """
}