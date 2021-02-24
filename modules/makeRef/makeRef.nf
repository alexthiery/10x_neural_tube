#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process makeRef {
    label 'process_medium'

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