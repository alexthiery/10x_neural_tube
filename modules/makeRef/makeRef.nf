#!/usr/bin/env nextflow

nextflow.preview.dsl=2

process makeRef {
    label 'mid_memory'

    publishDir "${params.outdir}",
        mode: "copy", overwrite: true

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