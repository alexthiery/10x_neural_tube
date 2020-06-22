#!/usr/bin/env nextflow

nextflow.preview.dsl=2

process makeRef {
    cpus = params.cpus

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
    --fasta=${fasta} \
    --genes=${filt_genome} \
    --nthreads=${params.cpus} \
    --memgb=${params.ram}

    """
}