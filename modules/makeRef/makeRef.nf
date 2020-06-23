#!/usr/bin/env nextflow

nextflow.preview.dsl=2

process makeRef {
    cpus = params.cpus

    publishDir "${params.outdir}",
        mode: "copy", overwrite: true

    input:
        path("$baseDir/"filt_genome)
        path("$baseDir/"fasta)

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