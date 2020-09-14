#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process modifyGTF {

    publishDir "${params.alignment_outDir}/modifiedGTF",
    mode: "copy", overwrite: true

    label 'low_memory'

    container 'alexthiery/10x_neural_tube:v1.0'


    input:
        path(pythonFile)
        path(gtf)
        path(geneNames)

    output:
        path('./modified.gtf')

    """
    #!/bin/bash
    python ${pythonFile} ${gtf} ${geneNames} ./modified.gtf

    """
}