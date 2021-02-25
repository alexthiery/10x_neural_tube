#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process modifyGTF {

    publishDir "${params.output}/modifiedGTF",
    mode: "copy", overwrite: true

    label 'process_low'

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