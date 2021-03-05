#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process runR {
    publishDir "${params.output}",
        mode: "copy", overwrite: true

    label 'process_high'

    input:
        path rFile
        path 'input/*'

    output:
        file '*'

    """
    Rscript ${rFile} --customFuncs ${params.customFuncs} --networkGenes ${params.networkGenes} --cores ${task.cpus} --runtype nextflow
    rm -r input
    """
}