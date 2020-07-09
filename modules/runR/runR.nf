#!/usr/bin/env nextflow

nextflow.preview.dsl=2

process runR {
    publishDir "${params.r_outDir}/test",
        mode: "copy", overwrite: true

    label 'high_memory'

    input:
        path(rFile)
        path(cellrangerOut)

    output:
        path("plots")
        path("RDS.files")

    """
    Rscript ${rFile} --customFuncs ${params.customFuncs} --networkGenes ${params.networkGenes} --cores ${task.cpus} --runtype nextflow
    """
}