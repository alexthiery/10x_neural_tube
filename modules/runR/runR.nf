#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process runR {
    publishDir "${params.outDir}",
        mode: "copy", overwrite: true

    label 'process_high'

    input:
        path(rFile)
        path(cellrangerOut)

    output:
        path("plots")
        path("RDS.files")
        path("antler.input")

    """
    Rscript ${rFile} --customFuncs ${params.customFuncs} --networkGenes ${params.networkGenes} --cores ${task.cpus} --runtype nextflow
    """
}