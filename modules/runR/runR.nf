#!/usr/bin/env nextflow

nextflow.preview.dsl=2

process runR {
    label 'high_memory'

    input:
        path(rFile)
        path(cellrangerOut)

    output:
        path("plots")
        path("RDS.files")

    """
    Rscript ${rFile} --customFuncs ${params.customFuncs} --networkGenes ${params.networkGenes} --countFiles ${cellrangerOut} --cores ${task.cpus} --runtype nextflow
    """
}