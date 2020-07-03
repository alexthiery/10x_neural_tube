#!/usr/bin/env nextflow

process runR {
    label 'high_memory'

    // split params.rFile and keep only the extension as name for outDir
    bits = params.rFile.take(params.rFile.lastIndexOf('.')).split("/")

    publishDir "${params.outDir}/${bits[bits.length-1]}",
        mode: "copy", overwrite: true

    input:
        path(rFile)
        path(cellrangerOut)

    output:
        path("plots")
        path("RDS.files")

    """
    Rscript ${rFile} --myfuncs ${params.customFunctions} --extraData ${params.extraDat} --countFiles ${cellrangerOut} --cores ${task.cpus} --runtype nextflow
    """
}