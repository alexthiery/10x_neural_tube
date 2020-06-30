#!/usr/bin/env nextflow

process run_seurat_full {
    label 'high_memory'

    // split params.rFile and keep only the extension as name for outDir
    bits = params.rFile.take(params.rFile.lastIndexOf('.')).split("/")

    publishDir "${params.outDir}/${bits[bits.length-1]}",
        mode: "copy", overwrite: true

    input:
        path samples from sample_ch
        path extraDat from extraData_ch

    output:
        path("plots")
        path("RDS.files")

    """
    Rscript ${params.rFile} --myfuncs ${params.customFunctions} --extraData ${extraDat} --samples ${samples} --cores ${task.cpus} --location CAMP
    """
}