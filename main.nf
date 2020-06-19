#!/usr/bin/env nextflow

params.rFile = "$baseDir/bin/R/1_seurat_full_nosexfilt.R"
params.customFunctions = "$baseDir/bin/R/my_functions"

sample_ch = Channel.fromPath(params.sampleDir)

process run_1_seurat_full_nosexfilt {
    
    // cpus determined by params in profile configs
    cpus 10

    // split params.rFile and keep only the extension as name for outDir
    bits = params.rFile.take(params.rFile.lastIndexOf('.')).split("/")

    publishDir "${params.outDir}/${bits[bits.length-1]}",
        mode: "copy", overwrite: true

    input:
        path samples from sample_ch

    output:
        path("plots")
        path("RDS.files")

    """
    Rscript ${params.rFile} --myfuncs ${params.customFunctions} --samples ${samples} --cores ${task.cpus} --location CAMP
    """
}
