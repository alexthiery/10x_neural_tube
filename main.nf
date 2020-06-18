#!/usr/bin/env nextflow
params.runName = '1_seurat_full'

params.rFile = "$baseDir/bin/R/1_seurat_full.R"
params.customFunctions = "$baseDir/bin/R/my_functions"

sample_ch = Channel.fromPath(params.sampleDir)

process run_1_seurat_full {
    
    cpus 16

    publishDir "${params.outDir}/${params.runName}",
        mode: "copy", overwrite: true

    input:
        path samples from sample_ch

    output:
        path("plots")
        path("processed_data")

    """
    Rscript ${params.rFile} ${params.customFunctions} ${samples} ${params.ncores}
    """
}

// params.runName = '2_seurat_full'

// process run_2_neural_subset {

//     publishDir "${params.outDir}/${params.runName}",
//         mode: "copy", overwrite: false

//     input:
//         path samples from sample_ch
//         path rFile from rFile_ch

//     output:
//         path("plots")
//         path("processed_data")

//     """
//     Rscript ${rFile} ${samples}
//     """
// }