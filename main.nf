#!/usr/bin/env nextflow

rFile_ch = Channel.fromPath('./bin/R/test.R')
customFunctions_ch = Channel.fromPath('./bin/R/my_functions')

process run_1_seurat_full {

    publishDir "${params.outDir}/${params.runName}",
        mode: "copy", overwrite: false

    input:
        //path samples from sample_ch
        path rFile from rFile_ch
        path customFunctions from customFunctions_ch

    output:
        path("plots")
        path("processed_data")

    """
    Rscript ${rFile} ${customFunctions}
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