#!/usr/bin/env nextflow

params.inputfile = '$baseDir/bin/R/test.R'
params.myFunc = '$baseDir/bin/R/my_functions'

//rFile_ch = Channel.fromPath(params.inputfile)
//customFunctions_ch = Channel.fromPath(params.myFunc)

params.runName = '1_seurat_full'

process run_1_seurat_full {

    publishDir "${params.outDir}/${params.runName}",
        mode: "copy", overwrite: false

    input:
    //    path rFile from rFile_ch
      //  path customFunctions from customFunctions_ch

    output:
        path("plots")
        path("processed_data")

    """
    Rscript ${params.inputfile} ${params.myFunc}
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