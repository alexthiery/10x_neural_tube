#!/usr/bin/env nextflow

nextflow.preview.dsl=2

params.rFile = "$baseDir/bin/R/seurat_full.R"
params.extraData = "$baseDir/bin/network_genes"
params.customFunctions = "$baseDir/bin/R/my_functions"

Channel
    .fromPath(params.rFile)
    .set { ch_rFile }

workflow {
    runR( ch_rFil )
}