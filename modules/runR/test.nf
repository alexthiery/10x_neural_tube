#!/usr/bin/env nextflow

nextflow.preview.dsl=2

params.rFile = "/Users/alex/dev/repos/10x_neural_tube/bin/R/seurat_full.R"
params.customFuncs = "/Users/alex/dev/repos/10x_neural_tube/bin/R/custom_functions"
params.networkGenes = "/Users/alex/dev/repos/10x_neural_tube/bin/network_genes"

include runR from "/Users/alex/dev/repos/10x_neural_tube/modules/runR/runR.nf"

Channel
    .fromPath(params.rFile)
    .set { ch_rFile }
    
workflow {
    runR( ch_rFile, params.cellrangerOut )
}



// nextflow run ./modules/runR/test.nf \
// -with-docker alexthiery/10x_neural_tube:dev \
// --rFile ~/dev/repos/10x_neural_tube/bin/R/seurat_full.R \
// --customFuncs ~/dev/repos/10x_neural_tube/bin/R/custom_functions \
// --networkGenes ~/dev/repos/10x_neural_tube/bin/network_genes \
// --cellrangerOut ~/dev/repos/10x_neural_tube/results/cellrangerCounts