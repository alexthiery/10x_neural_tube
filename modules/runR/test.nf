#!/usr/bin/env nextflow

nextflow.preview.dsl=2

// params.rFile = "$baseDir/bin/R/seurat_full.R"
// params.customFunctions = "$baseDir/bin/R/my_functions"
// params.extraData = "$baseDir/bin/network_genes"

Channel
    .fromPath(params.rFile)
    .set { ch_rFile }

process runR {
    label 'high_memory'

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


workflow {
    runR( ch_rFile, params.sampleDir )
}



// nextflow run ./modules/runR/test.nf -with-docker alexthiery/10x_neural_tube:dev --rFile ./bin/R/seurat_full.R --customFunctions ./bin/R/my_functions --extraData ./bin/network_genes --sampleDir ~/dev/repos/10x_neural_tube/results/cellrangerCounts