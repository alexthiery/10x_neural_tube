#!/usr/bin/env nextflow

nextflow.preview.dsl=2

include renameFeatures from "$baseDir/renameFeatures.nf"
// params.outdir = "./output2"

params.sampleDir = "/Users/alex/dev/repos/10x_neural_tube/testData/*"
params.gtf = "/Users/alex/dev/genomes/galgal6/Gallus_gallus.GRCg6a.97.gtf"


Channel
    .fromPath( params.sampleDir, type: 'dir'  )
    .view()
    .set { ch_features }

Channel
    .fromPath( params.gtf)
    .view()
    .set { ch_gtf }


workflow {
    renameFeatures( ch_features.combine(ch_gtf) )
}