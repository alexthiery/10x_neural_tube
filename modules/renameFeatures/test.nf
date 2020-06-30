#!/usr/bin/env nextflow

nextflow.preview.dsl=2

include renameFeatures from "$baseDir/renameFeatures.nf"

Channel
    .fromPath( params.sampleDir )
    .set { ch_features }

gtf = params.gtf

workflow {
    renameFeatures( ch_features, gtf )
}