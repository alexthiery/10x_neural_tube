#!/usr/bin/env nextflow

// Define DSL2
nextflow.preview.dsl=2

// Log
log.info ("Starting test pipeline for makeRef")

/*------------------------------------------------------------------------------------*/
/* Module inclusions 
--------------------------------------------------------------------------------------*/

include makeRef from "$baseDir/makeRef.nf"

/*------------------------------------------------------------------------------------*/
/* Params
--------------------------------------------------------------------------------------*/
params.outgenomename = "galgal6_filtered_ref_genome"
params.outdir = "./temp_out"

// Check params
if (!params.cpus) {
    params.cpus = 2
}
println params.cpus
if (!params.ram) {
    params.ram = 4
}

fa = params.fa
gtf = params.gtf

/*------------------------------------------------------------------------------------*/
/* Define input channels
--------------------------------------------------------------------------------------*/
Channel
  .from(fa)
  .set {ch_fa}

Channel
  .from(gtf)
  .set {ch_gtf}

workflow {
    makeRef(ch_gtf, ch_fa)
}