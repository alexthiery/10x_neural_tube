#!/usr/bin/env nextflow

nextflow.preview.dsl=2

/*-----------------------------------------------------------------------------------------------------------------------------
Include modules
-------------------------------------------------------------------------------------------------------------------------------*/

include projectHeader from "$baseDir/modules/projectHeader/projectHeader.nf"
include filterGTF from "$baseDir/modules/filterGTF/filterGTF.nf"
include makeRef from "$baseDir/modules/makeRef/makeRef.nf"
include cellrangerCount from "$baseDir/modules/cellrangerCount/cellrangerCount.nf"


/*-----------------------------------------------------------------------------------------------------------------------------
Pipeline params
-------------------------------------------------------------------------------------------------------------------------------*/
params.outgenomename = "galgal6_filtered_ref_genome"
params.outdir = "./temp_out"

testFastq = [['THI300A1', 'cellranger_count_hh4', params.testfile]]

// Check params
if (!params.gtf) {
    exit 1, "No gtf file provided"
}
if (!params.fa) {
    exit 1, "No fa file provided"
}
if (!params.cpus) {
    params.cpus = 2
}
if (!params.ram) {
    params.ram = 4
}

/*-----------------------------------------------------------------------------------------------------------------------------
Init
-------------------------------------------------------------------------------------------------------------------------------*/
// Show banner
log.info projectHeader()


/*------------------------------------------------------------------------------------*/
/* Define input channels
--------------------------------------------------------------------------------------*/

Channel
  .from(params.gtf)
  .set {ch_gtf}

Channel
  .from(params.fa)
  .set {ch_fa}

Channel
  .from(testFastq)
  .set {ch_test_fastq}

/*-----------------------------------------------------------------------------------------------------------------------------
Main workflow
-------------------------------------------------------------------------------------------------------------------------------*/

workflow {
    filterGTF(ch_gtf)
    makeRef(filterGTF.out, ch_fa)
    cellrangerCount(ch_test_fastq, makeRef.out)
}

