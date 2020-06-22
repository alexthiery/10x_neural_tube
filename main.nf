#!/usr/bin/env nextflow

nextflow.preview.dsl=2

/*-----------------------------------------------------------------------------------------------------------------------------
Include modules
-------------------------------------------------------------------------------------------------------------------------------*/

include projectHeader from "$baseDir/modules/projectHeader.nf"
include filterGTF from "$baseDir/modules/filterGTF.nf"
include makeRef from "$baseDir/modules/makeRef.nf"
include cellrangerCount from "$baseDir/modules/cellrangerCount.nf"


/*-----------------------------------------------------------------------------------------------------------------------------
Pipeline params
-------------------------------------------------------------------------------------------------------------------------------*/

params.outgenomename = "galgal6_filtered_ref_genome"
params.outdir = "./temp_out"


/*-----------------------------------------------------------------------------------------------------------------------------
Init
-------------------------------------------------------------------------------------------------------------------------------*/
// Show banner
log.info projectHeader()

// Check params
if (!params.cpus) {
    params.cpus = 2
}
if (!params.ram) {
    params.ram = 4
}

/*-----------------------------------------------------------------------------------------------------------------------------
Main workflow
-------------------------------------------------------------------------------------------------------------------------------*/
workflow {
    filterGTF(params.gtf)
    makeRef(filterGTF.out, params.fa)
    cellrangerCount(insert_tuple, makeRef.out)
}