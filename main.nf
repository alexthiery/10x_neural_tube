#!/usr/bin/env nextflow

nextflow.preview.dsl=2

/*-----------------------------------------------------------------------------------------------------------------------------
Pipeline params
-------------------------------------------------------------------------------------------------------------------------------*/
params.rFile = "$baseDir/bin/R/seurat_full.R"
params.customFunctions = "$baseDir/bin/R/my_functions"
params.extraData = "$baseDir/bin/network_genes"
// // Check params
// if (!params.gtf) {
//     exit 1, "No gtf file provided"
// }
// if (!params.fa) {
//     exit 1, "No fa file provided"
// }
// if (!params.cpus) {
//     params.cpus = 2
// }
// if (!params.ram) {
//     params.ram = 4
// }

/*-----------------------------------------------------------------------------------------------------------------------------
Include modules
-------------------------------------------------------------------------------------------------------------------------------*/

include projectHeader from "$baseDir/modules/projectHeader/projectHeader.nf"
include filterGTF from "$baseDir/modules/filterGTF/filterGTF.nf"
include makeRef from "$baseDir/modules/makeRef/makeRef.nf"
include cellrangerCount from "$baseDir/modules/cellrangerCount/cellrangerCount.nf"
// include renameFeatures from "$baseDir/modules/renameFeatures/renameFeatures.nf"
// include runR from "$baseDir/modules/runR/runR.nf"

/*-----------------------------------------------------------------------------------------------------------------------------
Init
-------------------------------------------------------------------------------------------------------------------------------*/
// Show banner
log.info projectHeader()

/*------------------------------------------------------------------------------------*/
/* Define input channels
--------------------------------------------------------------------------------------*/

Channel
    .from( params.gtf )
    .set {ch_gtf}

Channel
    .from( params.fa )
    .set {ch_fa}

Channel
    .fromPath( params.metadata )
    .splitCsv(header: ['sample_id', 'sample_name', 'dir1', 'dir2'], skip: 1 )
    .map { row -> [row.sample_id, row.sample_name, file(row.dir1), file(row.dir2)] }
    .set { ch_fastq }

// Channel
//     .fromPath(params.rFile)
//     .set { ch_rFile }

/*-----------------------------------------------------------------------------------------------------------------------------
Main workflow
-------------------------------------------------------------------------------------------------------------------------------*/

workflow {
    filterGTF( ch_gtf )
    makeRef( filterGTF.out, ch_fa )
    cellrangerCount( ch_fastq.combine(makeRef.out) )
    // renameFeatures( cellrangerCount.out.sampleName.combine(filterGTF.out), cellrangerCount.out.countFiles )
    // runR( ch_rFile, renameFeatures.out.collect() )
}