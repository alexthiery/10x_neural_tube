#!/usr/bin/env nextflow

nextflow.preview.dsl=2

/*-----------------------------------------------------------------------------------------------------------------------------
Include modules
-------------------------------------------------------------------------------------------------------------------------------*/

include projectHeader from "$baseDir/modules/projectHeader/projectHeader.nf"
include filterGTF from "$baseDir/modules/filterGTF/filterGTF.nf"
include makeRef from "$baseDir/modules/makeRef/makeRef.nf"
include cellrangerCount from "$baseDir/modules/cellrangerCount/cellrangerCount.nf"
// include renameFeatures from "$baseDir/renameFeatures.nf"

/*-----------------------------------------------------------------------------------------------------------------------------
Pipeline params
-------------------------------------------------------------------------------------------------------------------------------*/
params.outgenomename = "galgal6_filtered_ref_genome"
params.outdir = "./temp_out"

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
//     .fromPath( params.sampleDir )
//     .set { ch_features }

// Channel
//     .fromPath(params.extraData)
//     .set { ch_extraData }

// params.rFile = "$baseDir/bin/R/seurat_full.R"
// params.customFunctions = "$baseDir/bin/R/my_functions"

/*-----------------------------------------------------------------------------------------------------------------------------
Main workflow
-------------------------------------------------------------------------------------------------------------------------------*/

workflow {
    filterGTF( ch_gtf )
    makeRef( filterGTF.out, ch_fa )
    cellrangerCount( ch_fastq.combine(makeRef.out) )
    // renameFeatures( cellrangerCount.out, filterGTF.out )
    // runR( ch_extraData, renameFeatures.out )
}