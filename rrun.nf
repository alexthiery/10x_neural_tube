#!/usr/bin/env nextflow

nextflow.preview.dsl=2

/*-----------------------------------------------------------------------------------------------------------------------------
Pipeline params
-------------------------------------------------------------------------------------------------------------------------------*/
params.rFile = "$baseDir/bin/R/seurat_full.R"
params.customFuncs = "$baseDir/bin/R/custom_functions"
params.networkGenes = "$baseDir/bin/network_genes"

/*-----------------------------------------------------------------------------------------------------------------------------
Include modules
-------------------------------------------------------------------------------------------------------------------------------*/

include projectHeader from "$baseDir/modules/projectHeader/projectHeader.nf"
include runR from "$baseDir/modules/runR/runR.nf"

/*-----------------------------------------------------------------------------------------------------------------------------
Init
-------------------------------------------------------------------------------------------------------------------------------*/
// Show banner
log.info projectHeader()

// Header log info
def summary = [:]
summary['Run Name']               = custom_runName ?: workflow.runName
summary['Metadata File']          = params.metadata
summary['Fasta File']             = params.fa
summary['GTF File']               = params.gtf
summary['Max Resources']          = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
if (workflow.containerEngine)     summary['Container'] = "$workflow.containerEngine - $workflow.container"
summary['Output Dir']             = params.outdir
summary['Launch Dir']             = workflow.launchDir
summary['Working Dir']            = workflow.workDir
summary['Script Dir']             = workflow.projectDir
summary['User']                   = workflow.userName
summary['Config Profile']         = workflow.profile
log.info summary.collect { k,v -> "${k.padRight(20)}: $v" }.join('\n')
log.info "-\033[2m--------------------------------------------------\033[0m-"

/*------------------------------------------------------------------------------------*/
/* Define input channels
--------------------------------------------------------------------------------------*/

Channel
    .fromPath(params.counts)
    .collect()
    .set { ch_feat }

Channel
    .fromPath(params.rFile)
    .set { ch_rFile }

/*-----------------------------------------------------------------------------------------------------------------------------
Main workflow
-------------------------------------------------------------------------------------------------------------------------------*/

workflow {
    runR( ch_rFile, ch_feat )
}