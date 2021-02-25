#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*-----------------------------------------------------------------------------------------------------------------------------
Pipeline params
-------------------------------------------------------------------------------------------------------------------------------*/
params.rFile = "$baseDir/../bin/R/seurat_full_slim.R"
params.customFuncs = "$baseDir/../bin/R/custom_functions"
params.networkGenes = "$baseDir/../input_files/network_expression.csv"


/*-----------------------------------------------------------------------------------------------------------------------------
Include modules
-------------------------------------------------------------------------------------------------------------------------------*/

include {runR} from "$baseDir/../modules/runR/runR.nf"
include {projectHeader} from "$baseDir/../modules/projectHeader/projectHeader.nf"


// Show banner
log.info projectHeader()

// Header log info
def summary = [:]
summary['Run Name']               = workflow.runName
summary['Max Resources']          = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
summary['Output Dir']             = workflow.launchDir/params.output
summary['Launch Dir']             = workflow.launchDir
summary['Working Dir']            = workflow.workDir
summary['Script Dir']             = workflow.projectDir
summary['User']                   = workflow.userName
summary['Config Profile']         = workflow.profile
log.info summary.collect { k,v -> "${k.padRight(20)}: $v" }.join('\n')
log.info "-\033[2m--------------------------------------------------\033[0m-"

/*------------------------------------------------------------------------------------
Define input channels
--------------------------------------------------------------------------------------*/
Channel
    .fromPath( params.input )
    .collect()
    .set { ch_cellrangercounts }

/*-----------------------------------------------------------------------------------------------------------------------------
Main workflow
-------------------------------------------------------------------------------------------------------------------------------*/
workflow {
    runR( file(params.rFile), ch_cellrangercounts )
}
