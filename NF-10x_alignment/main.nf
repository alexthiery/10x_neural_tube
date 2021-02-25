#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*-----------------------------------------------------------------------------------------------------------------------------
Pipeline params
-------------------------------------------------------------------------------------------------------------------------------*/
params.wGenes = "$baseDir/../input_files/wGenes.csv"
params.pymodifyGTF = "$baseDir/../bin/python/modifyGTF.py"

/*-----------------------------------------------------------------------------------------------------------------------------
Include modules
-------------------------------------------------------------------------------------------------------------------------------*/

include {projectHeader} from "$baseDir/../modules/projectHeader/projectHeader.nf"
include {modifyGTF} from "$baseDir/../modules/modifyGTF/modifyGTF.nf"
include {filterGTF} from "$baseDir/../modules/filterGTF/filterGTF.nf"
include {makeRef} from "$baseDir/../modules/makeRef/makeRef.nf"
include {cellrangerCount} from "$baseDir/../modules/cellrangerCount/cellrangerCount.nf"

/*-----------------------------------------------------------------------------------------------------------------------------
Init
-------------------------------------------------------------------------------------------------------------------------------*/
// Show banner
log.info projectHeader()

// Header log info
def summary = [:]
summary['Run Name']               = workflow.runName
summary['Sample File']            = params.samplesheet
summary['Fasta File']             = params.fasta
summary['GTF File']               = params.gtf
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
    .fromPath( params.gtf )
    .set {ch_gtf}

Channel
    .fromPath( params.fasta )
    .set {ch_fasta}

Channel
    .fromPath( params.samplesheet )
    .splitCsv(header: ['sample_id', 'sample_name', 'dir1', 'dir2'], skip: 1 )
    .map { row -> [row.sample_id, row.sample_name, file(row.dir1), file(row.dir2)] }
    .set { ch_fastq }

/*-----------------------------------------------------------------------------------------------------------------------------
Main workflow
-------------------------------------------------------------------------------------------------------------------------------*/

workflow {
    modifyGTF( file(params.pymodifyGTF), ch_gtf, file(params.wGenes) )
    filterGTF( modifyGTF.out )
    makeRef( filterGTF.out, ch_fasta )
    cellrangerCount( ch_fastq.combine(makeRef.out) )
}