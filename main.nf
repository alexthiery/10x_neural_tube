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
temp
-------------------------------------------------------------------------------------------------------------------------------*/


def names1 = []
new File(params.folder1).eachDir() { file->
    names1 << "$baseDir/" + params.folder1.split('/')[-1] + "/" + file.getName()   
}

// def names2 = []
// new File(params.folder2).eachDir() { file->
//     names2 << "$baseDir/" + params.folder2.split('/')[-1] + "/" + file.getName()  
// }

names = names1

def index = [
    ['THI300A1', 'cellranger_count_hh4'],
    ['THI300A3', 'cellranger_count_hh4'],
    ['THI300A4', 'cellranger_count_hh4'],
    ['THI300A6', 'cellranger_count_hh4']
    ]

def testFastq = []
for (item in index) {
    testFastq << item + names.findAll{it -> it.contains(item[0])}
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

