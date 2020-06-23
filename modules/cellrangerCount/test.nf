#!/usr/bin/env nextflow

// Define DSL2
nextflow.preview.dsl=2

// Log
log.info ("Starting test pipeline for cellrangerCount")

/*------------------------------------------------------------------------------------*/
/* Module inclusions 
--------------------------------------------------------------------------------------*/

include cellrangerCount from './cellrangerCount.nf'

/*------------------------------------------------------------------------------------*/
/* Params
--------------------------------------------------------------------------------------*/

testFastq = [['THI300A1', 'cellranger_count_hh4', "$baseDir/testData/THI300A1"]]
ref_genome = [["$baseDir/testData/galgal6_filtered_reference_genome"]]

/*------------------------------------------------------------------------------------*/
/* Define input channels
--------------------------------------------------------------------------------------*/
Channel
  .from(ref_genome)
  .set {ch_ref_genome}

Channel
  .from(testFastq)
  .set {ch_test_fastq}

workflow{
  //Run cellrangerCount
  cellrangerCount(ch_test_fastq, ch_ref_genome)
}

