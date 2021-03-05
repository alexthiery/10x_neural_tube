#!/usr/bin/env nextflow

nextflow.enable.dsl=2


process cellrangerCount {

    publishDir "${params.output}/cellrangerCounts",
        mode: "copy", overwrite: true

    label 'process_high'

    input:
        tuple val(sample_id), val(sample_name), path(fastqs) path(reference_genome)

    output:
        val sample_name, emit: sampleName
        path sample_name, emit: countFiles

    """
    #!/bin/bash
    
    cellranger count --id="cellrangerOut_${sample_name}" \
    --fastqs=${fastqs} \
    --sample=${sample_id} \
    --transcriptome=${reference_genome}

    mkdir ${sample_name}

    mv cellrangerOut_${sample_name}/outs/filtered_feature_bc_matrix/*.gz ${sample_name}
    """
}