#!/usr/bin/env nextflow

nextflow.preview.dsl=2


process cellrangerCount {

    publishDir "${params.outDir}/cellrangerCounts",
        mode: "copy", overwrite: true

    label 'high_memory'

    input:
        tuple val(sample_id), val(sample_name), path('dir1/*'), path('dir2/*'), path(reference_genome)

    output:
        val sample_name, emit: sampleName
        path sample_name, emit: countFiles

    """
    #!/bin/bash
    
    cellranger count --id="cellrangerOut_${sample_name}" \
    --fastqs="dir1/${sample_id}","dir2/${sample_id}" \
    --sample=${sample_id} \
    --transcriptome=${reference_genome}

    mv cellrangerOut_${sample_name}/outs/filtered_feature_bc_matrix/*.gz ${sample_name}
    """
}