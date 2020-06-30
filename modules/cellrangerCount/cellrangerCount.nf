#!/usr/bin/env nextflow

nextflow.preview.dsl=2


process cellrangerCount {
    cpus = params.threads

    publishDir "${params.outdir}",
        mode: "copy", overwrite: true

    input:
        tuple val(sample_id), val(sample_name), path('dir1/*'), path('dir2/*')
        path(reference_genome)

    output:
        path("sample_name")

    """
    hostname
    ulimit
    cellranger count --id=${sample_name} \
    --fastqs="dir1/${sample_id}","dir2/${sample_id}" \
    --sample=${sample_id} \
    --transcriptome=${reference_genome}
    """
}