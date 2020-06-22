#!/usr/bin/env nextflow

nextflow.preview.dsl=2

process makeRef {
    cpus = params.cpus

    publishDir "${params.outdir}",
        mode: "copy", overwrite: true

    input:
        tuple val(sample_id), val(sample_name), path(fastq1), path(fastq2)
        path(reference_genome)

    output:
        path("sample_name")

"""
cellranger count --id=${sample_name} \
--fastqs=${fastq1},${fastq2} \
--sample=${sample_id} \
--transcriptome=${reference_genome}
"""



sample_id = THI300A1
sample_name = cellranger_count_hh4
fastq1 = ~/Documents/10x/Crick_output_raw/SC19040/190524_K00102_0343_AH7MYHBBXY/fastq/THI300A1
fastq2 = ~/Documents/10x/Crick_output_raw/SC19040/190614_K00102_0353_BH7T3WBBXY/fastq/THI300A1