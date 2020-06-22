//params.gtf = "~/dev/genomes/galgal6/Gallus_gallus.GRCg6a.97.gtf "
//params.fa = "galgal6/Gallus_gallus.GRCg6a.dna.toplevel.fa"

params.outgenomename = "galgal6_filtered_ref_genome"
params.outdir = "./temp_out"

gtf_ch = Channel.fromPath(params.gtf)

process filterGTF {

    input:
        path gtf from gtf_ch
    output:
        file 'filtered_genome.gtf' into filtered_genome

    """
    #!/bin/bash

    # this step filters out genes based on the gene biotypes listed in attributes.
    cellranger mkgtf ${gtf} filtered_genome.gtf \
    --attribute=gene_biotype:protein_coding \
    --attribute=gene_biotype:lncRNA \
    """
}

fa_ch = Channel.fromPath(params.fa)

process makeRef {
    cpus = params.threads

    publishDir "${params.outdir}",
        mode: "copy", overwrite: true

    input:
        file filt_genome from filtered_genome
        path fa from fa_ch

    output:
        path("galgal6_filtered_reference_genome")

    """
    #!/bin/bash

    # make reference
    cellranger mkref \
    --genome=galgal6_filtered_reference_genome \
    --fasta= ${fa} \
    --genes=${filt_genome} \
    --nthreads=${params.threads} \
    --memgb=${params.mem}
    """
}