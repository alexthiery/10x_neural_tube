#!/usr/bin/env nextflow

nextflow.preview.dsl=2

process renameFeatures {

    publishDir "${params.alignment_outDir}/cellrangerCounts_renamed",
    mode: "copy", overwrite: true

    input:
        tuple val(sample_name), path(gtf)
        path(cellrangerOut)

    output:
        path(cellrangerOut)
    
    """
    #!/usr/bin/env python
    # -*- coding: utf-8 -*-

    import re
    import gzip


    def matrix_edit(filt_gtf, feat):

        features = gzip.open(feat, 'rt')
        lab = ['MT', 'Z', 'W']

        output = {}
        with open(filt_gtf, 'rt') as f:
            for line in f.readlines():
                # split the words in each line of the gtf file and keep just the first word which is the chromosome ID
                chrm = line.split()[0]
                # if the chromosome ID is found in the labels list then continue
                if chrm in lab:
                    # within the line which this chromosome is on, gene_name is present then take the name and append it to output (dict key) with the chromosome ID (dict value)
                    if 'gene_name' in line:
                        gene_name = line.split('gene_name')[1].split(';')[0].replace('"', '').strip()
                        if gene_name not in output:
                            output[gene_name] = chrm
                    else:
                        gene_id = line.split('gene_id')[1].split(';')[0].replace('"', '').strip()
                        if gene_id not in output:
                            output[gene_id] = chrm

        # create new features file with gene names replaced if they are present in output
        new_feat = ''
        for line in features:
            gene = line.split()[1]
            if gene in output and not gene.startswith('MT-'):
                newline = re.sub("\t", "\t"+output[line.split()[1]]+"-", line, 1)
                new_feat += newline
            else:
                new_feat += line
                
        with gzip.open(feat, 'wt') as zipfile:
            zipfile.write(new_feat)

    matrix_edit(filt_gtf="${gtf}", feat="${cellrangerOut}/features.tsv.gz")
    """
 
}
