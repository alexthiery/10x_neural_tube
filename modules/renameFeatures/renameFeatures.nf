#!/usr/bin/env nextflow

nextflow.preview.dsl=2

process renameFeatures {

    publishDir "${params.outdir}",
    mode: "copy", overwrite: true

    input:
        tuple path(cellrangerOut), path(gtf)

    output:
        path(cellrangerOut)
    
   """
    #!/usr/bin/env python

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
                        # if there is no gene name - i.e. if there is no annotation, then take the gene ID instead so that we can still identify ensembl IDs within the downstream analysis
                        else:
                            gene_id = line.split('gene_id')[1].split(';')[0].replace('"', '').strip()
                            if gene_id not in output:
                                output[gene_id] = chrm

        # make new character vector which will be become the new ammended features file
        # split the features lines and take the second element [2] as this contains the gene names
        # if gene name in features file is in the output list we generated above, then keep append the line with the value from the dictionary, but only if MT- is not already there
        
        new_feat = ''
        for line in features:
            gene = line.split()[1]
            if gene in output and not gene.startswith('MT-'):
                newline = re.sub("\t", "\t"+output[line.split()[1]]+"-", line, 1)
                new_feat += newline
            else:
                new_feat += line

        # overwrite features file with edited MT- label features
        with gzip.open(feat, 'wt') as zipfile:
            zipfile.write(new_feat)

    matrix_edit(filt_gtf="${gtf}", feat="${cellrangerOut}/outs/filtered_feature_bc_matrix/features.tsv.gz")
    """
 
}



