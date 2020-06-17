import re


input_gtf = open("/Users/alex/dev/genomes/galgal6/Gallus_gallus.GRCg6a.97.gtf", "rt")

gtf = input_gtf.read()


test_files = ['/Users/alex/dev/data/10x_neural/neural_induction_related_genes_HH4.txt', '/Users/alex/dev/data/10x_neural/neural_induction_related_genes_HH6.txt']

for i in test_files:

	print('Searching file:', i, "\nGenes missing from galgal6 GTF:")

	with open(i, "rt") as gene_file:

		genes = gene_file.read().split('\n')

		for gene in genes:
			if bool(re.search(gene, gtf)) == False:
				print(gene)

	print("\n\n")