# HH4
norm.data <- readRDS("./results/R_results/RDS.files/norm.data.RDS")

temp <- "BACH1,CSRP2,DLX5,ELF1,ETS2,GATA2,GRHL3,HIC2,HOXA1,ID2,ID3,
KLF5,MEF2D,MSX1,MYC,NANOG,PDLIM5,PITX2,RASSF7,TFAP2A,VGLL1,
EPAS1,ATF3,CCND1,CEBPB,ETV1,ETV5,IRF7,KLF6,LIN28B,MGA,
NRIP1,OTX2,PRDM14,RUNX1T1,SNAI1,SOX4,TFAP2C,TGIF1,TRIM24,
PRDM1,CBFA2T2,CBX1,CBX2,CDR2,CITED2,CITED4,CRIP2,EHMT1,EOMES,ETV4,EZH2,
GATAD2B,GBX2,GLI2,HEY1,HIF1A,HIVEP3,HMGA1,ING5,IRX1,IRX2,KAT2B,
MAFA,MLXIP,MTA1,MYB,MYCN,NAB1,NFKB2,PDLIM4,PHTF1,
RNF14,SALL1,SETD2,STOX2,TBL1XR1,TCF7,TCF7L1,UBB,YEATS4,ZBTB1,ZHX2,ZIC3,ZNF277,
AATF,BMI1,CRTC1,DNMT3B,FOXA2,HES1,IRF8,IVNS1ABP,LHX5,PBX4,SAMD4A,
SIN3A,SMARCD1,SOX13,ZIC2,ZNF423,ZNF516,ZNF521,
JADE2,ENSGALG00000015705"

temp <- unlist(strsplit(gsub("\n", "", temp), ","))

temp[!temp %in% rownames(norm.data)]


HH4_missing <- c(
  'CDX2' = 'ENSGALG00000034983',
  'GRHL1' = 'ENSGALG00000038495',
  'DNMT3A' = 'ENSGALG00000041454',
  'ENC' = 'ENSGALG00000014930',
  'ZNF462' = 'ENSGALG00000015451',
  'BAZ1A' = 'ENSGALG00000010030',
  'SOX3' = 'ENSGALG00000054923',
  'SOX11' = 'ENSGALG00000016397',
  'CDCA7' = 'ENSGALG00000040833',
  'KDM4A' = 'ENSGALG00000010074',
  'MEIS2' = 'ENSGALG00000039118',
  'MYCL1' = 'ENSGALG00000024047',
  'RFX3' = 'ENSGALG00000010179',
  'ZNF469' = 'ENSGALG00000026534',
  'NSD1' = 'ENSGALG00000002971'
)
    

HH4_missing <- data.frame(gene_id = HH4_missing, gene_names = names(HH4_missing), row.names = NULL, stringsAsFactors = F)


HH4_genelist <- c(HH4_missing[HH4_missing$gene_id %in% norm.data@misc$geneIDs[,"gene_ID"],2], temp)

# write HH4 genelist
lapply(HH4_genelist, write, "./bin/network_genes/neural_induction_related_genes_HH4.txt", append = T)


# genes missing from merged data
HH4_missing[!HH4_missing$gene_id %in% merged.data@misc$geneIDs[,"gene_ID"],]

HH4_missing <- c(
"ENSGALG00000026534 ZNF469 - not in gtf",
"ENSGALG00000002971 NSD1 - pseudogene filtered in ref genome",
"ENSGALG00000034983 CDX2 - pseudogene filtered in ref genome",
"ENSGALG00000054923 SOX3 - pseudogene filtered in ref genome",
"ENSGALG00000016397 SOX11 - pseudogene filtered in ref genome")




# HH6
temp <- "ATF7IP,BRD8,FUBP1,HESX1,MBD3,MDFI,MID1,SOX2,SP4,TCF7L2,ZNF821,
HMGA2,LIN28A,LMO1,LMX1B,GMNN,PDCD4,SIX3,
CBY1,CNOT2,E2F3,FEZF2,MAML2,MEIS1,NKX6-2,SIX1,ZEB2"

temp <- unlist(strsplit(gsub("\n", "", temp), ","))

temp[!temp %in% rownames(norm.data)]


HH6_missing <- c(
  'E2F2' = 'ENSGALG00000051375',
  'FOXI3' = 'ENSGALG00000037457',
  'PTTG1' = 'ENSGALG00000001506',
  'HOXB1' = 'ENSGALG00000013454',
  'SMAD2' = 'ENSGALG00000036001'
)

HH6_missing <- data.frame(gene_id = HH6_missing, gene_names = names(HH6_missing), row.names = NULL, stringsAsFactors = F)

# names of genes in original dataset
HH6_genelist <- c(HH6_missing[HH6_missing$gene_id %in% norm.data@misc$geneIDs[,"gene_ID"],2], temp)

lapply(HH6_genelist, write, "./bin/network_genes/neural_induction_related_genes_HH6.txt", append = T)

# genes missing from merged dataset
HH6_missing[!HH6_missing$gene_id %in% norm.data@misc$geneIDs[,"gene_ID"],]

HH6_missing <- c("ENSGALG00000013454 HOXB1 - pseudogene")



