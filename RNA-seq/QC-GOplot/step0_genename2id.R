## transfer enemble transcript id to gene symbol
library('biomaRt')
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
#mart <- useDataset("mmulatta_gene_ensembl", useMart("ensembl"))

ts <- fread("results/DEgenes/rh.CPGZ.downinGZ.fld4.fdr001.gid.txt", header = F)
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("hgnc_symbol"),values=ts,mart= mart)

write.table(unique(G_list$hgnc_symbol), "results/DEgenes/rh.CPGZ.downinGZ.fld4.fdr001.gname.txt", row.names = F, col.names = F, quote = F, sep = "\t")


