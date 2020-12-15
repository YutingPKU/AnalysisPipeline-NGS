################  normalize RNA-Seq raw data: GC content, gene length, sequencing depth
########## read count table and gene models

gclen <- fread("code/ensemble_Mmu8_gene_GC_lengths.txt", header = F)
names(gclen) <- c("id", "len", "gc")

all = readRDS("data/monkey-exonbygene-summaried-expriment-rmd.rds")
mat <- assay(all)
names(mat) <- unlist(as.character(colData(all)[,1]))
#rownames(mat) <- unlist(as.character(rowData(all)[,1]))

cqn.subset = cqn(mat, lengths = gclen$len, x = gclen$gc, lengthMethod = "fixed")
RPKM.cqn <- cqn.subset$y + cqn.subset$offset
colnames(RPKM.cqn) <- unlist(as.character(colData(all)[,1]))

#####filter by 
boxplot(RPKM.cqn)
#write.table(RPKM.cqn, "results/monkey-RNASeq-cqnnormalized-FPKM-log2-rmd-exonbygene.txt", sep = '\t', row.names = T, col.names = T, quote = F)
write.xlsx(RPKM.cqn, "results/monkey-RNASeq-cqnnormalized-FPKM-log2-rmd-exonbygene.xlsx", colNames = T, rowNames = T)

RPM <- sweep(log2(mat + 1), 2, log2(colSums(mat)/10^6))
RPKM.std <- sweep(RPM, 1, log2(gclen$len/ 10^3))
write.xlsx(RPKM.std, "results/monkey-RNASeq-standward-FPKM-log2-rmd-exonbygene.xlsx", colNames = T, rowNames = T)


