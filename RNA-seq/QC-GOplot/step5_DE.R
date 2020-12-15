
load("data/orthologs.hg.rh.rds")
all <- orth.se
#rownames(all) <- unlist(orth.hg.rh)
#species <- c(rep("hg",18), rep("rh",18))
#sampleData <- colData(all)
#sampleData <- cbind(sampleData, data.frame(species))
#colData(all) <- sampleData
## construct dds object
dds <- DESeqDataSet(all, design = ~ region + species.X )

## pre filtering
nrow(dds)
keep <- rowSums(counts(dds)) >= 20
dds <- dds[keep,]
nrow(dds)

## vsd
vst <- vst(dds, blind = T)
## DE
dds <- DESeq(dds)

## results 
res <- results(dds, contrast = c("species.X","hg","rh"))

res.up <- res[which(res$log2FoldChange > 2 & res$padj < 0.01),]
res.down <- res[which(res$log2FoldChange < -2 & res$padj < 0.01), ]

de <- rbind(res.up, res.down)
de <- data.frame(de)
de$id <- rownames(de)
write.xlsx(de, "results/hg.rh.orthologs.DEG.2logFC.xlsx", rowNames = T, colNames = T)

res.up <- data.frame(res.up)
write.table(rownames(res.up), "results/gene.up.in.human.gname.txt", row.names = F, col.names = F, quote = F, sep = "\t")
res.down <- data.frame(res.down)
write.table(row.names(res.down), "results/gene.down.in.human.gname.txt", row.names = F, col.names = F, quote = F, sep = "\t")


fc <- apply(counts(dds),1, FUN = function(vec){
  log2(unlist(mean(vec[19:36]))/mean(unlist(vec[1:18])))
} )
