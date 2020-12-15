
all <- readRDS("data/monkey-exonbygene-summaried-expriment-rmd.rds")

#rownames(all) <- unlist(orth.hg.rh)
#species <- c(rep("hg",18), rep("rh",18))
#sampleData <- colData(all)
#sampleData <- cbind(sampleData, data.frame(species))
#colData(all) <- sampleData
## construct dds object
dds <- DESeqDataSet(all, design = ~ region)

## pre filtering
nrow(dds)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
nrow(dds)

## vsd
vst <- vst(dds, blind = T)
## DE
dds <- DESeq(dds)

## results 
res <- results(dds, contrast = c("region", "GZ","CP"))

res.up <- res[which(res$log2FoldChange > 1 & res$padj < 0.01),]
res.down <- res[which(res$log2FoldChange < -1 & res$padj < 0.01), ]


res.up <- data.frame(res.up)
write.table(rownames(res.up), "results/DEgenes/rh.CPGZ.upinGZ.fld4.fdr001.gid.txt", row.names = F, col.names = F, quote = F, sep = "\t")
res.down <- data.frame(res.down)
write.table(row.names(res.down), "results/DEgenes/rh.CPGZ.downinGZ.fld4.fdr001.gid.txt", row.names = F, col.names = F, quote = F, sep = "\t")


plotMA(res)
res.dat <- data.frame(res)
res.up <- res.dat[which(res.dat$log2FoldChange>=1 & res.dat$padj <= 0.01 ),]
res.up <- res.up[order( res.up$log2FoldChange, -log10(res.up$pvalue), decreasing = T),]
res.down <- res.dat[which(res.dat$log2FoldChange <= -1 & res.dat$padj <= 0.01),]
res.down <- res.down[order( abs(res.down$log2FoldChange), -log10(res.down$pvalue), decreasing = T),]


####### volcano plot
gene_list = res.dat[complete.cases(res.dat),]
#gene_list = gene_list[order(abs(gene_list$log2FoldChange), -log10(gene_list$pvalue), decreasing = T),]
gene_list$threshold = as.factor(abs(gene_list$log2FoldChange) > 2 & gene_list$padj <= 0.01)
gene_list$gene = rownames(gene_list)
loci.up = match(rownames(res.up), gene_list$gene)[1:20]
loci.down = match(rownames(res.down), gene_list$gene)[1:20]
g = ggplot(data=gene_list, aes(x=log2FoldChange, y=-log10(pvalue), colour=threshold)) +
  geom_point(size=1) +
  theme(legend.position = "none") +
  xlim(c(-8, 8)) + ylim(c(0, 30)) +
  xlab("log2 fold change") + ylab("-log10 p-value")+
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"))+
  theme(axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=18),  
        axis.title.x = element_text(size=28),
        axis.title.y = element_text(size=28),
        legend.position = "none",
        plot.title = element_text(size=28, hjust = 0.5))+
  scale_colour_manual(values = c("FALSE"= "black", "TRUE"="red"))+
  geom_text_repel(data=gene_list[c(loci.up, loci.down),], aes(x=log2FoldChange, y=-log10(pvalue),label=gene))

pdf("results/P0-DE-gene-volcano-top20.pdf", width = 12, height = 10)
g
dev.off()
