####### construct DESeq objects
setwd("~/lustrelyt/monkey-brain/merged-hic/RNA-Seq-new")
hg.count <- readRDS("data/human.feature.count.byXSAnno.genes.rds")
rh.count <- readRDS("data/macaque.feature.count.byXSAnno.genes.rds")

h.cp.count <- hg.count$counts[,1:9]
m.cp.count <- rh.count$counts[,c(1:3,7:9,13:15)]
#h.cp.count <- hg.count$counts[,10:18]
#m.cp.count <- rh.count$counts[,-c(1:3,7:9,13:15)]
allcount <- cbind(h.cp.count, m.cp.count)
colnames(allcount) <- c(paste0("hg", seq(1,9)), paste0("rh", seq(1,9)))

coldata <- cbind(id = colnames(allcount), species = rep(c("hg","rh"), each = 9))
coldata <- DataFrame(coldata)
rownames(coldata) <- colnames(allcount)

dds <- DESeqDataSetFromMatrix(countData = allcount, colData = coldata, design = ~species)

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
res <- results(dds, contrast = c("species","hg","rh"))
res.data <- data.frame(res)

res.up <- res.data[which(res$log2FoldChange > 2 & res$padj < 0.01),]
res.down <- res.data[which(res$log2FoldChange < -2 & res$padj < 0.01), ]


res.up <- data.frame(res.up)
write.table(rownames(res.up), "results/DEgenes/intra-species/up.DEgene.in.human.log2FC2.FDR001.txt", row.names = F, col.names = F, quote = F, sep = "\t")
res.down <- data.frame(res.down)
write.table(row.names(res.down), "results/DEgenes/intra-species/down.DEgene.in.human.log2FC2.FDR001.txt", row.names = F, col.names = F, quote = F, sep = "\t")
write.table(res.data, "results/DEgenes/intra-species/GZ.DEgene.human.vs.macaque.log2FC2.FDR001.txt", row.names = T,
            col.names = T, quote = F, sep = "\t")

plotMA(res)
res.dat <- data.frame(res)
res.up <- res.dat[which(res.dat$log2FoldChange>=2 & res.dat$padj <= 0.01 ),]
res.up <- res.up[order( res.up$log2FoldChange, -log10(res.up$pvalue), decreasing = T),]
res.down <- res.dat[which(res.dat$log2FoldChange <= -2 & res.dat$padj <= 0.01),]
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
  xlim(c(-10, 10)) + ylim(c(0, 300)) +
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
