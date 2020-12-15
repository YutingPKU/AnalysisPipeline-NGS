
setwd("~/lustrelyt/monkey-brain/merged-hic/RNA-Seq-new/XSAnno")


anno.hg <- fread("step2_Blat/blatFilter/blatFiltered.hg19TorheMac8.hg19.withname.bed", header = F, sep = "\t")
anno.rh <- fread("step2_Blat/blatFilter/blatFiltered.hg19TorheMac8.rheMac8.withname.bed", header = F, sep = "\t")

hgcount <- readRDS("step3_simNGS/hg19.simluated.RNASeq.featureCount.byexon.allMO.rds")
rhcount <- readRDS("step3_simNGS/rheMac8.simluated.RNASeq.featureCount.byexons.allMO.rds")
allcount <- cbind(hgcount$counts, rhcount$counts)
colnames(allcount) <- c(paste0("hg", seq(1,10)), paste0("rh", seq(1,10)))

coldata <- cbind(id = colnames(allcount), species = rep(c("hg","rh"), each = 10))
coldata <- DataFrame(coldata)
rownames(coldata) <- colnames(allcount)

dds <- DESeqDataSetFromMatrix(countData = allcount, colData = coldata, design = ~species)
dds <- DESeq(dds)
res <- results(dds, contrast = c("species", "hg","rh"))

# simFilter from blat
loci.exon <- which(res$padj < 0.05)
loci.gene <- which(res$padj < 0.05)
gene.filter <- rownames(res)[loci.gene]
exon.filter <- rownames(res)[loci.exon]




anno.hg.filter <- anno.hg[-loci.exon,]
anno.rh.filter <- anno.rh[-loci.exon, ]

loci.gene <- which(anno.hg.filter$V5 %in% gene.filter)
anno.hg.filter <- anno.hg.filter[-loci.gene, ]
anno.rh.filter <- anno.rh.filter[-loci.gene, ]

write.table(anno.hg.filter, "step3_simNGS/simFiltered.byDIM.hg19TorheMac8.hg19.bed", row.names = F, col.names = F, quote = F, sep = "\t")
write.table(anno.rh.filter, "step3_simNGS/simFiltered.byDIM.hg19TorheMac8.rheMac8.bed", row.names = F, col.names = F, quote = F, sep = "\t")


