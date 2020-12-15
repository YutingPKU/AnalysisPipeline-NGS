
load("data/orthologs.hg.rh.rds")
all <- orth.se
dds <- DESeqDataSet(all, design = ~ species.X + region)

## pre filtering
nrow(dds)
keep <- rowSums(counts(dds)) >= 20
dds <- dds[keep,]
nrow(dds)

## vsd
vst <- vst(dds, blind = T)
mat <- assay(vst)

## hgs gene
hgs <- fread("~/lustrelyt/monkey-brain/merged-hic/basic-data/TAD/Version4/data/GenomiFeatures/genes/TSS.upstream2kb.GENCODE-V19.inHgsBD.extend20kb.geneid.txt", header = F)
#names(hgs) <- c("geneid","genename","chr1","x1","x2","chr2","y1","y2")
#hgs$dist <- abs(hgs$y1-hgs$x1)

gname <- hgs$genename
loci <- match(gname, rownames(mat))
loci <- loci[which(!is.na(loci))]
fc <- lapply(loci, function(x){
  log2(mean(mat[x, 19:36])/mean(mat[x,1:18]))
})
fc <- unlist(fc)

hgs$log2FC <- fc

write.xlsx(hgs, file = "data/genename.in.human_CPGZ_common_loop_hgs.withloop_withFC.xlsx", colName = T)


hg.mean <- rowMeans(mat[,19:36])
rh.mean <- rowMeans(mat[,1:18])
boxplot(hg.mean[loci], rh.mean[loci])





hgs <- fread("../LMD-RNA-Seq/data/TSSinloops/human_rhesus_commloop_TSS2kb.inH3K27ac.genename.txt", header = F)
names(hgs) <- c("geneid","genename","chr1","x1","x2","chr2","y1","y2")
hgs$dist <- abs(hgs$y1-hgs$x1)

gname <- hgs$V1
loci <- match(gname, rownames(mat))
loci <- loci[which(!is.na(loci))]
case.mat <- data.frame(species = rep(colData(all)[,5], length(loci)), value = unlist(split(mat[loci,], seq(nrow(mat[loci,])))))

g12 <- ggplot(case.mat, aes(species, value, fill=species)) + geom_boxplot(outlier.shape = NA)+
  ylab("vsd normalized expression") +
  xlab("")+
  ggtitle(gname)+
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background = element_blank(),panel.border = element_blank())+
  theme(axis.text.x = element_text(size=24),axis.text.y = element_text(size=18),  
        axis.title.x = element_text(size=28),axis.title.y = element_text(size=18),
        axis.line.y = element_line(colour = "black", size = 1.5),axis.ticks.x=element_blank(),
        axis.line.x = element_line(colour = "black", size = 1.5),
        axis.ticks.y = element_line(size = 1.5),legend.text=element_text(size=18),
        axis.ticks.length = unit(.25, "cm"),
        plot.title = element_text(size=22, hjust = 0.5))
g12
