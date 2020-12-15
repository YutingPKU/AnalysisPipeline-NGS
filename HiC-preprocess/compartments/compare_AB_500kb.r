## color set: 蓝色 rgb(100,149,237)  绿色 "chartreuse3"   红色 rgb(238,106,80)  紫色 rgb(131,133,187) 黄色 rgb(242,189,62)
library(refGenome)
library(GenomicRanges)

CP <- fread("results/CP-500kb-ABcompartments-PC1-Raw.txt", header = T, stringsAsFactors = F)
GZ <- fread("results/GZ-500kb-ABcompartments-PC1-Raw.txt", stringsAsFactors = F)
MEF <- fread("/lustre/user/liclab/liuyt/monkey-brain/rheMac-brain/Basic-Data/ABcompartment/ABcompartments-500k-rawmat/10084A-MEF-ABcompartment-500k-rawmat.txt", header = T, stringsAsFactors = F)

mat <- cbind("seqnames" = CP$seqnames, "start"=CP$start, "end"=CP$end, "CP-PCA1" = CP$score, "CP-AB" = as.character( CP$ccompartments), 
             "GZ-PCA1" = GZ$score, "GZ-AB" = as.character(GZ$ccompartments), "MEF-PCA1" = MEF$score, "MEF-AB" = as.character(MEF$ccompartments) )
mat <- data.frame(mat)
mat$CP.PCA1 = as.numeric(as.character(mat$CP.PCA1))
mat$GZ.PCA1 = as.numeric(as.character(mat$GZ.PCA1))
mat$MEF.PCA1= as.numeric(as.character(mat$MEF.PCA1))


loci <- which(mat$CP.AB == mat$GZ.AB & mat$CP.AB == mat$MEF.AB)
diffmat <- mat[-loci,c(4,6,8)]
diffmat <- diffmat[complete.cases(diffmat), ]
diffmat$CP.PCA1 = as.numeric(as.character(diffmat$CP.PCA1))
diffmat$GZ.PCA1 = as.numeric(as.character(diffmat$GZ.PCA1))
diffmat$MEF.PCA1 = as.numeric(as.character(diffmat$MEF.PCA1))
diffmat <- diffmat[which(diffmat$CP.PCA1 <= 0.1 & diffmat$GZ.PCA1 <= 0.1 & diffmat$MEF.PCA1 <= 0.1),]
par(mar = c(2,2,2,2))
pdf("results/ABcompare/diffAB_heatmap_CPGZMEF_new_withcolorbar.pdf", width = 6, height = 9)
par(mar = c(2,2,2,2))
pheatmap(diffmat, cluster_cols = T, cluster_rows = T, show_rownames = F, show_colnames = T, labels_col = c("CP","GZ","MEF"),
       color = colorRampPalette(grDevices::topo.colors(6))(100))
pheatmap(diffmat, cluster_cols = T, cluster_rows = T, show_rownames = F, show_colnames = T, labels_col = c("CP","GZ","MEF"),
                 color = colorRampPalette(colorRamps::blue2red(10))(100))

dev.off()



tmat <- mat[complete.cases(mat),c(4,6,8)]
tmat <- apply(tmat,2,function(xx){
  xx <- as.numeric(as.character(xx))
})
tmat <- data.frame(tmat)
tmat <- tmat[which(abs(tmat$CP.PCA1) <= 0.1 & abs(tmat$GZ.PCA1) <= 0.1 & abs(tmat$MEF.PCA1) <= 0.1),]
pdf("plot/AB-compare-wholegenome.pdf", width = 4, height = 6)
pheatmap(tmat, cluster_rows = T, cluster_cols = T, show_rownames = F)
dev.off()
heatmap(as.matrix(tmat))



length(which(mat$CP.AB == "A" & mat$GZ.AB == "A"))
length(which(mat$CP.AB == "B" & mat$GZ.AB == "B"))
length(which(mat$CP.AB == "A" & mat$GZ.AB == "B"))
length(which(mat$CP.AB == "B" & mat$GZ.AB == "A"))
length(which(mat$CP.AB != mat$GZ.AB))/5683
length(which(mat$MEF.AB != mat$CP.AB))/5683
length(which(mat$MEF.AB != mat$GZ.AB))/5683

##################### bar plot of AB change in three cells
count <- data.frame("Region" = c("CP","CP","GZ","GZ","MEF","MEF"), "Compartments" = rep(c("A","B"),3), "Number" = c(2670,2900,2672,2898,2619,2950))
count <- data.frame("Region" = c(rep("CP",4), rep("GZ",4)),
                    "Compartments" = c("stable A","CP-A","CP-B","stable B","stable A","GZ-B","GZ-A","stable B"), 
                    "Number" = rep(c(2095,203,217,2810),2))

p <- ggplot(data=count, aes(y=Number, x=Region, fill=factor(Compartments))) + 
  geom_bar(stat="identity", colour="black") +
  theme_bw()+
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background = element_blank(),panel.border = element_blank())+
  scale_fill_manual(values = c(
    "stable A"=rgb(100,149,237, maxColorValue = 255),
    "CP-A" = rgb(238,106,80, maxColorValue = 255),
    "CP-B" = "cadetblue1",
    "stable B" = "chartreuse3",
    "GZ-A" = "cadetblue1",
    "GZ-B" = rgb(238,106,80, maxColorValue = 255)))
pdf("results/ABcompare/CPGZ_switchAB_barplot.pdf", height = 6, width = 6)
p
dev.off()





################### GO of genes in A/B switch regions
lociA2B = mat[which(mat$CP.AB == "A" & mat$GZ.AB =="B" ),1:3]
lociB2A = mat[which(mat$CP.AB == "B" & mat$GZ.AB == "A"),1:3]
lociStable = mat[which(mat$CP.AB == mat$GZ.AB),1:3]

A2Bbed = makeGRangesFromDataFrame(lociA2B)
B2Abed = makeGRangesFromDataFrame(lociB2A)
Stablebed = makeGRangesFromDataFrame(lociStable)

gtffile = "/lustre/user/liclab/liuyt/monkey-brain/ref_genome/ensemble/Mmul_8.0.1/Macaca_mulatta.Mmul_8.0.1.89.chr.gtf"
setwd("/lustre/user/liclab/liuyt/monkey-brain/ref_genome/ensemble/Mmul_8.0.1/")
ens <- ensemblGenome()
read.gtf(ens, "Macaca_mulatta.Mmul_8.0.1.89.chr.gtf")
my_gene <- getGenePositions(ens)
my_gr <- with(my_gene, GRanges(paste0("chr", seqid), IRanges(start, end), strand, id = gene_id, name = gene_name))

A2Bgene = findOverlaps(my_gr, A2Bbed, type = "within")
A2Bgene = mcols(my_gr)[queryHits(A2Bgene),1]
A2Bgene = A2Bgene[complete.cases(A2Bgene)]

B2Agene = findOverlaps(my_gr, B2Abed, type = "within")
B2Agene = mcols(my_gr)[queryHits(B2Agene),1]
B2Agene = B2Agene[complete.cases(B2Agene)]

Stablegene = findOverlaps(my_gr, Stablebed, type = "within")
Stablegene = mcols(my_gr)[queryHits(Stablegene),1]
Stablegene = Stablegene[complete.cases(Stablegene)]


write.table(A2Bgene, "results/table/CPA-to-GZB-compartments-500kb-genesname.txt", sep = "\t", row.names = F, col.names = F, quote = F)
write.table(B2Agene, "results/table/CPB-to-GZA-compartments-500kb-genesname.txt", sep = "\t", row.names = F, col.names = F, quote = F)



################  gene expression change for regions which switch compartments
#expmat = counts(dds, normalized = T)
#expmat = expmat
exp <- read.xlsx("../../RNA-Seq/data/monkey-RNASeq-cqnnormalized-FPKM-log2-rmd-exonbygene.xlsx", rowNames = T, colNames = T)
exp <- 2^exp
#se <- readRDS("/lustre/user/liclab/liuyt/monkey-brain/rna-seq/results/monkey-exonbygene-summaried-expriment-rmd.rds")
#dds <- DESeqDataSet(se, design = ~region)
#vst <- vst(dds, blind = T)
#exp <- assay(vst)
colnames(exp) <- as.character(colData(se)[,1])
expmat = exp[,c(grep("CP", colnames(exp)), grep("GZ", colnames(exp)))]

#expmat$CP = rowMeans(expmat[,1:3])
#expmat$GZ = rowMeans(expmat[,4:6])
ratio = rowMeans(expmat[,1:9])/rowMeans(expmat[,10:18])


sgene = cbind("compartment" = rep("Stable", nrow(expmat)), "ratio" = ratio)
sgene = data.frame(sgene)
rownames(sgene) <- rownames(expmat)
sgene <- sgene[match(c(A2Bgene, B2Agene, Stablegene), rownames(sgene)),]

sgene$compartment = as.character(sgene$compartment)
sgene$compartment[match(A2Bgene, rownames(sgene))] = "A2B"
sgene$compartment[match(B2Agene, rownames(sgene))] = "B2A"
sgene$compartment = as.factor(sgene$compartment)
sgene$ratio = as.numeric(as.character(sgene$ratio))

#sgene$logr = log2(sgene$ratio)

library(ggplot2)
ylim1 = boxplot.stats(sgene$ratio)$stats[c(1, 5)]
p1 = p0 + coord_cartesian(ylim = ylim1*1.05)
g = ggplot(sgene, aes(x = compartment, y = ratio, fill = compartment))+
    geom_boxplot(notch = F, outlier.shape = NA, show.legend = F)+
  stat_boxplot(geom ='errorbar')+
   #stat_summary(fun.y=mean, geom="point", shape=23, size=4)+
  coord_cartesian(ylim = c(0,2))+
  xlab("Compartments") +
  ylab("log2 FPKM (CP/GZ)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_blank())+
  theme(axis.text.x = element_text(size=24),
        axis.text.y = element_text(size=24),  
        axis.title.x = element_text(size=28),
        axis.title.y = element_text(size=28),
        axis.line.y = element_line(colour = "black", size = 1.5),
        axis.ticks.x=element_blank(),
        axis.ticks.y = element_line(size = 1.5),
        legend.text=element_text(size=18),
        plot.title = element_text(size=22, hjust = 0.5))+
  #theme(legend.position="none")+
  scale_x_discrete(labels=c("A to B","B to A","Stable")) +
  scale_y_continuous(breaks = seq(-0.8,0.8,0.4), labels = seq(-0.8,0.8,0.4))

pdf("plot/ABswitch-CPGZ-geneexp-boxplot.pdf", width = 6, height = 6)
g
dev.off()



########################  correlation btw PCA1 value and histone marker signal enrichment
test = cbind("loci" = CP$start, "PC1" = CP$score, "Genedesity" = CP$genedens)
test = data.frame(test)
test = test[complete.cases(test),]
par(mar=c(2, 4, 4, 6) + 0.1)
barplot(test$PC1, col = rep("firebrick2", length(test$PC1)), border = NA, space=0)
par(new=TRUE)
barplot(test$Genedesity,space = 0,axes=FALSE, ylim = c(0,60))
