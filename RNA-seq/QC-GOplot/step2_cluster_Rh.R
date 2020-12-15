library(RColorBrewer)
library(ggrepel)

all <- readRDS("data/monkey-summaried-expriment-rmd.rds")
sampData <- colData(all)
sampData <- cbind(sampData, batch = c(rep(1,12), rep(2,6)))
sampData$batch <- as.factor(sampData$batch)
colData(all) <- DataFrame(sampData)
#all <- orth.se
## construct dds object
dds <- DESeqDataSet(all, design = ~1)

## pre filtering
nrow(dds)
keep <- rowSums(counts(dds)) >= 20
dds <- dds[keep,]
nrow(dds)


## normalization
dds <- estimateSizeFactors(dds)
head(counts(dds),3)
vsd <- vst(dds, blind = T)
head(assay(vsd),3)
rld <- rlog(dds, blind = T)
head(assay(rld),3)


sampleDists.rld <- dist(t(assay(rld)))
sampleDistMatrix.rld <- as.matrix( sampleDists.rld )
rownames(sampleDistMatrix.rld) <- vsd$id
colnames(sampleDistMatrix.rld) <- NULL

sampleDists.vsd <- dist(t(assay(vsd)))
sampleDistMatrix.vsd <- as.matrix( sampleDists.vsd)
rownames(sampleDistMatrix.vsd) <- vsd$id
colnames(sampleDistMatrix.vsd) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

re.vsd <- pheatmap(sampleDistMatrix.vsd,
                   clustering_distance_rows = sampleDists.vsd,
                   clustering_distance_cols = sampleDists.vsd,
                   col = colors, main = "vsd")
re.rld <- pheatmap(sampleDistMatrix.rld,
                   clustering_distance_rows = sampleDists.rld,
                   clustering_distance_cols = sampleDists.rld,
                   col = colors, main = "rld")


cp <- cor(assay(vsd), method = "pearson")
pheatmap(cp)
hc <- hclust(d = dist(1-cp))
plot(hc)
cp <- cor(assay(rld), method = "pearson")
pheatmap(cp)
hc <- hclust(d = dist(1-cp))
plot(hc)

mat <- assay(vsd)
loci <- grep("VRK2", rownames(vsd))
boxplot(unlist(mat[loci,1:4]), unlist(mat[loci,5:22]), unlist(mat[loci,23:40]))

mat <- assay(rld)
cutvar <- function(mat, cutvalue){
  gvar <- apply(mat, 1, var)
  gcut <- quantile(gvar, cutvalue)
  loci <- which(gvar >= gcut)
  mcut <- mat[loci,]
  mcut <- as.matrix(mcut)
  colnames(mcut) = colnames(mat)
  return(mcut)
}

mat.filter <- cutvar(mat, 0.9)

res.pca <- prcomp(t(mat.filter), scale = TRUE)
eig <- (res.pca$sdev)^2
variance <- eig*100/sum(eig)
cumvar <- cumsum(variance)
eig.decathlon2.active <- data.frame(eig = eig, variance = variance,
                                    cumvariance = cumvar)
head(eig.decathlon2.active)

sampleData = data.frame(colData(all))
g2 <- ggplot(data.frame(res.pca$x[,1:2]), aes(x = PC1, y = PC2, color = sampleData$region, shape = sampleData$individual),
             label = sampleData$id) +
  geom_point(size =6) +
  geom_label_repel(aes(label = sampleData$id),
                   box.padding   = 0.35, 
                   point.padding = 0.5)+
  xlab(paste0("PC1 (", round(eig.decathlon2.active[1,2]), "% variance)")) +
  ylab(paste0("PC2 (", round(eig.decathlon2.active[2,2]), "% variance)")) +
  ggtitle(paste0("befcombat PCA (", nrow(mat.filter), "genes) "))+
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"))+
  theme(axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=18),  
        axis.title.x = element_text(size=28),
        axis.title.y = element_text(size=28),
        legend.text=element_text(size=28),
        plot.title = element_text(size=28, hjust = 0.5))
pdf("results/hg.rh.orthologs.vst.PCA.pdf", width = 12, height = 10)
g2
dev.off()



cp = cor(mat.filter, method = "pearson")
c = cor(mat.filter, method = "spearman")
#annote = as.data.frame(cbind("region" = c("CP","GZ","CP","GZ"), "individual" = c("07456A","07456A","07456B","07456B")))
annote = as.data.frame(colData(all)[,2:3])
annote$region = as.factor(annote$region)
annote$individual = as.factor(annote$individual)
#annote$replicate = as.factor(annote$replicate)
rownames(annote) <- colnames(cp)
pheatmap(as.matrix(cp),  cluster_rows = T, fontsize = 28, cluster_cols = T, annotation_col  = annote,
         colorRampPalette(c("white","red"))(100), show_colnames = F, clustering_distance_rows = "correlation", 
         clustering_distance_cols = "correlation")

pheatmap(as.matrix(cp),  cluster_rows = T, fontsize = 28, cluster_cols = T, annotation_col  = annote,
         colorRampPalette(c("blue","red"))(100), show_colnames = F, show_rownames = F, clustering_distance_rows = "correlation", 
         clustering_distance_cols = "correlation")




