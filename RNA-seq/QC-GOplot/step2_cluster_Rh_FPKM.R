library(RColorBrewer)
library(ggrepel)
cutvar <- function(mat, cutvalue){
  gvar <- apply(mat, 1, var)
  gcut <- quantile(gvar, cutvalue)
  loci <- which(gvar >= gcut)
  mcut <- mat[loci,]
  mcut <- as.matrix(mcut)
  colnames(mcut) = colnames(mat)
  return(mcut)
}


all <- readRDS("data/monkey-exonbygene-summaried-expriment-rmd.rds")
sampleData <- colData(all)

mat <- read.xlsx("results/monkey-RNASeq-cqnnormalized-combat-removebatch-FPKM-log2-rmd.xlsx", colNames = T, rowNames = T)


mat.filter <- cutvar(mat, 0.9)

res.pca <- prcomp(t(mat.filter), scale = TRUE)
eig <- (res.pca$sdev)^2
variance <- eig*100/sum(eig)
cumvar <- cumsum(variance)
eig.decathlon2.active <- data.frame(eig = eig, variance = variance,
                                    cumvariance = cumvar)
head(eig.decathlon2.active)



g2 <- ggplot(data.frame(res.pca$x[,1:2]), aes(x = PC1, y = PC2, color = sampleData$region, shape = sampleData$individual),
             label = sampleData$id) +
  geom_point(size =6) +
  geom_label_repel(aes(label = sampleData$id),
                   box.padding   = 0.35, 
                   point.padding = 0.5)+
  xlab(paste0("PC1 (", round(eig.decathlon2.active[1,2]), "% variance)")) +
  ylab(paste0("PC2 (", round(eig.decathlon2.active[2,2]), "% variance)")) +
  ggtitle(paste0("befcombat PCA (", nrow(mat.filter), "genes) "))+
  #scale_color_manual(values=c("CP" = "cornflowerblue","GZ" = "chartreuse3"))+
  scale_color_manual(values=c("CP" = rgb(131,133,187, maxColorValue = 255),"GZ" = rgb(242,189,62, maxColorValue = 255)))+
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
g1 <- ggplot(data.frame(res.pca$x[,1:2]), aes(x = PC1, y = PC2, color = sampleData$region, shape = sampleData$individual),
             label = sampleData$id) +
  geom_point(size =2) +
  #geom_label_repel(aes(label = sampleData$id),
  #                 box.padding   = 0.35, 
  #                 point.padding = 0.5)+
  xlab(paste0("PC1 (", round(eig.decathlon2.active[1,2]), "% variance)")) +
  ylab(paste0("PC2 (", round(eig.decathlon2.active[2,2]), "% variance)")) +
  ggtitle(paste0("cqnFPKM combat PCA (", nrow(mat.filter), "genes) "))+
  scale_color_manual(values=c("CP" = rgb(131,133,187, maxColorValue = 255),"GZ" = rgb(242,189,62, maxColorValue = 255)))+
  #scale_color_manual(values=c("CP" = "cornflowerblue","GZ" = "chartreuse3"))+
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", size = linesize))+
  theme(legend.position = "none",
    axis.text.x = element_text(size=fontsize),
        axis.text.y = element_text(size=fontsize),  
        axis.title.x = element_text(size=fontsize),
        axis.title.y = element_text(size=fontsize),
        legend.text=element_text(size=fontsize),
        plot.title = element_text(size=fontsize, hjust = 0.5))


pdf("results/clusterPlot/monkey.allsamples.cqnFPKM.top5vargenes.PCA.newcolor.pdf", width = 8, height = 6)
g2
dev.off()
pdf("results/clusterPlot/monkey.allsamples.cqnFPKM.top5vargenes.PCA.newcolor.smallsize.pdf", width = 2, height = 2)
g1
dev.off()

cp = cor(mat.filter, method = "pearson")
c = cor(mat.filter, method = "spearman")
#annote = as.data.frame(cbind("region" = c("CP","GZ","CP","GZ"), "individual" = c("07456A","07456A","07456B","07456B")))
annote = as.data.frame(colData(all)[,2:3])
annote$region = as.factor(annote$region)
annote$individual = as.factor(annote$individual)
#annote$replicate = as.factor(annote$replicate)
rownames(annote) <- colnames(cp)

Var1        <- c(rgb(131,133,187, maxColorValue = 255), rgb(242,189,62, maxColorValue = 255))
names(Var1) <- c("CP", "GZ")
Var2 <- c("cornflowerblue","chartreuse3","coral3")
names(Var2) <- c("07456A","07456B","11002B")
anno_colors <- list(region = Var1, individual = Var2)
breakslist <- seq(0.7,1,0.002)

fontsize = 6
linesize = 0.4
pdf("results/clusterPlot/monkey.allsamples.cqnFPKM.top5vargenes.pearson.correlation.heatmap.pdf", width = 4, height = 4)
pheatmap(as.matrix(cp),  cluster_rows = T,  cluster_cols = T, annotation_col  = annote, annotation_colors = anno_colors,
         colorRampPalette(c("white","red"))(100), show_colnames = F, clustering_distance_rows = "correlation", 
         clustering_distance_cols = "correlation", border_color = NA, fontsize= fontsize)
pheatmap(as.matrix(cp),  cluster_rows = T, fontsize = fontsize, cluster_cols = T, annotation_col  = annote, annotation_colors = anno_colors,
         colorRampPalette(c("white","red"))(100), show_colnames = F, show_rownames = F,clustering_distance_rows = "correlation", 
         clustering_distance_cols = "correlation", border_color = NA)
pheatmap(as.matrix(cp),  cluster_rows = T, fontsize = fontsize, cluster_cols = T, annotation_col  = annote, annotation_colors = anno_colors,
         colorRampPalette(c("white","red"))(100), show_colnames = F, show_rownames = F,clustering_distance_rows = "correlation", 
         clustering_distance_cols = "correlation", border_color = NA, annotation_legend = F, annotation_names_col = F)
pheatmap(as.matrix(cp),  cluster_rows = T, fontsize = fontsize, cluster_cols = T, annotation_col  = annote, annotation_colors = anno_colors,
         colorRampPalette(colorRamps::blue2red(10))(152), breaks = breakslist, show_colnames = F,show_rownames = F, clustering_distance_rows = "correlation", 
         clustering_distance_cols = "correlation", border_color = NA)
pheatmap(as.matrix(cp),  cluster_rows = T, fontsize = fontsize, cluster_cols = T, annotation_col  = annote, annotation_colors = anno_colors,
         colorRampPalette(colorRamps::blue2red(10))(152), breaks = breakslist, show_colnames = F,show_rownames = F, clustering_distance_rows = "correlation", 
         clustering_distance_cols = "correlation", border_color = NA, annotation_legend = F, annotation_names_col = F)
dev.off()
pdf("results/clusterPlot/monkey.allsamples.cqnFPKM.top5vargenes.pearson.correlation.heatmap.smallsize.pdf", width = 3, height = 3)
pheatmap(as.matrix(cp),  cluster_rows = T, fontsize = fontsize, cluster_cols = T, annotation_col  = annote, annotation_colors = anno_colors,
         colorRampPalette(c("white","red"))(100), show_colnames = F, show_rownames = F,clustering_distance_rows = "correlation", 
         clustering_distance_cols = "correlation", border_color = NA, annotation_legend = F, annotation_names_col = F)
pheatmap(as.matrix(cp),  cluster_rows = T, fontsize = fontsize, cluster_cols = T, annotation_col  = annote, annotation_colors = anno_colors,
         colorRampPalette(colorRamps::blue2red(10))(152), breaks = breakslist, show_colnames = F,show_rownames = F, clustering_distance_rows = "correlation", 
         clustering_distance_cols = "correlation", border_color = NA, annotation_legend = F, annotation_names_col = F)
dev.off()
pheatmap(as.matrix(cp),  cluster_rows = T, fontsize = 28, cluster_cols = T, annotation_col  = annote,
         colorRampPalette(c("white","red"))(100), show_colnames = F, show_rownames = F, clustering_distance_rows = "correlation", 
         clustering_distance_cols = "correlation")



