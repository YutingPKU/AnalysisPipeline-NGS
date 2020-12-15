
####### construct DESeq objects
setwd("~/lustrelyt/monkey-brain/merged-hic/RNA-Seq-new")
hg.count <- readRDS("data/human.feature.count.byXSAnno.genes.rds")
rh.count <- readRDS("data/macaque.feature.count.byXSAnno.genes.rds")

#h.cp.count <- hg.count$counts[,1:9]
#m.cp.count <- rh.count$counts[,c(1:3,7:9,13:15)]
allcount <- cbind(hg.count$counts, rh.count$counts)
colnames(allcount) <- c(paste0("hg", seq(1,18)), paste0("rh", seq(1,18)))

coldata <- cbind(id = colnames(allcount), species = rep(c("hg","rh"), each = 18), 
                 regions = c(rep("CP",9), rep("GZ",9), rep(c("CP","GZ"), each = 3), rep(c("CP","GZ"), each = 3),rep(c("CP","GZ"), each = 3)))
coldata <- DataFrame(coldata)
rownames(coldata) <- colnames(allcount)

dds <- DESeqDataSetFromMatrix(countData = allcount, colData = coldata, design = ~species)


####### vsd and rld normalize
## pre filtering
nrow(dds)
keep <- rowSums(counts(dds)) >= 20
dds <- dds[keep,]
nrow(dds)


## normalization
head(counts(dds),3)
vsd <- vst(dds, blind = T)
head(assay(vsd),3)
rld <- rlog(dds, blind = T)
head(assay(rld),3)
dds <- estimateSizeFactors(dds)

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
cp <- cor(assay(vsd), method = "spearman")
pheatmap(cp)


mat <- assay(vsd)
loci <- grep("VRK2", rownames(vsd))
boxplot(unlist(mat[loci,1:4]), unlist(mat[loci,5:22]), unlist(mat[loci,23:40]))

mat <- assay(vsd)
cutvar <- function(mat, cutvalue){
  gvar <- apply(mat, 1, var)
  gcut <- quantile(gvar, cutvalue)
  loci <- which(gvar >= gcut)
  mcut <- mat[loci,]
  mcut <- as.matrix(mcut)
  colnames(mcut) = colnames(mat)
  return(mcut)
}

mat.filter <- cutvar(mat, 0.8)

res.pca <- prcomp(t(mat.filter), scale = TRUE)
eig <- (res.pca$sdev)^2
variance <- eig*100/sum(eig)
cumvar <- cumsum(variance)
eig.decathlon2.active <- data.frame(eig = eig, variance = variance,
                                    cumvariance = cumvar)
head(eig.decathlon2.active)

sampleData = coldata
g2 <- ggplot(data.frame(res.pca$x[,1:2]), aes(x = PC1, y = PC2, color = sampleData$species, shape = sampleData$regions),
             label = sampleData$id) +
  geom_point(size =6) +
  geom_label_repel(aes(label = sampleData$id),
                   box.padding   = 0.35, 
                   point.padding = 0.5)+
  xlab(paste0("PC1 (", round(eig.decathlon2.active[1,2]), "% variance)")) +
  ylab(paste0("PC2 (", round(eig.decathlon2.active[2,2]), "% variance)")) +
  ggtitle(paste0("vsd PCA (", nrow(mat.filter), "genes) "))+
  #scale_color_manual(values=c("CP" = "cornflowerblue","GZ" = "chartreuse3"))+
  scale_color_manual(values=c("hg" = "cornflowerblue","rh"="chartreuse3"))+
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
g1 <- ggplot(data.frame(res.pca$x[,1:2]), aes(x = PC1, y = PC2, color = sampleData$species, shape = sampleData$regions),
             label = sampleData$id) +
  geom_point(size =2) +
  #geom_label_repel(aes(label = sampleData$id),
  #                 box.padding   = 0.35, 
  #                 point.padding = 0.5)+
  xlab(paste0("PC1 (", round(eig.decathlon2.active[1,2]), "% variance)")) +
  ylab(paste0("PC2 (", round(eig.decathlon2.active[2,2]), "% variance)")) +
  ggtitle(paste0("vsd PCA (", nrow(mat.filter), "genes) "))+
  scale_color_manual(values=c("hg" = "cornflowerblue","rh"="chartreuse3"))+
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


pdf("results/clusterPlot/monkey.human.vsd.top20vargenes.PCA.newcolor.pdf", width = 8, height = 6)
g2
dev.off()
pdf("results/clusterPlot/monkey.human.vsd.top20vargenes.PCA.newcolor.smallsize.pdf", width = 2, height = 2)
g1
dev.off()








#################### plot hcluster 
# Customized colors
hc <- hclust(dist(1-cor(mat, method = "pearson")))
hcd <- as.dendrogram(hc)
#hc.re <- rev(hcd)
hc.re <- as.dendrogram(hcd)

pdf("results/clusterPlot/monkey.human.top100genes.vsd.pearson.correlation.cluster.pdf", height = 3, width =6)
par(mar = c(2,1,1,4))
hc.re %>%  set("branches_lwd", 0.4) %>% 
  set("labels_colors", k=2, values = c("chartreuse3","cornflowerblue")) %>% set("labels_cex", c(0.8,0.8)) %>% set("leaves_cex", 1) %>%
  set("leaves_pch", 19) %>% set("leaves_col", c(rep(rgb(242,189,62, maxColorValue = 255),9), rep(rgb(131,133,187, maxColorValue = 255),9),
                                                rep(rgb(131,133,187, maxColorValue = 255),9),rep(rgb(242,189,62, maxColorValue = 255),9))) %>% 
  plot( horiz = F)
dev.off()



cp = cor(mat.filter, method = "pearson")
c = cor(mat.filter, method = "spearman")
#annote = as.data.frame(cbind("region" = c("CP","GZ","CP","GZ"), "individual" = c("07456A","07456A","07456B","07456B")))
annote = as.data.frame(coldata[,2:3])
annote$regions = as.factor(annote$regions)
annote$species = as.factor(annote$species)
#annote$replicate = as.factor(annote$replicate)
rownames(annote) <- colnames(cp)
Var1        <- c(rgb(131,133,187, maxColorValue = 255), rgb(242,189,62, maxColorValue = 255))
names(Var1) <- c("CP", "GZ")
Var2 <- c("cornflowerblue","chartreuse3")
names(Var2) <- c("hg","rh")
anno_colors <- list(regions = Var1, species = Var2)
breakslist <- seq(0.4,1,0.002)

fontsize = 6
linesize = 0.4
pdf("results/clusterPlot/monkey.human.top20genes.vsd.pearson.correlation.heatmap.pdf", width = 4, height = 4)
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
         colorRampPalette(colorRamps::blue2red(10))(300), breaks = breakslist, show_colnames = F,show_rownames = F, clustering_distance_rows = "correlation", 
         clustering_distance_cols = "correlation", border_color = NA)
pheatmap(as.matrix(cp),  cluster_rows = T, fontsize = fontsize, cluster_cols = T, annotation_col  = annote, annotation_colors = anno_colors,
         colorRampPalette(colorRamps::blue2red(10))(300), breaks = breakslist, show_colnames = F,show_rownames = F, clustering_distance_rows = "correlation", 
         clustering_distance_cols = "correlation", border_color = NA, annotation_legend = F, annotation_names_col = F)
dev.off()
pdf("results/clusterPlot/monkey.human.top20genes.vsd.pearson.correlation.heatmap.smallsize.pdf", width = 3, height = 3)
pheatmap(as.matrix(cp),  cluster_rows = T, fontsize = fontsize, cluster_cols = T, annotation_col  = annote, annotation_colors = anno_colors,
         colorRampPalette(c("white","red"))(300), show_colnames = F, show_rownames = F,clustering_distance_rows = "correlation", 
         clustering_distance_cols = "correlation", border_color = NA, annotation_legend = F, annotation_names_col = F)
pheatmap(as.matrix(cp),  cluster_rows = T, fontsize = fontsize, cluster_cols = T, annotation_col  = annote, annotation_colors = anno_colors,
         colorRampPalette(colorRamps::blue2red(10))(300), breaks = breakslist, show_colnames = F,show_rownames = F, clustering_distance_rows = "correlation", 
         clustering_distance_cols = "correlation", border_color = NA, annotation_legend = F, annotation_names_col = F)
dev.off()