#################################
# plot pearson corelation heatmap for cp, gz, mef and liver
# plot A/B compartments
## color set: 蓝色 rgb(100,149,237)  绿色 "chartreuse3"   红色 rgb(238,106,80)  紫色 rgb(131,133,187) 黄色 rgb(242,189,62)

#################################
library(HiTC)
#### read input matrix and compartments
setwd("~/lustrelyt/monkey-brain/merged-hic/basic-data/ABcompartment")
cp <- "data/CP_500000.matrix"
gz <- "data/GZ_500000.matrix"
xgi <- ygi <- "data/CP_500000_abs.bed"
cp.hic <- importC(cp, xgi, ygi)
gz.hic <- importC(gz, xgi, ygi)

i <- 1

cp.Comp <- fread("results/CP-500kb-ABcompartments-PC1-Raw.txt", header = T)
gz.Comp <- fread("results/GZ-500kb-ABcompartments-PC1-Raw.txt", header = T)
cp.comp <- cp.Comp[which(cp.Comp$seqnames == paste0("chr",i)),]$score
gz.comp <- gz.Comp[which(gz.Comp$seqnames == paste0("chr",i)),]$score

cp.gend <- cp.Comp[which(cp.Comp$seqnames == paste0("chr",i)),]$genedens
gz.gend <- gz.Comp[which(gz.Comp$seqnames == paste0("chr",i)),]$genedens

h3k27ac <- fread("data/H3K27ac.peaks.density.500kb.bed", header = F)
h3k27ac <- h3k27ac[which(h3k27ac$V1 == paste0("chr",i)),]$V4

h3k4me2 <- fread("data/H3K4me2.peaks.density.500kb.bed", header = F)
h3k4me2 <- h3k4me2[which(h3k4me2$V1 == paste0("chr",i)),]$V4

x <- seq(1, length(cp.comp))
#### plot function def
library(ggplot2)
library(gridExtra)
library(ComplexHeatmap)
library(fields)
source("/lustre/user/liclab/liuyt/monkey-brain/rheMac-brain/Basic-Data/ABcompartment/plotAB.R")
plotpearson <- function(i, hic.intra, case_name){
  chrom <- paste0("chr",i,"chr",i)
  chr.hic <- hic.intra[[chrom]]
  chr.hic <- forceSymmetric(chr.hic)
  chr.p <- getPearsonMap(chr.hic)
  chr.cor <- intdata(chr.p)
  par(mar = c(0,8,0,2))
  g <- pheatmap(chr.cor, cluster_rows = F, cluster_cols = F, show_rownames = F, show_colnames = F, 
           color = colorRampPalette(colors = c("navy", "red3"))(400), main = case_name, legend = F, silent = T)
  

}
plotgrap <- function(i, hic.intra, compartments, case_name){
  chrom <- paste0("chr",i,"chr",i)
  chr.hic <- hic.intra[[chrom]]
  chr.hic <- forceSymmetric(chr.hic)
  mat <- intdata(chr.hic)
  comp <- compartments[which(compartments$seqnames == paste0("chr", i)), ]
  casename <- paste0(case_name, " chr",i)
  cp.g <- plotAB(mat, casename, comp)
}
##


mat <- as.matrix(chr.cor)
#pdf("results/CP-chr1-ABcompartment-with-markerSingal.pdf", height = 8, width = 8)
layout(matrix(c(1:4),4,1,byrow = TRUE),heights = c(6,0.5,0.5,0.5), widths = rep(6, times = 4))
#layout.show()
par(mar = c(1,10,0,4))
#image(x=0:(nrow(mat)-1),y=0:(ncol(mat)-1), mat, col = colorRampPalette(colorRamps::blue2red(10))(100),  
#      axes = F, xlab = "", ylab = "")
image(x=0:(nrow(mat)-1),y=0:(ncol(mat)-1), mat, col = colorRampPalette(c("cadetblue1","green","yellow"))(100),  
      axes = F, xlab = "", ylab = "")
#image.plot(mat, legend.only = T, legend.shrink = 0.2, col = colorRampPalette(colorRamps::blue2red(10))(100), add = F)
#plotpearson(i, cp.hic, paste0("CP chr",i))

#add.image(0,0, mat, col = colorRampPalette(c("navy","red3"))(400), adj.x = 0.5, adj.y = 0.5, image.height = 2, image.width = 2 )

par(mar = c(0,8,0,2))
y <- gz.comp
plot(x, y, type = 'n',xaxt = "n",frame.plot = F,font.lab =4,ylab = "PCA1vec",xlab = "")
segments(x,rep(0,10),x,y,ljoin = 0,lwd = 1.2,col =ifelse(y>0,rgb(255,12,0, maxColorValue = 255), rgb(0,37,255, maxColorValue = 255)))

par(mar = c(0,8,0,2))
y <- gz.gend
plot(x,y,type = 'n',xaxt = "n",frame.plot = F,font.lab =4,ylab = "GeneDensity",xlab = "")
segments(x,rep(0,10),x,y,ljoin = 0,lwd = 1.2,col =rgb(100,149,237, maxColorValue = 255))

par(mar = c(0,8,0,2))
y <- h3k27ac
plot(x,y,type = 'n',xaxt = "n",frame.plot = F,font.lab =4,ylab = "H3K27ac",xlab = "")
segments(x,rep(0,10),x,y,ljoin = 0,lwd = 1.2,col ="chartreuse3")


#dev.off()

