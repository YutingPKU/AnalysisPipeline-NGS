#################################
# plot pearson corelation heatmap for cp, gz, mef and liver
# plot A/B compartments
#################################
#### read input matrix and compartments
setwd("~/lustrelyt/monkey-brain/human-brain/Basic-Data/ABcompartment/")
cp <- "HiCPro-matrix/CP_500000.matrix"
gz <- "HiCPro-matrix/GZ_500000.matrix"
hc <- "HiCPro-matrix/HC_500000.matrix"
co <- "HiCPro-matrix/CO_500000.matrix"
hesc <- "HiCPro-matrix/hESC_500000.matrix"
lg <- "HiCPro-matrix/LG_500000.matrix"

xgi <- ygi <- "HiCPro-matrix/CP_500000_abs.bed"
cp.hic <- importC(cp, xgi, ygi)
gz.hic <- importC(gz, xgi, ygi)
hc.hic <- importC(hc, xgi, ygi)
co.hic <- importC(co, xgi, ygi)
hesc.hic <- importC(hesc, xgi, ygi)
lg.hic <- importC(lg, xgi, ygi)
cp.comp <- fread("compartments-500k-byrawmat/CP_500k_wholegenome_compartments.txt", header = T)
gz.comp <- fread("compartments-500k-byrawmat/GZ_500k_wholegenome_compartments.txt", header = T)
hesc.comp <- fread("compartments-500k-byrawmat/hESC_500k_wholegenome_compartments.txt", header = T)
lg.comp <- fread("compartments-500k-byrawmat/LG_500k_wholegenome_compartments.txt", header = T)
co.comp <- fread("compartments-500k-byrawmat/CO_500k_wholegenome_compartments_corrected.txt", header = T)
hc.comp <- fread("compartments-500k-byrawmat/HC_500k_wholegenome_compartments.txt", header = T)

#### plot function def
library(ggplot2)
library(gridExtra)
source("/lustre/user/liclab/liuyt/monkey-brain/rheMac-brain/Basic-Data/ABcompartment/plotAB.R")
plotpearson <- function(i, hic.intra, case_name){
  chrom <- paste0("chr",i,"chr",i)
  chr.hic <- hic.intra[[chrom]]
  chr.hic <- forceSymmetric(chr.hic)
  chr.p <- getPearsonMap(chr.hic)
  chr.cor <- intdata(chr.p)
  g <- pheatmap(chr.cor, cluster_rows = F, cluster_cols = F, show_rownames = F, show_colnames = F, 
                color = colorRampPalette(colors = c("navy", "red3"))(400), main = case_name, silent = T)
  return(g)
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
pdf("ABplot/observed.mat.public.pdf", height = 6, width = 28)
nf <- layout(matrix(c(1:12),3,4),heights = c(4, 0.5, 0.6))
#layout.show(nf)
for(i in c(1:22, "X")){
  #plotgrap(i, cp.hic, cp.comp, "CP")
  #plotgrap(i, gz.hic, gz.comp, "GZ")
  plotgrap(i, hesc.hic, hesc.comp, "hESC")
  plotgrap(i, lg.hic, lg.comp, "LG")
  plotgrap(i, co.hic, co.comp, "CO")
  plotgrap(i, hc.hic, hc.comp, "HC")
  
}
dev.off()

pdf("ABplot/pearson.cor.heatmap.public.pdf", height = 10, width = 40)
#nf <- layout(matrix(c(1:4),1,4))
for(i in c(1:22, "X")){
  plot_list=list()
  #g1 <- plotpearson(i, cp.hic, paste0("CP ","chr",i))
  #g2 <- plotpearson(i, gz.hic, paste0("GZ ","chr",i))
  g1 <- plotpearson(i, hesc.hic, paste0("hESC ", "chr", i))
  g2 <- plotpearson(i, lg.hic ,paste0("LG ","chr", i))
  g3 <- plotpearson(i, co.hic, paste0("CO chr", i))
  g4 <- plotpearson(i, hc.hic, paste0("HC chr", i))
  plot_list[[1]] = g1[[4]]
  plot_list[[2]] = g2[[4]]
  plot_list[[3]] = g3[[4]]
  plot_list[[4]] = g4[[4]]
  m <- marrangeGrob(plot_list, nrow=1, ncol=4, top = "")
  grid.draw(m)
}
dev.off()




