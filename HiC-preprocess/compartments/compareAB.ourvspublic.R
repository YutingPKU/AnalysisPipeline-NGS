##########################################
# compare our results with Ren's results
##########################################

setwd("/lustre/user/liclab/liuyt/SV-3Dgenome/brain-hicdata/Compartment_primary_cohort/")
h1 <- fread("h1.pc.bw.bedGraph")
lg <- fread("LG.pc.bw.bedGraph")
co <- fread("CO.pc.bw.bedGraph")
hc <- fread("HC.pc.bw.bedGraph")

setwd("~/lustrelyt/monkey-brain/human-brain/Basic-Data/")
h1.o <- fread("ABcompartment/compartments-500k-byrawmat/hESC_1M_wholegenome_compartments.txt", header = T)
lg.o <- fread("ABcompartment/compartments-500k-byrawmat/LG_1M_wholegenome_compartments.txt", header = T)
co.o <- fread("ABcompartment/compartments-500k-byrawmat/CO_1M_wholegenome_compartments.txt", header = T)
hc.o <- fread("ABcompartment/compartments-500k-byrawmat/HC_1M_wholegenome_compartments.txt", header = T)

plotcomp <- function(comp, i, name){
  names(comp) <- c("chr","start","end","value")
  chrom <- paste0("chr", i)
  comp <- comp[which(comp$chr == chrom), ]
  comp$value <- comp$value/max(comp$value)
  #comp$score[is.na(comp$score)] = 0
  len <- comp$end[length(comp$end)]/1000000
  x <- seq(1, len)
  y <- rep(0, times = len)
  
  par(mar = c(0,6,2,2))
  plot(x,y,type = 'n',xaxt = "n",frame.plot = F,font.lab =4, ylab = name, xlab = "")
  #segments(x,rep(0,10),x,y,ljoin = 0,lwd = 1.2,col =ifelse(y>0,"red", "blue")) # plot A/B
  rect(xleft = comp$start/1000000, ybottom = rep(0, times = length(comp$start)), xright = comp$end/1000000, ytop = comp$value, 
       col =ifelse(comp$value>0,rgb(219,110,49, max = 255), 
                   rgb(72,111,167, max = 255)), border = NA)
}

plotcomp.o <- function(comp, i){
  chrom <- paste0("chr", i)
  comp <- comp[which(comp$seqnames == chrom), ]
  comp$score[is.na(comp$score)] = 0
  y = as.numeric(comp$score)
  x<- seq(1:(length(y)))
  d <- as.numeric(comp$genedens)
  par(mar = c(0,6,2,2))
  plot(x,y,type = 'n',xaxt = "n",frame.plot = F,font.lab =4, ylab = "PCA1vec",xlab = "")
  #segments(x,rep(0,10),x,y,ljoin = 0,lwd = 1.2,col =ifelse(y>0,"red", "blue")) # plot A/B
  rect(xleft = comp$start/1000000, ybottom = rep(0, length(y)), xright = comp$end/1000000, ytop = comp$score, 
       col =ifelse(y>0,rgb(219,110,49, max = 255), rgb(72,111,167, max = 255)), border = NA)
  
  par(mar = c(1,6,0,2))
  plot(1:length(y),d,type = 'n',xaxt = "n",frame.plot = F,font.lab =4,ylab = "Gene density",xlab = "")
  segments(x,rep(0,10),x,d,ljoin = 0,lwd = 1.2,col = 'black')
 
}


pdf("ABcompartment/ABplot/compare.ourvspublic.ABresults.1M.pdf", height = 20, width = 26)
nf <- layout(matrix(c(1:12),12,1,byrow = TRUE),heights = rep(4, times = 12))
layout.show(nf)
for(i in c(1:22, "X")){
  plotcomp(h1, i, "H1")
  plotcomp.o(h1.o, i)
  plotcomp(lg, i, "LG")
  plotcomp.o(lg.o, i)
  plotcomp(co, i, "CO")
  plotcomp.o(co.o, i)
  plotcomp(hc, i, "HC")
  plotcomp.o(hc.o, i)
}
dev.off()