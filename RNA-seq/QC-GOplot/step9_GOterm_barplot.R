#"CP" = rgb(131,133,187, maxColorValue = 255),"GZ" = rgb(242,189,62, maxColorValue = 255)
setwd("~/lustrelyt/monkey-brain/merged-hic/RNA-Seq-new")
#### GO barplot
gob2a = read.table("results/DEgenes/chart_Rh_CPGZ_downinGZ_GOterm.txt", sep = "\t", header = T)
gob2a <- gob2a[which(gob2a$PValue < 0.05),]
gob2a <- gob2a[order(gob2a$PValue), ]
#cc = gob2a[which(gob2a$Category == "GOTERM_CC_DIRECT"),]
#bp = gob2a[which(gob2a$Category == "GOTERM_BP_DIRECT"),]
#mf = gob2a[which(gob2a$Category == "GOTERM_MF_DIRECT"),]

#gob2a = gob2a[which(gob2a$Category == "GOTERM_CC_DIRECT" | gob2a$Category == "GOTERM_BP_DIRECT" |
#                      gob2a$Category == "GOTERM_MF_DIRECT"),]
#gob2a$Term = substr(gob2a$Term, 12, 100)
#gob2a <- gob2a[c(1,2,3,5,6,13,16,17,20,29),]
gob2a <- gob2a[c(1:3,6,9,11,14,21,22,29),]

pdf("results/DEgenes/rh.CPGZ.downinGZ.GO.term.barplot.pdf", width = 26, height = 16)
par(mar = c(4,100,2,2))
barplot(-log10(gob2a$PValue), horiz = T, names.arg = gob2a$Term, las = 1, 
        col = rgb(242,189,62, maxColorValue = 255),
        cex.names = 3.0, cex.axis = 2.0, xlab = "-log10(P value)", cex.main = 4, axes = T, cex.lab = 1.5)
#axis(side = 1, at = seq(0,50,10), labels = seq(0,50,10), tick = T, lwd = 4, cex.axis = 2)
#segments(x0 = -log10(0.05), y0 = 0, x1 = -log10(0.05), y1 = 40, col = 'red', lwd = 4)
dev.off()
pdf("results/DEgenes/rh.CPGZ.downinGZ.GO.term.barplot.smallsize.pdf", width = 1.5 , height = 2.5)
par(mar = c(2,0,0,1))
barplot(-log10(gob2a$PValue), horiz = T, names.arg = rep("",10), las = 1, border = NA,
        col =rgb(131,133,187, maxColorValue = 255),
        cex.names = 1, cex.axis = 1.0, xlab = "", 
        cex.main = 1, axes = F, cex.lab = 1, xlim = c(0,26))
axis(side = 1, at = seq(0,25,5), labels = seq(0,25,5), tick = T, lwd = 0.5, cex.axis = 0.2)
#segments(x0 = -log10(0.05), y0 = 0, x1 = -log10(0.05), y1 = 40, col = 'red', lwd = 4)
dev.off()




pdf("results/GenomicFeatures/hgsBD.GO.term.barplot.pdf", width = 26, height = 12)
par(mar = c(4,100,2,2))
barplot(-log10(gob2a$PValue), horiz = T, names.arg = gob2a$Term, las = 1, 
        col = c( rep("chocolate3",5),rep("chartreuse4",4), rep("royalblue4",4)),
        cex.names = 3.0, cex.axis = 2.0, xlab = "-log10(P value)", cex.main = 4, axes = T, cex.lab = 1.5)
#axis(side = 1, at = seq(0,50,10), labels = seq(0,50,10), tick = T, lwd = 4, cex.axis = 2)
#segments(x0 = -log10(0.05), y0 = 0, x1 = -log10(0.05), y1 = 40, col = 'red', lwd = 4)
dev.off()
pdf("results/GenomicFeatures/hgsBD.GO.term.barplot.smallsize.pdf", width = 1 , height = 2)

par(mar = c(2,0,0,0))
barplot(-log10(gob2a$PValue), horiz = T, names.arg = rep("",13), las = 1, border = NA,
        col = ,
        cex.names = 1, cex.axis = 1.0, xlab = "", 
        cex.main = 1, axes = F, cex.lab = 1)
axis(side = 1, at = seq(0,2.5,0.5), labels = seq(0,2.5,0.5), tick = T, lwd = 0.5, cex.axis = 0.2)
#segments(x0 = -log10(0.05), y0 = 0, x1 = -log10(0.05), y1 = 40, col = 'red', lwd = 4)
dev.off()
