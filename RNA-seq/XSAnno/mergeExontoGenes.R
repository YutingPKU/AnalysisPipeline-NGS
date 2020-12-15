
hg <- fread("XSAnno.byself.Ensemble.hg19TorheMac8.hg19.withname.bed", header = F, sep = "\t")
names(hg) <- c("chr","start","end","eid","gid","tid","gn","tn","type")
hg.unique <- hg[!duplicated(hg$gid),]
table(hg.unique$type)

pdf("../results/XSAnno.byself.Ensemble.hg19.category.piechart.pdf", width = 1.5, height = 1.5)
linesize = 0.5
par(mai= c(0.2,0.2,0.2,0.2))
pie(c(16983,3220,2636,2360,1523), labels =c("protein coding16983","lincRNA3220","antisense2636","pseudogene2360","other1523"), main = "", 
    col = rev(brewer.pal(5,"Blues")), lwd = linesize, cex = 0.5)
par(mai= c(0.1,0.1,0.1,0.1))
pie(c(16983,3220,2636,2360,1523), labels =NA, main = "", 
    col = rev(brewer.pal(5,"Blues")), lwd = linesize)
dev.off()