####### correlation for PCA1 bwt replicates
## color set: 蓝色 rgb(100,149,237)  绿色 rgb(102,205,0)   红色 rgb(238,106,80)  紫色 rgb(131,133,187) 黄色 rgb(242,189,62)


########  correlation with gene density, gene expression, H3K27ac and H3K4me2

cp = fread("results/CP-500kb-ABcompartments-PC1-Raw.txt", header = T)
score = cp$score
score[which(is.na(score))] = 0
cor(score, cp$genedens, method = "spearman")

ac.den = fread("data/DNase.CP.peaks.density.500kb.bed")
ac.den = ac.den[which(ac.den$V1 != "chrY" & ac.den$V1 != "chrM"), ]
ac.den$chr = substring(ac.den$V1, 4,8)
ac.den$chr[which(ac.den$chr == "X")] = 21
ac.den$chr = as.numeric(as.character(ac.den$chr))
ac.den = ac.den[order(ac.den$chr),]
acd = ac.den$V5
cor(score, acd, method = "spearman")

#pdf("results/CP-PC1cor-genedensity-h3k27acdensity-withoutlable-new.pdf", width = 12, height = 6)
pdf("results/CP-PC1cor-genedensity-DNasedensity-withoutlable-new.pdf", width = 12, height = 6)
layout(matrix(c(1,2), nrow = 2, byrow = T))
x = seq(1,5659)
par(mar = c(4,8,4,8))
mp <- barplot(score, bty = "n", ylab = "PC1", cex.axis = 2, cex.lab = 2)
par(new=TRUE)
#barplot(cp$genedens, axes = FALSE, bty = "n", xlab = "", ylab = "", border = "red")
barplot(cp$genedens, axes = FALSE, bty = "n", xlab = "", ylab = "", border = rgb(100,149,237, maxColorValue = 255))
axis(side=4, at = pretty(range(cp$genedens)), cex.axis = 2, lwd.ticks = 2)
axis(side = 1, at = mp[seq(1,5001,1000)], labels =  c(rep("",6)) , cex.axis = 2 )
#mtext("Gene Density", side=4, line=3, col = "red", cex = 2)
mtext("Gene Density", side=4, line=3, col = rgb(100,149,237, maxColorValue = 255), cex = 2)


ac.den = fread("data/DNase.CP.peaks.density.500kb.bed")
ac.den = ac.den[which(ac.den$V1 != "chrY"), ]
ac.den$chr = substring(ac.den$V1, 4,8)
ac.den$chr[which(ac.den$chr == "X")] = 21
ac.den$chr = as.numeric(as.character(ac.den$chr))
ac.den = ac.den[order(ac.den$chr),]

acd = ac.den$V5
cor(score, acd, method = "spearman")

par(mar = c(4,8,4,8))
barplot(score, bty = "n",ylab = "PC1", cex.axis = 2, cex.lab = 2)
par(new=TRUE)
#mp <- barplot(acd, axes = FALSE, bty = "n", xlab = "", ylab = "", border = "chartreuse3")
mp <- barplot(acd, axes = FALSE, bty = "n", xlab = "", ylab = "", border = "chartreuse3")
axis(side=4, at = pretty(range(acd)), cex.axis = 2)
axis(side = 1, at = mp[seq(1,5001,1000)], labels =  c(rep("",6)) , cex.axis = 2)
#mtext("H3K27ac", side=4, line=3, col = "chartreuse3", cex = 2)
mtext("DNase", side=4, line=3, col = "chartreuse3", cex = 2)
dev.off()



