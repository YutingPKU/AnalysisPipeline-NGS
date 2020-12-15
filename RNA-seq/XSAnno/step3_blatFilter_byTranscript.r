#############
## chimp ##
#############
source("/lustre/user/liclab/liuyt/monkey-brain/merged-hic/RNA-Seq-new/XSAnno/XSAnno/bin/Functions_BlatFilter_byTranscript.r")
setwd("~/lustrelyt/monkey-brain/merged-hic/RNA-Seq-new/XSAnno/step2_Blat/blat/")
sp12sp1="blat.exon.hg19Tohg19.filtered.txt"
sp22sp2="blat.exon.rheMac8TorheMac8.filtered.txt"
sp12sp2="blat.exon.hg19TorheMac8.filtered.txt"
sp22sp1="blat.exon.rheMac8Tohg19.filtered.txt"

sp1="hg19"
sp2="rheMac8"

blat.Sp1ToSp1 <- read.table(sp12sp1, as.is = T, header=T)
blat.Sp1ToSp2 <- read.table(sp12sp2, as.is = T,  header=T)
blat.Sp2ToSp1 <- read.table(sp22sp1, as.is = T, header=T)
blat.Sp2ToSp2 <- read.table(sp22sp2,as.is = T, header=T)

# the vectors of PID and PL where the threshold to choose from
IDs <- seq(0.8, 0.999, 0.01)
PLs <- seq(0.8, 0.999, 0.05)

setwd("~/lustrelyt/monkey-brain/merged-hic/RNA-Seq-new/XSAnno/step2_Blat/blatFilter")
# plot the number of exons against the PID and PL used
pdf("exonNum_vs_blat_shresholds_interSp.pdf", 8,6)
	exonNumTable.inter <- chooseThreshold.inter(blat.Sp1ToSp2, blat.Sp2ToSp1, IDs, PLs)
dev.off()

pdf("exonNum_vs_blat_shresholds_intraSp.pdf", 8,6)
	exonNumTable.intra <- chooseThreshold.intra(blat.Sp1ToSp1, blat.Sp2ToSp2, IDs, PLs)
dev.off()





# choose interID, interPL, intraID and intraPL, when maximum exon number reached.
interID <- 0.9
interPL <- 0.9
intraID <- 0.97
intraPL <- 0.95

blatFiltered <- blatFilter(blat.Sp1ToSp1, blat.Sp1ToSp2, blat.Sp2ToSp1, blat.Sp2ToSp2, interID, interPL, intraID, intraPL, sp1Name=sp1, sp2Name=sp2)
write.table(blatFiltered, paste("blatFiltered.", sp1, "To", sp2, ".txt", sep=""), quote=F, sep="\t", col.names=T, row.names=F)
write.table(cbind(blatFiltered[,c(2,4,5,1)], rep(1, nrow(blatFiltered)), blatFiltered[,3]), paste("blatFiltered.", sp1, "To", sp2, ".", sp1, ".bed", sep=""), quote=F, sep="\t", col.names=F, row.names=F)
write.table(cbind(blatFiltered[,c(6,8,9,1)], rep(1, nrow(blatFiltered)), blatFiltered[,7]),paste("blatFiltered.", sp1, "To", sp2, ".", sp2, ".bed", sep=""), quote=F, sep="\t", col.names=F, row.names=F)




interID <- 0.9
interPL <- 0.95
intraID <- 0.97
intraPL <- 0.95

blatFiltered <- blatFilter(blat.Sp1ToSp1, blat.Sp1ToSp2, blat.Sp2ToSp1, blat.Sp2ToSp2, interID, interPL, intraID, intraPL, sp1Name=sp1, sp2Name=sp2)
write.table(blatFiltered, paste("blatFiltered.interPL095.", sp1, "To", sp2, ".txt", sep=""), quote=F, sep="\t", col.names=T, row.names=F)
write.table(cbind(blatFiltered[,c(2,4,5,1)], rep(1, nrow(blatFiltered)), blatFiltered[,3]), paste("blatFiltered.", sp1, "To", sp2, ".", sp1, ".bed", sep=""), quote=F, sep="\t", col.names=F, row.names=F)
write.table(cbind(blatFiltered[,c(6,8,9,1)], rep(1, nrow(blatFiltered)), blatFiltered[,7]),paste("blatFiltered.", sp1, "To", sp2, ".", sp2, ".bed", sep=""), quote=F, sep="\t", col.names=F, row.names=F)
