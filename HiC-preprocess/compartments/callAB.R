## Yuting Liu
## call AB compartments for matrix from HiC-Pro
## step1: construct HiTC list object from whole genome contact matrix in triple list format
## step2: call A/B by pca.hic
## step3: correct A/B assign results
## step4: test A/B assign results by gene density
## step5: test A/B assign results by plotting heatmap and coresspongding compartmetns, gene density

## step1 construct HiTC object 
source("/lustre/user/liclab/liuyt/monkey-brain/rheMac-brain/Basic-Data/ABcompartment/plotAB.R")
library(data.table)
library(HiTC)
library(refGenome)
setwd("/lustre/user/liclab/liuyt/monkey-brain/merged-hic/basic-data/ABcompartment/")
matrix <- "/lustre/user/liclab/liuyt/monkey-brain/merged-hic/basic-data/ABcompartment/GZ_500000.matrix"
xgi <- ygi <- "/lustre/user/liclab/liuyt/monkey-brain/07456AB/basicdata/ABcompartments/A-CP_500000_abs.bed"
hic <- importC(matrix, xgi, ygi)


## step2 extract individual chromosomes interaction maps and call A/B comparments
gtffile = "/lustre/user/liclab/liuyt/monkey-brain/ref_genome/ensemble/Mmul_8.0.1/Macaca_mulatta.Mmul_8.0.1.89.chr.gtf"
setwd("/lustre/user/liclab/liuyt/monkey-brain/ref_genome/ensemble/Mmul_8.0.1/")
ens <- ensemblGenome()
read.gtf(ens, "Macaca_mulatta.Mmul_8.0.1.89.chr.gtf")
my_gene <- getGenePositions(ens)
my_gr <- with(my_gene, GRanges(paste0("chr", seqid), IRanges(start, end), strand, id = gene_id))
setwd("/lustre/user/liclab/liuyt/monkey-brain/merged-hic/basic-data/ABcompartment/")



hic.intra <- hic[isIntraChrom(hic)] ## get all cis maps
hic.intra
hic.intra <- hic.intra[1:21] # contain chromosome 1-20 and X

pca.intra <- lapply(hic.intra, function(xx){
  #xx.iced <- normICE_our(xx, max_iter = 100, eps = 1e-6)
  #pca <- pca.hic_our(xx, normPerExpected = TRUE, npc = 1, gene.gr = my_gr)
  pca <- pca.hic(xx, normPerExpected = T, npc = 1, gene.gr = my_gr)
  return(pca)
})

compartments <- data.frame(group=character(), group_name=character(), seqnames=character(), start=integer(), end=integer(), width=integer(), strand=character(),
                           score=double(), genedens=integer(), ccompartments=character(),stringsAsFactors=FALSE)

for(i in c(1:20,"X")){
  chrom <- paste0("chr",i,"chr",i)
  pca.chrom <- data.frame(pca.intra[chrom])
  pca.chrom <- pca.chrom[which(pca.chrom[,1] == "1"),]
  names(pca.chrom) <- c("group","group_name","seqnames","start","end","width","strand","score","genedens","ccompartments")
  pos.den <- mean(pca.chrom$genedens[which(pca.chrom$score >= 0)])
  neg.den <- mean(pca.chrom$genedens[which(pca.chrom$score < 0)])
  print(paste0(chrom, " pos ", pos.den, " neg ", neg.den))
  if( pos.den < neg.den){
    pca.chrom$score = -pca.chrom$score   # correct A/B assign: HiTC pca.hic may make an error because it counts gene sum not gene density
  }
  pca.chrom$ccompartments[which(pca.chrom$score >= 0)] = "A"
  pca.chrom$ccompartments[which(pca.chrom$score < 0)] = "B"
  compartments <- rbind(compartments, pca.chrom)
}



## step4 save results
write.table(compartments, "results/GZ-500kb-ABcompartments-PC1-Raw.txt",
            sep = '\t', row.names = F, col.names = T, quote = F)