library(Rsubread)
library(data.table)

setwd("~/lustrelyt/monkey-brain/merged-hic/RNA-Seq-new/XSAnno")
bam.files <- list.files(path = "step3_simNGS/bamfiles/human/", full.names = T, pattern = "bam")

anno <- fread("step2_Blat/blatFilter/blatFiltered.hg19TorheMac8.hg19.withname.bed", header = F, sep = "\t")
anno <- anno[,c(5,4,1:3,11)]
names(anno) <- c("GeneID","ExonID","Chr","Start","End", "Strand")
#anno$Chr <- substr(anno$Chr,4,20)

count <- featureCounts(bam.files, annot.ext = anno, isGTFAnnotationFile = F, maxFragLength = 1000,
                       allowMultiOverlap = T,
                       isPairedEnd = T, requireBothEndsMapped = T, nthreads = 5, useMetaFeatures = F)
saveRDS(count, file = "step3_simNGS/hg19.simluated.RNASeq.featureCount.byexon.allMO.rds")

