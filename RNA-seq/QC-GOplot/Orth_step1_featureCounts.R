library(Rsubread)
library(data.table)

setwd("~/lustrelyt/monkey-brain/merged-hic/RNA-Seq-new/")
bam.files <- list.files(path = "data/bamfiles/human/", full.names = T, pattern = "bam")

anno <- fread("XSAnno/XSAnno.byself.Ensemble.hg19TorheMac8.hg19.withname.bed", header = F, sep = "\t")
anno <- anno[,c(5,4,7,1:4)]
names(anno) <- c("GeneID","ExonID","GeneSy","Chr","Start","End", "Strand")
anno$Strand <- "."
#anno$Chr <- substr(anno$Chr,4,20)

count <- featureCounts(bam.files, annot.ext = anno, isGTFAnnotationFile = F, maxFragLength = 1000,
                       isPairedEnd = T, requireBothEndsMapped = T, nthreads = 6, useMetaFeatures = T)
saveRDS(count, file = "data/human.feature.count.byXSAnno.genes.rds")




