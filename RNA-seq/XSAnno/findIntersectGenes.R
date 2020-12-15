
ref <- fread("../public/Gencode.v10.fullTransExon.hg19TorheMac2.hg19.bed", header = F, sep = "\t")
gid <- substr(ref$V4,1,15)
gid <- unique(gid)


our <- fread("step3_simNGS/simFiltered.hg19TorheMac8.rheMac8.withname.bed", header = F)
loci <- our$V5 %in% gid
our.intersect <- our[loci, ]


write.table(our.intersect, "XSAnno.byself.Ensemble.hg19TorheMac8.rheMac8.withname.bed", 
            row.names = F, col.names = F, quote = F, sep = "\t")
