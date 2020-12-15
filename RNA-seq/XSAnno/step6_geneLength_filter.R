
anno.hg <- fread("step3_simNGS/simFiltered.byDIM.hg19TorheMac8.hg19.bed", header = F)
anno.rh <- fread("step3_simNGS/simFiltered.byDIM.hg19TorheMac8.rheMac8.bed", header = F)
names(anno.hg) <- names(anno.rh) <- c("chr","start","end","eid","gid","tid","gn","tn","type","id","strand")

anno.hg.gr <- makeGRangesFromDataFrame(anno.hg, keep.extra.columns = T)
anno.rh.gr <- makeGRangesFromDataFrame(anno.rh, keep.extra.columns = T)
anno.hg.gr.ls <- split(anno.hg.gr, as.factor(anno.hg.gr$gid))
anno.rh.gr.ls <- split(anno.rh.gr, as.factor(anno.rh.gr$gid))

hg.glen <- lapply(anno.hg.gr.ls, function(gr){
  sum(width(gr))
})
rh.glen <- lapply(anno.rh.gr.ls, function(gr){
  sum(width(gr))
})



anno.pre.hg <- fread("step1_LiftOver/Ensembl.GRCh37.exons.new.withname.bed", header = F, sep = "\t")
names(anno.pre.hg) <-  c("chr","start","end","eid","gid","tid","gn","tn","type")
anno.pre.hg <- anno.pre.hg[which(anno.pre.hg$gid %in% anno.hg$gid), ]

anno.pre.hg.gr <- makeGRangesFromDataFrame(anno.pre.hg, keep.extra.columns = T)
anno.pre.hg.gr.ls <- split(anno.pre.hg.gr, as.factor(anno.pre.hg.gr$gid))

hg.pre.glen <- lapply(anno.pre.hg.gr.ls, function(gr){
  sum(width(gr))
})
hg.glen.ratio <- unlist(hg.glen) / unlist(hg.pre.glen)

loci.hg.short <- which(hg.glen <= 1000)
loci.rh.short <- which(rh.glen <= 1000)
loci.hg.change <- which(hg.glen.ratio <= 0)

loci.gene.filter <- unique(c(loci.hg.short, loci.rh.short, loci.hg.change))
gene.filter <- names(anno.hg.gr.ls)[loci.gene.filter]


anno.hg.filter <- anno.hg[!which(anno.hg$gid %in% gene.filter), ]
anno.rh.filter <- anno.rh[!which(anno.rh$gid %in% gene.filter), ]

