
orth <- fread("public/Gencode.v10.fullTransExon.hg19TorheMac2.hg19.bed", header = F, sep = "\t")
orth.gid <- unique(substr(orth$V4,1,15))
orth.gsm <- fread("public/Gencode.v10.orthlogous.gene.symbol.txt", header = F)
orth.rh <- fread("public/Gencode.v10.fullTransExon.hg19TorheMac2.rheMac2.bed", header = F, sep = "\t")
orth.gid.rh <- unique(orth.rh$V4,1,15)


######## gene annotation data
load("data/Homo_sapiens.Ensembl.GRCh37.gene.GRanges.rds")
hg.ref <- my_gr
load("data/Macaca_mulatta.Mmul_8.0.1.89.chr.gene.GRanges.rds")
rh.ref <- my_gr

######## select orthlogous genes
comm.hg <- intersect(orth.gsm, hg.ref$gn)
comm.rh <- intersect(orth.gsm, rh.ref$gn)
comm <- intersect(comm.hg, comm.rh)

hg.mat <- readRDS("data/human.exonbygene.summarizeOverlaps.rds")
rh.mat <- readRDS("data/monkey-exonbygene-summaried-expriment-rmd.rds")
