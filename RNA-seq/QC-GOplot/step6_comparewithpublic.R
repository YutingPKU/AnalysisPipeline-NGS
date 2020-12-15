
pub <- read.xlsx("results/aan3456-TablesS1-S11.xlsx", startRow = 3, rowNames = F, colNames = T, sheet = 3)
gn <- pub$X1[-1]
gn <- lapply(gn, function(chr){
  substr(chr,17,40)
})
gn <- unlist(gn)

our <- read.xlsx("results/hg.rh.orthologs.DEG.xlsx", rowNames = T, colNames = T)
gn.o <- our$id

comm <- intersect(gn, gn.o)
comm.our <- our[match(comm, our$id),]



load("data/Homo_sapiens.Ensembl.GRCh37.gene.GRanges.rds")
loci <- match(our$id, my_gr$gn)
our$gid <- my_gr$id[loci]
our.up <- our[which(our$log2FoldChange > 0),]
our.dw <- our[which(our$log2FoldChange < 0),]
