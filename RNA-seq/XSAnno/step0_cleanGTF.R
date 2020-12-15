
ref <- fread("exons.gtf", header = F)
gid <- lapply(ref$V9, function(vec){
  vec <- gsub(" ", '', vec)
  vec <- noquote(vec)
  ls <- unlist(strsplit(vec, ";"))
  
  eid <- ls[grep("exon_id", ls)]
  eid <- noquote(eid)
  eid <- unlist(strsplit(eid,"\""))[2]
  
  gid <- ls[grep("gene_id", ls)]
  gid <- noquote(gid)
  gid <- unlist(strsplit(gid,"\""))[2]
  
  gn <- ls[grep("gene_name", ls)]
  gn <- noquote(gn)
  gn <- unlist(strsplit(gn,"\""))[2]
  
  tid <- ls[grep("transcript_id", ls)]
  tid <- noquote(tid)
  tid <- unlist(strsplit(tid,"\""))[2]
  
  tn <- ls[grep("transcript_name", ls)]
  tn <- noquote(tn)
  tn <- unlist(strsplit(tn,"\""))[2]
  
  type <- ls[grep("gene_biotype", ls)]
  type <- noquote(type)
  type <- unlist(strsplit(type,"\""))[2]
  
  ann <- paste(eid,gid,tid,gn,tn,type, sep = "|")
  return(ann)
})
ref$ann <- unlist(gid)

ebed <- ref[,c(1,4,5,10)]
names(ebed) <- c("chr","start","end","ann")
ebed$chr <- paste0("chr", ebed$chr)
ebed.gr <- makeGRangesFromDataFrame(ebed, keep.extra.columns = T)

write.table(ebed, "Ensembl.GRCh37.exons.new.bed", row.names = F, col.names = F, quote = F, sep = "\t")

query <- fread("H2R/hg19.hg19_rheMac8.liftOver.bed", header = F, sep = "\t")
names(query) <- c("chr","start","end","ann")
query.gr <- makeGRangesFromDataFrame(query, keep.extra.columns = T)

comm <- findOverlaps(query.gr, ebed.gr, type = "equal")
