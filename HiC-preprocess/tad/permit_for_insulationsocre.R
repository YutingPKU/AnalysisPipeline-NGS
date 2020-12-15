library(data.table)
bed = fread("/lustre/user/liclab/liuyt/monkey-brain/11002B/Basicdata/matrix/CP/raw/40000/CP_40000_abs.bed", sep = "\t", header = F)
names(bed) = c("chr","start","end","id")
setwd("/lustre/user/liclab/liuyt/monkey-brain/07456AB/basicdata/matrix/07456B-GZ/iced/40000/")
for(i in c(1:20, "X")){
  file = paste0("modify-07456B-GZ_40000_iced_chr",i,"_dense.matrix")
  mat = fread(file)
  #mat = mat[, -c(1:3)]
  mat$V1 = NULL
  mat$V2 = NULL
  mat$V3 = NULL
  print(dim(mat))
  cbed = bed[which(bed$chr == paste0("chr",i)),]
  name = lapply("num", paste, "|rheMac8|", as.character(cbed$chr), ":", cbed$start, "-", cbed$end, sep = "")[[1]]
  name1 = name
  name1[1] = paste0("\t",name1[1],collapse = "")
  rownames(mat) = name
  colnames(mat) = name1
  print(i)
  write.table(mat, paste0("/lustre/user/liclab/liuyt/monkey-brain/07456AB/basicdata/matrix/07456B-GZ/modify_insulation_mat/",file,"_insulation.mat"), sep = "\t", quote = F)
  
}
