############################
# plot public cell line heatmap to identify A/B results
############################
source("/lustre/user/liclab/liuyt/monkey-brain/rheMac-brain/Basic-Data/ABcompartment/plotAB.R")
setwd("/lustre/user/liclab/lirf/Project/hic/disease/data/GM12878_combined/matrix_500kb")
#setwd("/lustre/user/liclab/lirf/Project/hic/disease/data/IMR90/matrix_500k")
setwd("/lustre/user/liclab/lirf/Project/hic/disease/data/HMEC/matrix_500kb")
pdf("~/lustrelyt/monkey-brain/human-brain/Basic-Data/ABcompartment/ABplot/heatmap.HMEC.pdf")
for(i in c(1:23)){
  mat = fread(paste0("chr",i, "_500kb_rawmatrix.txt"))
  #mat = fread(paste0("chr", i, "_chr", i, "_500kb_normalmatrix.txt"))
  mat = mat[,-1]
  mat[is.na(mat)] = 0
  mat = as.matrix(mat)
  mat[mat >= quantile(mat, 0.95)] = quantile(mat, 0.95) # set maximum contact value
  resolution = 500000
  case_name = paste0("HMEC chr", i)
  plotheatmap(mat, case_name, resolution) # plot hic heatmap
}

dev.off()

