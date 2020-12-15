FindAregion <- function(comp, bed, type){
  region <- data.frame(chrom = character(), start = integer(), end = integer(), compartments = character(), stringsAsFactors=FALSE) 
  h <- 1
  t <- 1
  while ( h < length(comp)){
    if(!is.na(comp[h]) && comp[h] == type){
      t <- h
      while ( !is.na(comp[h+1]) && comp[h+1] == type){
        h <- h + 1
      }
      aregion <- data.frame(chrom = bed$chrom[t], start = bed$start[t], end = bed$end[h], compartments = type)
      region <- rbind(region, aregion)
    } 
    h <- h + 1
  }
  return(region)
}

regionlen <- function(comp, celltype){
  names(comp) = c("chrom","start","end","pc")
  comp$cc = rep("NA", length(comp$chr))
  comp$cc[which(comp$pc > 0)] = "A"
  comp$cc[which(comp$pc < 0)] = "B"
  regionCP <- data.frame(chrom = character(), start = integer(), end = integer(), stringsAsFactors=FALSE) 
  for(i in c(1:22, "X")){
    chr <- paste0("chr", i)
    compt = as.character(comp$cc[which(comp$chrom == chr)])
    bed = comp[which(comp$chrom == chr), 1:3]
    regionA <- FindAregion(compt, bed, "A")
    regionB <- FindAregion(compt, bed, "B")
    regionCP <- rbind(regionCP, regionA)
    regionCP <- rbind(regionCP, regionB)
  }
  cp <- makeGRangesFromDataFrame(regionCP, keep.extra.columns = T)
  cpmat <- data.frame("CellType" = rep(celltype, length(cp)), "Length" = width(cp), "Compartments" = mcols(cp)[,1])
  return(cpmat)
}

setwd("/lustre/user/liclab/liuyt/SV-3Dgenome/brain-hicdata/Compartment_primary_cohort/")
files <- list.files("/lustre/user/liclab/liuyt/SV-3Dgenome/brain-hicdata/Compartment_primary_cohort/", pattern = "*bedGraph")

width.mat <- data.frame(CellType = character(), Length = integer(), Compartments = character(), stringsAsFactors=FALSE) 
for(i in 1:length(files)){
  celltype <- substr(files[i],1,(nchar(files[i])-15))
  comp <- fread(files[i])
  mat <- regionlen(comp, celltype)
  width.mat <- rbind(width.mat, mat)
}
setwd("~/lustrelyt/monkey-brain/human-brain/Basic-Data/")
GZ <- fread("ABcompartment/compartments-500k-byrawmat/GZ_500k_wholegenome_compartments.txt", stringsAsFactors = F)
CP <- fread("ABcompartment/compartments-500k-byrawmat/CP_500k_wholegenome_compartments.txt", stringsAsFactors = F)
cp <- CP[,c(3:5,8)]
cp.mat <- regionlen(cp, "CP")
gz <- GZ[,c(3:5,8)]
gz.mat <- regionlen(gz, "GZ")
width.mat <- rbind(width.mat, cp.mat, gz.mat)

width.mat$Length <- as.numeric(as.character(width.mat$Length))
width.mat$Length <- log10(width.mat$Length)

g = 
  ggplot(width.mat)+
  geom_boxplot(data = width.mat, aes(x=factor(Compartments), y=Length, fill=factor(CellType)), notch = T, outlier.shape = NA, lwd = 1.5) +
  xlab("Compartments") +
  ylab("log10 (Length)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_blank())+
  theme(axis.text.x = element_text(size=24),
        axis.text.y = element_text(size=24),  
        axis.title.x = element_text(size=28),
        axis.title.y = element_text(size=28),
        axis.line.y = element_line(colour = "black", size = 1.5),
        axis.ticks.x=element_blank(),
        axis.ticks.y = element_line(size = 1.5),
        legend.text=element_text(size=18),
        plot.title = element_text(size=22, hjust = 0.5))+
  scale_x_discrete(labels=c("A","B")) 

g
