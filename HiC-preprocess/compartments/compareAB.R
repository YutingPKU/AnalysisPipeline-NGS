############################
# compare A/B btw liver, MEF, CP and GZ
# step1: heatmap and statistical comparsion of A/B switch ratio
# step2: statistical A/B border
############################

setwd("~/lustrelyt/monkey-brain/human-brain/Basic-Data/")
#hESC <- fread("ABcompartment/compartments-500k-byrawmat/hESC_500k_wholegenome_compartments.txt", stringsAsFactors = F)
LG <- fread("ABcompartment/compartments-500k-byrawmat/LG_500k_wholegenome_compartments.txt", stringsAsFactors = F)
GZ <- fread("ABcompartment/compartments-500k-byrawmat/GZ_500k_wholegenome_compartments.txt", stringsAsFactors = F)
CP <- fread("ABcompartment/compartments-500k-byrawmat/CP_500k_wholegenome_compartments.txt", stringsAsFactors = F)
CO <- fread("ABcompartment/compartments-500k-byrawmat/CO_500k_wholegenome_compartments_corrected.txt", stringsAsFactors = F)
HC <- fread("ABcompartment/compartments-500k-byrawmat/HC_500k_wholegenome_compartments.txt", stringsAsFactors = F)

## step1 heatmap of PC1
mat <- cbind("seqnames" = CP$seqnames, "start"=CP$start, "end"=CP$end, "CP-PCA1" = CP$score, "CP-AB" = as.character( CP$ccompartments), 
             "GZ-PCA1" = GZ$score, "GZ-AB" = as.character(GZ$ccompartments), "hESC-PCA1" = hESC$score,"hESC-AB" = as.character(hESC$ccompartments), 
             "LG-PCA1" = LG$score, "LG-AB" = as.character(LG$ccompartments), "CO-PCA1" = CO$score, "CO-AB" = as.character(CO$ccompartments), 
             "HC-PCA1" = HC$score, "HC-AB" = as.character(HC$ccompartments))
mat <- data.frame(mat)
mat$CP.PCA1 = as.numeric(as.character(mat$CP.PCA1))
mat$GZ.PCA1 = as.numeric(as.character(mat$GZ.PCA1))
mat$MEF.PCA1= as.numeric(as.character(mat$MEF.PCA1))
mat$Liver.PCA1 = as.numeric(as.character(mat$Liver.PCA1))


############################################# step2 statistical compare of A/B border
AbBoder <- function(comp){
  num <- 0
  comp <- as.character(comp)
  for(i in 1:length(comp)){
    if (!is.na(comp[i]) && !is.na(comp[i+1])) {
      if ((comp[i] == "A" && comp[i+1] == "B") ||  (comp[i] == "B" && comp[i+1] == "A")) {
        num <- num + 1
      }
    }
  }
  return(num) 
}
cp <- AbBoder(mat$CP.AB)
gz <- AbBoder(mat$GZ.AB)
hesc <- AbBoder(mat$hESC.AB)
lg <- AbBoder(mat$LG.AB)
co <- AbBoder(mat$CO.AB)
hc <- AbBoder(mat$HC.AB)

########################################## step3 find A/B regions length and number change
## define A region: a region contains continuous A compartments 

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
      aregion <- data.frame(chrom = bed$seqnames[t], start = bed$start[t], end = bed$end[h], compartments = type)
      region <- rbind(region, aregion)
    } 
    h <- h + 1
  }
  return(region)
}

regionCP <- data.frame(chrom = character(), start = integer(), end = integer(), stringsAsFactors=FALSE) 
for(i in c(1:22, "X")){
  chr <- paste0("chr", i)
  comp = as.character(mat$CP.AB[which(mat$seqnames == chr)])
  bed = mat[which(mat$seqnames == chr), 1:3]
  regionA <- FindAregion(comp, bed, "A")
  regionB <- FindAregion(comp, bed, "B")
  regionCP <- rbind(regionCP, regionA)
  regionCP <- rbind(regionCP, regionB)
}

regionGZ <- data.frame(chrom = character(), start = integer(), end = integer(), stringsAsFactors=FALSE) 
for(i in c(1:22, "X")){
  chr <- paste0("chr", i)
  comp = as.character(mat$GZ.AB[which(mat$seqnames == chr)])
  bed = mat[which(mat$seqnames == chr), 1:3]
  regionA <- FindAregion(comp, bed, "A")
  regionB <- FindAregion(comp, bed, "B")
  regionGZ <- rbind(regionGZ, regionA, regionB) 
}

regionhESC <- data.frame(chrom = character(), start = integer(), end = integer(), stringsAsFactors=FALSE) 
for(i in c(1:22, "X")){
  chr <- paste0("chr", i)
  comp = as.character(mat$hESC.AB[which(mat$seqnames == chr)])
  bed = mat[which(mat$seqnames == chr), 1:3]
  regionA <- FindAregion(comp, bed, "A")
  regionB <- FindAregion(comp, bed, "B")
  regionhESC <- rbind(regionhESC, regionA, regionB) 
}

regionLG <- data.frame(chrom = character(), start = integer(), end = integer(), stringsAsFactors=FALSE) 
for(i in c(1:22, "X")){
  chr <- paste0("chr", i)
  comp = as.character(mat$LG.AB[which(mat$seqnames == chr)])
  bed = mat[which(mat$seqnames == chr), 1:3]
  regionA <- FindAregion(comp, bed, "A")
  regionB <- FindAregion(comp, bed, "B")
  regionLG <- rbind(regionLG, regionA, regionB) 
}

regionCO <- data.frame(chrom = character(), start = integer(), end = integer(), stringsAsFactors=FALSE)
for(i in c(1:20, "X")){
  chr <- paste0("chr", i)
  comp = as.character(mat$CO.AB[which(mat$seqnames == chr)])
  bed = mat[which(mat$seqnames == chr), 1:3]
  regionA <- FindAregion(comp, bed, "A")
  regionB <- FindAregion(comp, bed, "B")
  regionCO <- rbind(regionCO, regionA, regionB) 
}

cp <- makeGRangesFromDataFrame(regionCP, keep.extra.columns = T)
cpmat <- data.frame("CellType" = rep("CP", length(cp)), "Length" = width(cp), "Compartments" = mcols(cp)[,1])
gz <- makeGRangesFromDataFrame(regionGZ, keep.extra.columns = T)
gzmat <- data.frame("CellType" = rep("GZ", length(gz)), "Length" = width(gz), "Compartments" = mcols(gz)[,1])
mef <- makeGRangesFromDataFrame(regionhESC, keep.extra.columns = T)
mefmat <- data.frame("CellType" = rep("hESC", length(mef)), "Length" = width(mef), "Compartments" = mcols(mef)[,1])
liver <- makeGRangesFromDataFrame(regionLG, keep.extra.columns = T)
livermat <- data.frame("CellType" = rep("Lung", length(liver)), "Length" = width(liver), "Compartments" = mcols(liver)[,1])
co <- makeGRangesFromDataFrame(regionCO, keep.extra.columns = T)
comat <- data.frame("CellType" = rep("CO", length(co)), "Length" = width(co), "Compartments" = mcols(co)[,1])

ggplot(data, aes(x=variety, y=note, fill=treatment)) + 
  geom_boxplot()

width.mat <- rbind(livermat,  gzmat, cpmat)
width.mat$Length <- as.numeric(as.character(width.mat$Length))
width.mat$Length <- log10(width.mat$Length)


pdf("ABplot/ABregion_length_change.pdf", width = 10, height = 6)
g = #ggplot(width.mat, aes(x=Compartments, y=Length, fill=CellType))+
  ggplot(width.mat)+
  geom_boxplot(data = width.mat, aes(x=factor(Compartments), y=Length, fill=factor(CellType)), notch = T, outlier.shape = NA, lwd = 1.5) +
  #geom_boxplot(notch = F, outlier.shape = NA, lwd = 2, show.legend = T)+
  #stat_boxplot(geom ='errorbar')+
  #stat_summary(fun.y=mean, geom="point", shape=23, size=4)+
  #coord_cartesian(ylim = c(4.5,8))+
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
  #theme(legend.position="none")+
  scale_x_discrete(labels=c("A","B")) 
#scale_fill_manual(values=c("dodgerblue3", "forestgreen", "darkorange2", "firebrick2"))
#scale_y_continuous(breaks = seq(4.5,8,0.5), labels = seq(4.5,8,0.5))
g
dev.off()


length(which(cp$compartments == "A"))
length(which(cp$compartments == "B"))
length(which(gz$compartments == "A"))
length(which(gz$compartments == "B"))