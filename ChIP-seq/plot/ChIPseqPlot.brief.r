##plot ChIP-Seq read depth around TSS and genebody

library(ggplot2)
library(reshape2)
library(pheatmap)
library(RColorBrewer)

############### TSS plot 
readTable = function(path){
  P1_Input = read.table(path, stringsAsFactors = F , header = F, blank.lines.skip = TRUE)
  rownames(P1_Input) = make.names(P1_Input$V1, unique = T)
  P1_Input = P1_Input[, -1]
  P1_Input = colMeans(P1_Input)
  P1_Input = as.vector(P1_Input)
  return(P1_Input)
}

CP_H3K27ac = readTable("10084A-CP-H3K27AC-lib2_R1.fq.gz.sam.sort.rmd.bam/tss.datatable")
CP_H3K4me3 = readTable("10084A-CP-H3K4me3-lib2_R1.fq.gz.sam.sort.rmd.bam/tss.datatable")
CP_H3K27me3 = readTable("10084A-CP-H3K27me3-lib2_R1.fq.gz.sam.sort.rmd.bam/tss.datatable") #inactive

GZ_H3K4me3 = readTable("10084A-GZ-H3K4me3-lib2_R1.fq.gz.sam.sort.rmd.bam/tss.datatable")

MEF_H3K27ac = readTable("10084A-MEF-H3K27AC-lib2_R1.fq.gz.sam.sort.rmd.bam/tss.datatable")
MEF_H3K4me3 = readTable("10084A-MEF-H3K4me3-lib2_R1.fq.gz.sam.sort.rmd.bam/tss.datatable")
MEF_H3K27me3 = readTable("10084A-MEF-H3K27me3-lib2_R1.fq.gz.sam.sort.rmd.bam/tss.datatable")#inactive

PFC_H3K27ac = readTable("public/tss.datatable")
CB_H3K27ac = readTable("public2/tss.datatable")
TN_H3k27ac = readTable("public3/tss.datatable")

#mat = cbind(CP_H3K27ac,CP_H3K27me3, CP_H3K4me3)
mat = cbind(CP_H3K27ac, MEF_H3K27ac)
mat = cbind(CP_H3K4me3, GZ_H3K4me3, MEF_H3K4me3)
mat = cbind(CP_H3K27me3, MEF_H3K27me3)
mat = cbind(CP_H3K27ac,CP_H3K27me3, CP_H3K4me3)
mat = cbind(MEF_H3K27ac, MEF_H3K27me3, MEF_H3K4me3)
mat = cbind(GZ_H3K4me3)
mat = cbind(MEF_H3K4me3)
mat = cbind(PFC_H3K27ac,CB_H3K27ac,TN_H3k27ac)

mat.melt = melt(mat)
colnames(mat.melt) = c('x', 'Type', 'RPM')
mat.melt$Sample = factor(mat.melt$Type)

#png("10084A_TSS_4k.png", width = 700, height = 500)
p <- ggplot(mat.melt, aes(x=x, y=RPM, col=Sample)) + geom_line(size=2) + 
  geom_vline(xintercept = 101, colour="#E69F00", size=1.2,alpha=0.8) +
  theme_classic(base_size = 20) +
  theme( text = element_text(size=20),
         #panel.background=element_blank(),
         axis.ticks.length = unit(.25, "cm"),
         #legend.position = "top",
         legend.key = element_rect(fill = "white"))+
  labs (y='reads per million') +
  scale_colour_brewer(palette="Dark2")
p + scale_x_continuous(name="", limits=c(0, 201),breaks=c(0, 101, 201), labels=c('-4000','TSS','+4000'))
#dev.off()





#################### heatmap #####################
readTable_hm = function(path){
  P1_Input = read.table(path, stringsAsFactors = F , header = F, blank.lines.skip = TRUE)
  rownames(P1_Input) = make.names(P1_Input$V1, unique = T)
  P1_Input = P1_Input[, -1]
  return(P1_Input)
}

input = readTable_hm("10084A-CP-H3K4me3-lib2_R1.fq.gz.sam.sort.rmd.bam/tss.datatable")
input = input[order(rowMeans(input),decreasing=T),]
png("CP_H3K4me3.pheatmap.png", width = 400, height = 1000)
pheatmap(log(input+ 0.1), cluster_rows = F, cluster_cols = F, border_color = NA,
         show_rownames = F, show_colnames = F,
         color = colorRampPalette(brewer.pal(n = 7, name = "Greens"))(100))
dev.off()


