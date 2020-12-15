# CHiP seq plot
setwd("/lustre/user/liclab/zhangc/Taolab/guan/ChIP-seq/plot2")

library(ggplot2)
library(reshape2)

readTable = function(path){
  P1_Input = read.table(path, stringsAsFactors = F , header = F, blank.lines.skip = TRUE)
  rownames(P1_Input) = make.names(P1_Input$V1, unique = T)
  P1_Input = P1_Input[, -1]
  P1_Input = colMeans(P1_Input)
  P1_Input = as.vector(P1_Input)
  return(P1_Input)
}

P1_Input = readTable("MEF-P1-Input.TSSplot/tss.datatable")

P4_Input = readTable("MEF-P4-Input.TSSplot/tss.datatable")

P7_Input = readTable("MEF-P7-Input.TSSplot/tss.datatable")

P1_H3K4me1 = readTable("MEF-P1-H3K4me1.TSSplot/tss.datatable")

P4_H3K4me1 = readTable("MEF-P4-H3K4me1.TSSplot/tss.datatable")

P7_H3K4me1 = readTable("MEF-P7-H3K4me1.TSSplot/tss.datatable")

P1_H3K27Ac = readTable("MEF-P1-H3K27Ac.TSSplot/tss.datatable")

P4_H3K27Ac = readTable("MEF-P4-H3K27Ac.TSSplot/tss.datatable")

P7_H3K27Ac = readTable("MEF-P7-H3K27Ac.TSSplot/tss.datatable")




#mat = cbind(P1_Input, P4_Input, P7_Input, P1_H3K4me1, P4_H3K4me1, P7_H3K4me1, P1_H3K27Ac, P4_H3K27Ac, P7_H3K27Ac)
mat = cbind(P1_H3K4me1, P4_H3K4me1, P7_H3K4me1)
mat = cbind(P1_H3K27Ac, P4_H3K27Ac, P7_H3K27Ac)

mat.melt = melt(mat)
colnames(mat.melt) = c('x', 'Type', 'RPM')
mat.melt$Sample = factor(mat.melt$Type)

png("H3K4me1_TSS_10k.png", width = 700, height = 500)
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
p + scale_x_continuous(name="", limits=c(0, 201),breaks=c(0, 101, 201), labels=c('-10k','TSS','10k'))
dev.off()




readTable = function(path){
  P1_Input = read.table(path, stringsAsFactors = F , header = F, blank.lines.skip = TRUE)
  rownames(P1_Input) = make.names(P1_Input$V1, unique = T)
  P1_Input = P1_Input[, -1]
  P1_Input = colMeans(P1_Input)
  P1_Input = as.vector(P1_Input)
  return(P1_Input)
}



############ gene body ###########

P1_Input = readTable("MEF-P1-Input.TSSplot/genebody.datatable")

P4_Input = readTable("MEF-P4-Input.TSSplot/genebody.datatable")

P7_Input = readTable("MEF-P7-Input.TSSplot/genebody.datatable")

P1_H3K4me1 = readTable("MEF-P1-H3K4me1.TSSplot/genebody.datatable")

P4_H3K4me1 = readTable("MEF-P4-H3K4me1.TSSplot/genebody.datatable")

P7_H3K4me1 = readTable("MEF-P7-H3K4me1.TSSplot/genebody.datatable")

P1_H3K27Ac = readTable("MEF-P1-H3K27Ac.TSSplot/genebody.datatable")

P4_H3K27Ac = readTable("MEF-P4-H3K27Ac.TSSplot/genebody.datatable")

P7_H3K27Ac = readTable("MEF-P7-H3K27Ac.TSSplot/genebody.datatable")




#mat = cbind(P1_Input, P4_Input, P7_Input, P1_H3K4me1, P4_H3K4me1, P7_H3K4me1, P1_H3K27Ac, P4_H3K27Ac, P7_H3K27Ac)
mat = cbind(P1_H3K4me1, P4_H3K4me1, P7_H3K4me1)
mat = cbind(P1_H3K27Ac, P4_H3K27Ac, P7_H3K27Ac)


mat.melt = melt(mat)
colnames(mat.melt) = c('x', 'Sample', 'RPM')
mat.melt$Sample = factor(mat.melt$Sample)

pdf("H3K4me1_geneBody.pdf", width = 700, height = 500)
p <- ggplot(mat.melt, aes(x=x, y=RPM, col=Sample)) + geom_line(size=1.2) +
  geom_vline(xintercept = 50, colour="gray", size=1.2,alpha=0.8) +
  geom_vline(xintercept = 150, colour="gray", size=1.2,alpha=0.8) +
  theme_bw()+
  theme( text = element_text(size=40),
         #panel.background=element_blank(),
         axis.ticks.length = unit(.25, "cm"),
         #legend.position = "none",
         legend.key = element_rect(fill = "white"))+
  labs (y='reads per million')
p + scale_x_continuous(name="", limits=c(0, 200),breaks=c(0, 50, 150, 200), labels=c('-2k','TSS', 'TTS', '2k'))
p
dev.off()



#################### heatmap #####################
readTable = function(path){
  P1_Input = read.table(path, stringsAsFactors = F , header = F, blank.lines.skip = TRUE)
  rownames(P1_Input) = make.names(P1_Input$V1, unique = T)
  P1_Input = P1_Input[, -1]
  return(P1_Input)
}

P1_Input = readTable("MEF-P1-Input.TSSplot/tss.datatable")

P4_Input = readTable("MEF-P4-Input.TSSplot/tss.datatable")

P7_Input = readTable("MEF-P7-Input.TSSplot/tss.datatable")

P1_H3K4me1 = readTable("MEF-P1-H3K4me1.TSSplot/tss.datatable")

P4_H3K4me1 = readTable("MEF-P4-H3K4me1.TSSplot/tss.datatable")

P7_H3K4me1 = readTable("MEF-P7-H3K4me1.TSSplot/tss.datatable")

P1_H3K27Ac = readTable("MEF-P1-H3K27Ac.TSSplot/tss.datatable")

P4_H3K27Ac = readTable("MEF-P4-H3K27Ac.TSSplot/tss.datatable")

P7_H3K27Ac = readTable("MEF-P7-H3K27Ac.TSSplot/tss.datatable")


P1_H3K4me1 = P1_H3K4me1[order(rowMeans(P1_H3K4me1),decreasing=T),]
genes = rownames(P1_H3K4me1)
P4_H3K4me1 = P4_H3K4me1[genes, ]
P7_H3K4me1 = P7_H3K4me1[genes, ]

P1_H3K27Ac = P1_H3K27Ac[genes, ]
P4_H3K27Ac = P4_H3K27Ac[genes, ]
P7_H3K27Ac = P7_H3K27Ac[genes, ]

png("H3K4me1_heatmap_P1.png", width = 400, height = 1000)
pheatmap(log(P1_H3K4me1 + 0.1), cluster_rows = F, cluster_cols = F, border_color = NA,
         show_rownames = F, show_colnames = F,
         color = colorRampPalette(brewer.pal(n = 7, name = "YlOrRd"))(100))
dev.off()
png("H3K4me1_heatmap_P4.png", width = 400, height = 1000)
pheatmap(log(P4_H3K4me1 + 0.1), cluster_rows = F, cluster_cols = F, border_color = NA,
         show_rownames = F, show_colnames = F,
         color = colorRampPalette(brewer.pal(n = 7, name = "YlOrRd"))(100))
dev.off()
png("H3K4me1_heatmap_P7.png", width = 400, height = 1000)
pheatmap(log(P7_H3K4me1 + 0.1), cluster_rows = F, cluster_cols = F, border_color = NA,
         show_rownames = F, show_colnames = F,
         color = colorRampPalette(brewer.pal(n = 7, name = "YlOrRd"))(100))
dev.off()

png("H3K27Ac_heatmap_P1.png", width = 400, height = 1000)
pheatmap(log(P1_H3K27Ac + 0.1), cluster_rows = F, cluster_cols = F, border_color = NA,
         show_rownames = F, show_colnames = F,
         color = colorRampPalette(brewer.pal(n = 7, name = "YlOrRd"))(100))
dev.off()
png("H3K27Ac_heatmap_P4.png", width = 400, height = 1000)
pheatmap(log(P4_H3K27Ac + 0.1), cluster_rows = F, cluster_cols = F, border_color = NA,
         show_rownames = F, show_colnames = F,
         color = colorRampPalette(brewer.pal(n = 7, name = "YlOrRd"))(100))
dev.off()
png("H3K27Ac_heatmap_P7.png", width = 400, height = 1000)
pheatmap(log(P7_H3K27Ac + 0.1), cluster_rows = F, cluster_cols = F, border_color = NA,
         show_rownames = F, show_colnames = F,
         color = colorRampPalette(brewer.pal(n = 7, name = "YlOrRd"))(100))
dev.off()


pheatmap(P7_H3K27Ac, cluster_rows = F, cluster_cols = F, border_color = NA,
         show_rownames = F, show_colnames = F,
         color = colorRampPalette(brewer.pal(n = 7, name = "YlOrRd"))(100))


############ repeat tRNA ##########


P1_H3K4me1 = readTable("MEF-P1-H3K4me1.tRNA_repeat/tss.datatable")

P4_H3K4me1 = readTable("MEF-P4-H3K4me1.tRNA_repeat/tss.datatable")

P7_H3K4me1 = readTable("MEF-P7-H3K4me1.tRNA_repeat/tss.datatable")

P1_H3K27Ac = readTable("MEF-P1-H3K27Ac.tRNA_repeat/tss.datatable")

P4_H3K27Ac = readTable("MEF-P4-H3K27Ac.tRNA_repeat/tss.datatable")

P7_H3K27Ac = readTable("MEF-P7-H3K27Ac.tRNA_repeat/tss.datatable")

mat = cbind(P1_H3K4me1, P4_H3K4me1, P7_H3K4me1)
mat = cbind(P1_H3K27Ac, P4_H3K27Ac, P7_H3K27Ac)

