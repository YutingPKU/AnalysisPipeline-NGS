
library(VennDiagram)
res.cp <- fread("results/DEgenes/intra-species/CP.DEgene.human.vs.macaque.log2FC2.FDR001.txt", header = T)
res.de.cp <- res.cp[which(abs(res.cp$log2FoldChange) >2 & res.cp$padj < 0.001),]
de.cp <- res.cp$gene[which(abs(res.cp$log2FoldChange) >2 & res.cp$padj < 0.001)]
res.gz <- fread("results/DEgenes/intra-species/GZ.DEgene.human.vs.macaque.log2FC2.FDR001.txt", header = T)
res.de.gz <- res.gz[which(abs(res.gz$log2FoldChange) >2 & res.gz$padj < 0.001),]
de.gz <- res.gz$gene[which(abs(res.gz$log2FoldChange) >2 & res.gz$padj < 0.001)]

pb.mat <- read.xlsx("public/aat8077_Tables-S7-S23.xlsx", sheet = 5, startRow = 3, colNames = T)
pre.mat <- pb.mat[2:8952,1:17]
colnames(pre.mat) <- c("geneid", pb.mat[1,2:17])
index.pre  <- apply(pre.mat[,2:17], 1, function(vec){
  if(all(vec == "FALSE")){
    return(0)
  }else{
    return(1)
  }
})
pre.mat <- pre.mat[which(index.pre == 1), ]

de.pre <- substr(pre.mat$geneid[which(index.pre == 1)],1,15)
de.pb<- substr(pre.mat$geneid,1,15)
length(intersect(de.gz, de.cp))
length(intersect(de.pre, de.cp))
length(intersect(de.pre, de.gz))



########## #regions which deg overlap
commtype <- function(gid){
  index <- grep(gid, pre.mat$geneid)
  if(length(index) == 0){
    return(0)
  }else{
    return(length(which(pre.mat[index,2:17] == "TRUE")))
  }
}
comm.type.cp <- lapply(de.cp, commtype)
comm.type.cp <- unlist(comm.type.cp)
table(comm.type.cp)
comm.type.gz <- lapply(de.gz, commtype)
comm.type.gz <- unlist(comm.type.gz)
table(comm.type.gz)

q1.cp <- which(comm.type.cp ==0)
q2.cp <- which(comm.type.cp > 0 & comm.type.cp <= 4 )
q3.cp <- which(comm.type.cp > 4 & comm.type.cp <= 8 )
q4.cp <- which(comm.type.cp > 9 & comm.type.cp <= 12 )
q5.cp <- which(comm.type.cp > 12 & comm.type.cp <= 16 )

boxplot(-log10(res.de.cp$padj)[q1.cp], -log10(res.de.cp$padj)[q2.cp],-log10(res.de.cp$padj)[q3.cp],
        -log10(res.de.cp$padj)[q4.cp], -log10(res.de.cp$padj)[q5.cp], ylim = c(0,300))
mat <- cbind(type = c(rep("0",length(q1.cp)), rep("1-4", length(q2.cp)), 
                      rep("5-8", length(q3.cp)), rep("9-12", length(q4.cp)), rep("13-16", length(q5.cp))), 
             value = c(-log10(res.de.cp$padj)[q1.cp], -log10(res.de.cp$padj)[q2.cp],-log10(res.de.cp$padj)[q3.cp],
                       -log10(res.de.cp$padj)[q4.cp], -log10(res.de.cp$padj)[q5.cp]))
q1.gz <- which(comm.type.gz ==0)
q2.gz <- which(comm.type.gz > 0 & comm.type.gz <= 4 )
q3.gz <- which(comm.type.gz > 4 & comm.type.gz <= 8 )
q4.gz <- which(comm.type.gz > 9 & comm.type.gz <= 12 )
q5.gz <- which(comm.type.gz > 12 & comm.type.gz <= 16 )

boxplot(-log10(res.de.cp$padj)[q1.cp], -log10(res.de.cp$padj)[q2.cp],-log10(res.de.cp$padj)[q3.cp],
        -log10(res.de.cp$padj)[q4.cp], -log10(res.de.cp$padj)[q5.cp], ylim = c(0,300))
mat <- cbind(type = c(rep("0",length(q1.gz)), rep("1-4", length(q2.gz)), 
                      rep("5-8", length(q3.gz)), rep("9-12", length(q4.gz)), rep("13-16", length(q5.gz))), 
             value = c(-log10(res.de.gz$padj)[q1.gz], -log10(res.de.gz$padj)[q2.gz],-log10(res.de.gz$padj)[q3.gz],
                       -log10(res.de.gz$padj)[q4.gz], -log10(res.de.gz$padj)[q5.gz]))
mat <- data.frame(mat)
mat$value <- as.numeric(as.character(mat$value))
mat$type <- factor(mat$type, levels = c("0","1-4","5-8","9-12","13-16"))
linesize = 0.4
fontsize = 2
pdf("results/DEgenes/intra-species/compare.with.public.prenatal.boxplot.GZ.DEgene.withdiffernet.overlapration.pdf", width = 2, height = 2)
ggplot(mat, aes(type, value, fill=type)) + geom_boxplot(outlier.shape = NA, notch = T)+stat_boxplot(geom ='errorbar')+
  ylab("-log10 P.adj") +
  xlab("overlap ratio in 16 regions")+
  scale_fill_manual(values =cbPalette <- c("azure3", rgb(131,133,187, maxColorValue = 255),"cornflowerblue", rgb(242,189,62, maxColorValue = 255), 
                                           "chartreuse3"))+
  scale_y_continuous(breaks = seq(0,150,50),  
                     limits = c(0, 150))+
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background = element_blank(),panel.border = element_blank())+
  theme(
        axis.line.y = element_line(colour = "black", size = linesize),axis.ticks.x=element_blank(),
        axis.line.x = element_line(colour = "black", size = linesize),
        axis.ticks.y = element_line(size = 1.5),legend.text=element_blank(),
        axis.ticks.length = unit(.05, "cm"),
        plot.title = element_text(size=fontsize, hjust = 0.5))
ggplot(mat, aes(type, value, fill=type)) + geom_boxplot(outlier.shape = NA, notch = T)+stat_boxplot(geom ='errorbar')+
  ylab("-log10 P.adj") +
  xlab("overlap ratio in 16 regions")+
  scale_fill_manual(values =cbPalette <- c("azure3", rgb(131,133,187, maxColorValue = 255),"cornflowerblue", rgb(242,189,62, maxColorValue = 255), 
                                           "chartreuse3"))+
  scale_y_continuous(breaks = seq(0,150,50), name = "", labels = c("","","",""), 
                     limits = c(0, 150))+
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background = element_blank(),panel.border = element_blank())+
  theme(axis.text.x = element_blank(),axis.text.y = element_blank(),  
        axis.title.x = element_blank(),axis.title.y = element_blank(),
        axis.line.y = element_line(colour = "black", size = linesize),axis.ticks.x=element_blank(),
        axis.line.x = element_line(colour = "black", size = linesize),
        axis.ticks.y = element_line(size = linesize),legend.text=element_blank(),
        axis.ticks.length = unit(.05, "cm"),
        plot.title = element_text(size=fontsize, hjust = 0.5), 
        legend.position = "none")
dev.off()


########### overlap genes per regions
de.region <- apply(pre.mat[,2:17], 2, function(vec){
  loci <- which(unlist(vec) == "TRUE")
  return(substr(pre.mat$geneid,1,15)[loci])
})

comm.cp <- lapply(de.region, function(vec){
  return(length(intersect(vec, de.cp)))
})
comm.gz <- lapply(de.region, function(vec){
  return(length(intersect(vec, de.gz)))
})

all=unlist(lapply(de.region, length))
comm.ratio <- data.frame(cbind(regions = colnames(pre.mat)[-1],
                               cp=unlist(comm.cp)/all, gz = unlist(comm.gz)/all))
names(comm.ratio) <- c("regions","cp","gz")
comm.ratio$cp <- as.numeric(as.character(comm.ratio$cp))
comm.ratio$gz <- as.numeric(as.character(comm.ratio$gz))
#comm.ratio$rcp <- comm.ratio$cp/comm.ratio$public
#comm.ratio$rgz <- comm.ratio$gz/comm.ratio$public
#rownames(comm.ratio) <- colnames(pre.mat)[-1]



mel.mat <- melt(comm.ratio, id = "regions")
pdf("results/DEgenes/intra-species/compare.with.public.prenatal.barplot.pdf", width = 2, height = 1)
ggplot(mel.mat, aes(regions, value, fill=variable)) + 
  geom_bar(stat="identity", position=position_dodge(), colour=NA, width = 0.6)+
  ylab("#overlap/#public") +
  ggtitle("cp:5621 gz:5101 pb:5029")+
  ylim(0,0.4)+
  xlab("")+
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background = element_blank(),panel.border = element_blank())+
  theme(axis.text.x = element_text(size=2),axis.text.y = element_text(size=2),  
        axis.title.x = element_text(size=2),axis.title.y = element_text(size=2),
        axis.line.y = element_line(colour = "black", size = .5),axis.ticks.x=element_blank(),
        axis.line.x = element_line(colour = "black", size = .5),
        axis.ticks.y = element_line(size = .5),legend.text=element_text(size=2),
        axis.ticks.length = unit(.05, "cm"),
        axis.ticks = element_line(size = 0.5),
        plot.title = element_text(size=2, hjust = 0.5))+
  scale_fill_manual(values = c(cp = rgb(131,133,187, maxColorValue = 255), gz = rgb(242,189,62, maxColorValue = 255)))

dev.off()

plotven.withoulabel.triple <- function(Vdemo){
  w <- compute.Venn(Vdemo)
  gp <- VennThemes(w)
  gp$FaceText$`001`$fontsize <- 0
  gp$FaceText$`011`$fontsize <- 0
  gp$FaceText$`010`$fontsize <- 0
  gp$FaceText$`001`$fontsize <- 0
  gp$FaceText$`101`$fontsize <- 0
  gp$FaceText$`100`$fontsize <- 0
  gp$FaceText$`111`$fontsize <- 0
  gp$FaceText$`110`$fontsize <- 0
  gp$SetText$Set3$fontsize <- 0
  gp$SetText$Set1$fontsize <- 0
  gp$SetText$Set2$fontsize <- 0
  gp$Set$Set1$lty <- 0
  gp$Set$Set2$lty <- 0
  gp$Set$Set3$lty <- 0
  grid.newpage() 
  plot(w, gp=gp, show=list(Universe=FALSE))
}
plotven.withoulabel.triple.col <- function(Vdemo, col011, col010, col001, col101, col100, col111, col110){
  w <- compute.Venn(Vdemo)
  gp <- VennThemes(w)
  gp$FaceText$`001`$fontsize <- 0
  gp$FaceText$`011`$fontsize <- 0
  gp$FaceText$`010`$fontsize <- 0
  gp$FaceText$`001`$fontsize <- 0
  gp$FaceText$`101`$fontsize <- 0
  gp$FaceText$`100`$fontsize <- 0
  gp$FaceText$`111`$fontsize <- 0
  gp$FaceText$`110`$fontsize <- 0
  gp$SetText$Set3$fontsize <- 0
  gp$SetText$Set1$fontsize <- 0
  gp$SetText$Set2$fontsize <- 0
  gp$Set$Set1$lty <- 0
  gp$Set$Set2$lty <- 0
  gp$Set$Set3$lty <- 0
  gp$Face$`011`$fill <- col011
  gp$Face$`010`$fill <- col010
  gp$Face$`001`$fill <- col001
  gp$Face$`101`$fill <- col101
  gp$Face$`100`$fill <- col100
  gp$Face$`111`$fill <- col111
  gp$Face$`110`$fill <- col110
  grid.newpage() 
  plot(w, gp=gp, show=list(Universe=FALSE))
}
draw.triple.venn(area1 = 5621, area2 = 5101, area3 = 5029, n12 = 3695, n13 = 977, n23 = 895, n123 = 711, 
                 category = c("cp", "gz", "public-prenatal"), euler.d = T, scaled = T, 
                 fill = c("cornflowerblue", "cyan4", "lightpink2"),
                 cex = 3,cat.cex = 2,cat.col = c( "red", "green", "blue"))
VenTest <- Venn(SetNames = c("cp", "gz", "public-prenatal"), 
                Weight = c(`100` = 1660, `010`=1222, `001`=3868, `110`=2984,
                           `101`=266, `011`=184, `111`=711))
pdf("results/DEgenes/intra-species/compare.with.public.prenatal.vennplot.pdf", width = 2, height = 2)
par(mai = c(0.1,0.1,0.1,0.1))
plotven.withoulabel.triple.col(VenTest, col100 = rgb(131,133,187, maxColorValue = 255),
                               col010 = rgb(242,189,62, maxColorValue = 255),
                               col001 = rgb(141,210,199, maxColorValue = 255),
                               col110 = rgb(128,177,211, maxColorValue = 255), 
                               col101 = rgb(251,128,114, maxColorValue = 255),
                               col011 = rgb(253,180,98, maxColorValue = 255), 
                               col111 = rgb(179,222,105, maxColorValue = 255))
dev.off()
pdf("results/DEgenes/intra-species/compare.with.public.prenatal.vennplot.withlabel.pdf", width = 3, height = 3)
plot(VenTest)
dev.off()