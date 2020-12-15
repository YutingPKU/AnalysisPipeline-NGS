
all <- readRDS("data/monkey-exonbygene-summaried-expriment-rmd.rds")
sampleData <- colData(all)
mat <- read.xlsx("results/monkey-RNASeq-cqnnormalized-FPKM-log2-rmd-exonbygene.xlsx", colNames = T, rowNames = T)
mat <- 2^mat
load("data/Macaca_mulatta.Mmul_8.0.1.89.chr.gene.GRanges.rds")


#gname <- c("NES","SOX2","NOTCH1","HES1","HES3","PAX6","FOXG1")
gname <- c("NEUROD1","TBR1","MAP2","FOXA2","FEZF2","NEUROD2","STMN1","SYP","DLG4","HRNBP3")
loci <- match(gname, my_gr$gn)
loci <- match(my_gr$id[loci], rownames(mat))

caseplot <- function(gname){
  loci <- match(gname, my_gr$gn)
  loci <- match(my_gr$id[loci], rownames(mat))
  case.mat <- data.frame(region = colData(all)[,2], value = unlist(mat[loci,]))
  data_summary <- function(data, varname, groupnames){
    require(plyr)
    summary_func <- function(x, col){
      c(mean = mean(x[[col]], na.rm=TRUE),
        sd = sd(x[[col]], na.rm=TRUE))
    }
    data_sum<-ddply(data, groupnames, .fun=summary_func,
                    varname)
    data_sum <- rename(data_sum, c("mean" = varname))
    return(data_sum)
  }
  df2 <- data_summary(case.mat, varname="value", 
                      groupnames=c("region"))
  # Convert dose to a factor variable
  df2$region=as.factor(df2$region)
  
  linesize = 0.4
  fontsize = 6
  ggplot(df2, aes(x=region, y=value, fill=region)) + 
    geom_bar(stat="identity", 
             position=position_dodge(), width = .5) +
    geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.2,
                  position=position_dodge(.9)) +
    labs(title=gname,  y = "FPKM")+
    theme_classic() +
    scale_fill_manual(values=c("CP" = rgb(131,133,187, maxColorValue = 255),"GZ" = rgb(242,189,62, maxColorValue = 255)))+
    #scale_fill_manual(values=c("CP" = "cornflowerblue","GZ" = "chartreuse3"))+
    theme_bw() +
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background = element_blank(),panel.border = element_blank())+
    theme(axis.text.x = element_text(size=fontsize),axis.text.y = element_text(size=fontsize),  
          axis.title.x = element_text(size=fontsize),axis.title.y = element_text(size=fontsize),
          axis.line.y = element_line(colour = "black", size = linesize),axis.ticks.x=element_blank(),
          axis.line.x = element_line(colour = "black", size = linesize),
          axis.ticks.y = element_line(size = linesize),legend.text=element_text(size=fontsize),
          axis.ticks.length = unit(.05, "cm"),
          plot.title = element_text(size=fontsize, hjust = 0.5))
  ggplot(df2, aes(x=region, y=value, fill=region)) + 
    geom_bar(stat="identity", 
             position=position_dodge(), width = .5) +
    geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.2,
                  position=position_dodge(.9)) +
    labs(title=gname,  y = "FPKM")+
    theme_classic() +
    scale_fill_manual(values=c("CP" = rgb(131,133,187, maxColorValue = 255),"GZ" = rgb(242,189,62, maxColorValue = 255)))+
    theme_bw() +
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background = element_blank(),panel.border = element_blank())+
    theme(legend.position = "none",
          axis.text.x = element_text(size=fontsize),axis.text.y = element_text(size=fontsize),  
          axis.title.x = element_text(size=fontsize),axis.title.y = element_text(size=fontsize),
          axis.line.y = element_line(colour = "black", size = linesize),axis.ticks.x=element_blank(),
          axis.line.x = element_line(colour = "black", size = linesize),
          axis.ticks.y = element_line(size = linesize),legend.text=element_text(size=fontsize),
          axis.ticks.length = unit(.05, "cm"),
          plot.title = element_text(size=fontsize, hjust = 0.5))
  
 
  
  
}

pdf("results/markerGeneBarplot/CP.marker.pdf", width = 1.5, height = 1.5)
lapply(gname, caseplot)
dev.off()
 

