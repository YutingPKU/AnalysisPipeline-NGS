
################# plot heatmap for hic data
#########################heatmap
#case_name <- "5534N:Raw:genome"
#输入矩阵路径

plotheatmap <- function(mat, case_name,resolution){
  raw_matrix <- as.matrix(mat) #raw_matrix就是Hi-C最终的原始矩阵
  
  sum_col <- colSums(raw_matrix,na.rm = T)
  real_col_num <- length(sum_col[sum_col!=0]) #计算列加和不为0的数量，即非N区域
  real_col_num
  
  #par(mar=c(0,5,2,5))
  
  row_num <- dim(raw_matrix)[1]
  col_num <- dim(raw_matrix)[2]
  start = 0
  end = row_num
  
  #生成坐标vector
  x1 <- rep(c(0:(col_num-1)),row_num)
  x2 <- x1 + 1
  y1 <- rep(c(0:(-col_num+1)),each=row_num)
  y2 <- y1 -1 
  par(mar = c(2,8,4,2))
  #plot(x=0:col_num,y=-row_num:0,type="n",frame.plot = F,cex.axis=1,xaxt = "n",
  #cex.lab=1.5,ylab = "Contact matrix",xlab = "",main = case_name,cex.main = 2)
  plot(x=0:col_num,y=-row_num:0,type="n",frame.plot = F,cex.axis=2,axes = FALSE,
       cex.lab=2.5,ylab = "",xlab = "",main = case_name,cex.main = 2)
  nn = 25 #调整画图时标注染色体位置的线，分为几段
  #调整画图时标注的数字：染色体位置
  #axis(side = 1,cex = 2, at = seq(0,col_num,by = 2*nn),labels = unlist(lapply(prettyNum((seq(start,end,by = 2*nn))*resolution/1000000,big.mark=","),paste0,"Mb")))
  
  #生成颜色vector
  region <- quantile(raw_matrix,probs = c(0,0.9,1),na.rm = T)
  raw_vector <- as.vector(raw_matrix)
  color_vector <- raw_vector / floor(region[2])
  color_vector[color_vector > 1 ] <- 1
  red_vector <- rep(1,length(raw_vector))
  green_vector <- 1 - color_vector
  blue_vector <- 1 - color_vector
  
  rect(x1,y1,x2,y2,col = rgb(red_vector,green_vector,blue_vector),border = NA)
}


#################### plot heatmap and AB compartments
pdf("results/CP-chr1-ABcompartment.pdf", width = 3.5, height = 3.5)
layout(matrix(c(1:4),4,1,byrow = TRUE),heights = c(0.1,6,0.5,0.5), widths = rep(6, times = 4))
layout.show()

case_name <- "CP chr1 Resolution 500kb"
mat = fread("../matrix/merged-CP/iced/500000/CP_500000_iced_chr1_dense.matrix", header = F, sep = "\t")
#mat = mat[,-c(1:3)]
mat = as.matrix(mat)
mat[mat >= quantile(mat, 0.95)] = quantile(mat, 0.95)
nn = 25 
col_num = end = dim(mat)[1]
start = 0
resolution = 500000

par(mar = c(2,4,2,6))
image(x=0:(nrow(mat)-1),y=0:(ncol(mat)-1), mat, col = colorRampPalette(c("navy","red3"))(100),  axes = F, xlab = "", ylab = "")
#mtext(text=unlist(lapply(prettyNum((seq(start,end,by = 2*nn))*resolution/1000000,big.mark=","),paste0,"Mb")), side=1, line=0.3, at = seq(0,col_num,by = 2*nn), las=1, cex=0.8)
#axis(side = 1,cex = 2, at = seq(0,col_num,by = 2*nn),labels = unlist(lapply(prettyNum((seq(start,end,by = 2*nn))*resolution/1000000,big.mark=","),paste0,"Mb")))
#image.plot( mat, legend.only = T, legend.shrink = 0.2, col = colorRampPalette(c("navy","red3"))(100) )

plotheatmap(mat, case_name, 500000)

comp = fread("compartments/MEF_HiTC_ABcomparments_500kb_wholegenome.txt", header = T, sep = "\t")
comp = comp[which(comp$seqnames == "chr1"),c(3:5,8:9)]
comp$score[is.na(comp$score)] = 0
y = as.numeric(comp$score)
x<- seq(0:(length(y)-1))
d = as.numeric(comp$genedens)

par(mar = c(0,8,0,2))
plot(x,y,type = 'n',xaxt = "n",frame.plot = F,font.lab =4,
     ylab = "PCA1vec",xlab = "")

segments(x,rep(0,10),x,y,ljoin = 0,lwd = 1.2,col =ifelse(y>0,"red", "blue"))
#segments(x,rep(0,10),x,y,ljoin = 0,lwd = 1,col =ifelse(y>0,rgb(219,110,49, max = 255), 
                                                       #rgb(72,111,167, max = 255)) )
par(mar = c(0,8,0,2))
plot(1:length(y),d,type = 'n',xaxt = "n",frame.plot = F,font.lab =4,ylab = "Gene density",xlab = "")

segments(x,rep(0,10),x,d,ljoin = 0,lwd = 1.2,col = 'black')

dev.off()





