library(fields)



matrix <- fread("data/GZ_chr1_500kb_KRnorm_OE_matrix.txt", header = F)


sub.rmat <- log2(matrix+0.0001)
#sub.rmat[is.na(sub.rmat)] <- 0
nrow <- nrow(sub.rmat)
cutvalue.pos <- 3.5
cutvalue.neg <- -3.5
sub.rmat[which(sub.rmat > cutvalue.pos)] = cutvalue.pos
sub.rmat[which(sub.rmat < cutvalue.neg)] = cutvalue.neg


breakslist <- seq(cutvalue.neg, cutvalue.pos, 0.01)
cbreaks <- length(breakslist) - 1
myColor <- colorRampPalette(c(rgb(0,0,1), rgb(.1,.1,1),rgb(.2,.2,1), rgb(.3,.3,1),rgb(.4, .4,1)  ,rgb(.8,.8,1),rgb(1,1,1),
                              rgb(1,.8,.8), rgb(1, .4,.4),rgb(1,.3,.3) ,rgb(1, .2,.2), rgb(1, .1,.1), rgb(1, 0,0)))(cbreaks)

myColor <- colorRampPalette(c(rgb(0,0, 255, maxColorValue = 255), rgb(54,41, 255, maxColorValue = 255),
                              rgb(104,99,255, maxColorValue = 255), rgb(180,180, 255, maxColorValue = 255),
                              rgb(255, 255, 255, maxColorValue = 255),
                              rgb(255, 180,180, maxColorValue = 255), rgb(255,95,95, maxColorValue = 255),
                              rgb(255,33,22, maxColorValue = 255),
                              rgb(255, 0,0, maxColorValue = 255)))(cbreaks)

sub.rmat <- as.matrix(sub.rmat)
pdf("results/GZ-chr1-ABcompartment-KR-OEmatrix.pdf", width = 2, height = 2)
par(mar = c(1,1,1,1))
image(x=0:(nrow(sub.rmat)-1),y= -(ncol(sub.rmat)-1):0, sub.rmat[,nrow:1],  col = myColor,
      axes = F, xlab = "", ylab = "", breaks = breakslist)

#axis(2, at= seq(-nrow,0, 50), labels = seq(0, nrow, 50) , lwd = 0.2, lwd.ticks = 0.2)


image(x=0:(nrow(sub.rmat)-1),y=-(ncol(sub.rmat)-1):0, sub.rmat[,nrow:1],  col = myColor,
      axes = F, xlab = "", ylab = "", breaks = breakslist)
#axis(2, at= seq(-nrow,0, 50), labels = rep("", length(seq(0, nrow, 50)) ), lwd = 0.2, lwd.ticks = 0.2)

image.plot(sub.rmat, legend.only = F, legend.shrink = 0.2, col = myColor, add = F,
           breaks = breakslist, legend.cex = .1)
dev.off()
########## change color palette
par(mar = c(1,1,1,1))
image(x=0:(nrow(sub.rmat)-1),y= -(ncol(sub.rmat)-1):0, sub.rmat[,nrow:1],  col = colorRampPalette(colorRamps::blue2red(10))(cbreaks),
      axes = F, xlab = "", ylab = "", breaks = breakslist)
axis(2, at= -seq(0, nrow, 50), labels = seq(0, nrow, 50) , lwd = 0.2, lwd.ticks = 0.2)
addloop(start, end, resolution, loop.start = loop$x2 , loop.end = loop$y2 , col = "black", linesize = 0.3)

image(x=0:(nrow(sub.rmat)-1),y=-(ncol(sub.rmat)-1):0, sub.rmat[,nrow:1],  col = colorRampPalette(colorRamps::blue2red(10))(cbreaks),
      axes = F, xlab = "", ylab = "", breaks = breakslist)
axis(2, at= -seq(0, nrow, 50), labels = rep("", length(seq(0, nrow, 50)) ), lwd = 0.2, lwd.ticks = 0.2)


image.plot(sub.rmat, legend.only = F, legend.shrink = 0.2, col = colorRampPalette(colorRamps::blue2red(10))(cbreaks), add = F,
           breaks = breakslist, legend.cex = .1)
dev.off()


