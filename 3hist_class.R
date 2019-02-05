options(stringsAsFactors = FALSE)

library(colorspace)

v5<-read.csv("rels.csv")

estimator<-"MEst"

v5$IBD<-v5[,estimator]

catcol<-rainbow_hcl(12, alpha=0.5)[c(1,3,5,8,10)]
borcol<-rainbow_hcl(12, alpha=0.9)[c(1,3,5,8,10)]
lcol<-rainbow_hcl(12)[c(1,3,5,8,10)]

breks<-seq(min(v5$IBD), 0.6, 0.01)

# windows()
pdf(file="Figures/Figure_3.pdf", width = 4.40945, height= 4.40945 )
par(mar=c(3.5,3.4,0,0))
hist(v5$IBD[which(v5$expected==0)], probability = FALSE, 
      col=catcol[1],
      yaxt="n",
      xlim=c(min(v5$IBD),0.6), 
      # ylim=c(0,65),
      xlab="", 
      ylab="",
      main="",
      border=borcol[1], 
      breaks=breks)
  hist(v5$IBD[which(v5$expected==0.25)], probability = FALSE, 
       col=catcol[4], add=TRUE, border=borcol[4], breaks=breks)
  hist(v5$IBD[which(v5$expected==0.5 & v5$RelCat=="PO")], probability = FALSE, 
       col=catcol[5], add=TRUE, border=borcol[5], breaks=breks)
  hist(v5$IBD[which(v5$expected==0.125)], probability = FALSE, 
       col=catcol[3], add=TRUE, border=borcol[3], breaks=breks)
  hist(v5$IBD[which(v5$expected==0.0625)], probability = FALSE, 
       col=catcol[2], add=TRUE, border=borcol[2], breaks=breks)
  axis(2, las=1)

  segments(0, 0, 0.6, 0)

mtext("Frequency", side = 2, line = 2.5)
mtext("Genomic Relatedness", side = 1, line = 2.5)


w<-v5$IBD[which(v5$expected==0.0)]
x<-v5$IBD[which(v5$expected==0.5 & v5$RelCat=="PO")]
y<-v5$IBD[which(v5$expected==0.25)]
z<-v5$IBD[which(v5$expected==0.125)]
a<-v5$IBD[which(v5$expected==0.0625)]

w1<-density(w, adjust=1)
x1<-density(x, adjust=2.4)
y1<-density(y, adjust=0.8)
z1<-density(z, adjust=0.8)
a1<-density(a, adjust=0.8)

lines(w1$x, w1$y, lwd=2, col=lcol[1])
lines(x1$x, x1$y, lwd=2, col=lcol[5])
lines(y1$x, y1$y, lwd=2, col=lcol[4])
lines(z1$x, z1$y, lwd=2, col=lcol[3])
lines(a1$x, a1$y, lwd=2, col=lcol[2])

lines(w1$x, w1$y, lty=2)
lines(x1$x, x1$y, lty=2)
lines(y1$x, y1$y, lty=2)
lines(z1$x, z1$y, lty=2)
lines(a1$x, a1$y, lty=2)

legtext <- c("0", "0.0625", "0.125", "0.25", "0.5")
xcoords <- c(0, 0.0625, 0.125, 0.25, 0.5)
secondvector <- (1:length(legtext))-1
textwidths <- xcoords/secondvector # this works for all but the first element
textwidths[1] <- 0 


legend(x=0, y=60, legend=c("0", "0.0625", "0.125", "0.25", "0.5"), pch=22, pt.cex=1, title="", col=borcol, pt.bg = catcol, bty="n", horiz=TRUE, x.intersp = 0.5, 
       text.width = c(0.02, 0.005, 0.035, 0.045, 0.09), xpd=TRUE, cex=0.8)
       
text(x=0.3, y=58,"Pedigree Relatedness", cex=0.9)

b1<-0.0362
b2<-0.0935
b3<-(max(v5$IBD[which(v5$expected==0.125)])+min(v5$IBD[which(v5$expected==0.25)]))/2
b4<-(max(v5$IBD[which(v5$expected==0.25)])+min(v5$IBD[which(v5$expected==0.5)]))/2

segments(x0=b1, y0=0, x1=b1, y1=50, lty=5, lwd=1)
segments(x0=b2, y0=0, x1=b2, y1=50, lty=5, lwd=1)
segments(x0=b3, y0=0, x1=b3, y1=50, lty=5, lwd=1)
segments(x0=b4, y0=0, x1=b4, y1=50, lty=5, lwd=1)


dev.off()


#With false negatives and true positives for slides
#Use brute force to find appropriate thresholds and divide 

#For plotting
# text(x = 0.021, y=40, "97%", cex=0.8)
# text(x = 0.065, y=40, "85%", cex=0.8)
# text(x = 0.13, y=40, "96%", cex=0.8)
# text(x = 0.25, y=40, "100%", cex=0.8)
# text(x = 0.45, y=40, "100%", cex=0.8)
# 
# text(x = 0.021, y=37, "5%", col="red", cex=0.8) #false positive
# text(x = 0.065, y=37, "5%", col="red", cex=0.8) #false positive
# text(x = 0.13, y=37, "5%", col="red", cex=0.8) #false positive
# text(x = 0.25, y=37, "0%", col="red", cex=0.8) #false positive
# text(x = 0.45, y=37, "0%", col="red", cex=0.8) #false positive
