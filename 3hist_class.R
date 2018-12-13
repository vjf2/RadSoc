options(stringsAsFactors = FALSE)

library(colorspace)

# rm(list=setdiff(ls(), "rels"))

estimator<-"MEst"

v5<-rels #or rels for simulation
v5$IBD<-v5[,estimator]
# v5$expected<-v5$biparental

catcol<-rainbow_hcl(12, alpha=0.5)[c(1,3,5,8,10)]
borcol<-rainbow_hcl(12, alpha=0.9)[c(1,3,5,8,10)]
lcol<-rainbow_hcl(12)[c(1,3,5,8,10)]

breks<-seq(min(v5$IBD), 0.6, 0.01)

#add negative values to breaks
#implement backward search

windows()
# pdf(file="rcthresholds.pdf")
par(mar=c(5,4.1,0,0))
{hist(v5$IBD[which(v5$expected==0)], probability = FALSE, 
      col=catcol[1],
      yaxt="n",
      xlim=c(min(v5$IBD),0.55), 
      # ylim=c(0,65),
      xlab="Genomic Relatedness", 
      main=" ",
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
}

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

legend(x=0.3, y=60, legend=c("0", "0.0625", "0.125", "0.25", "0.5"), pch=22, pt.cex=2, title="Pedigree\nRelatedness", col=borcol, pt.bg = catcol, bty="n")

#Use brute force to find appropriate thresholds and divide 

