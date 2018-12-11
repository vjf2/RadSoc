#Figures


windows();plot(restable$nind, restable$pvalue, type="n")
text(restable$nind, restable$pvalue, labels=restable$nsnps)
abline(h=0.05, col="red", lty=2)

windows();plot(restable$nind, restable$coef, type="n")
text(restable$nind, restable$coef, labels=restable$nsnps)
abline(h=coef(lqap[[57]])[2], col="red", lty=2)
abline(h=coef(lqap[[57]])[2]+lqap[[56]]$se[1], col="grey", lty=2)
abline(h=coef(lqap[[57]])[2]-lqap[[56]]$se[1], col="grey", lty=2)


windows();plot(restable$nind, restable$rsqrel, type="n")
text(restable$nind, restable$rsqrel, labels=restable$nsnps)
abline(h=restable$rsqrel[nrow(restable)], col="red", lty=2)

#sample size, allele frequency, number of markers

#add order that samples were obtained?

#use best estimates to determine sample size
#change allele freq to determine accuracy
#get snp code to make nicer graph


library(colorspace)
xc<-rainbow_hcl(8)

windows();plot(restable$nind, restable$rsqrel, type="n", ylim=c(0,0.6))

points(rsqrel~nind, data=restable[restable$nsnps==50,], col=xc[1], pch=1)
points(rsqrel~nind, data=restable[restable$nsnps==100,], col=xc[2], pch=2)
points(rsqrel~nind, data=restable[restable$nsnps==200,], col=xc[3], pch=3)
points(rsqrel~nind, data=restable[restable$nsnps==400,], col=xc[4], pch=4)
points(rsqrel~nind, data=restable[restable$nsnps==800,], col=xc[5], pch=5)
points(rsqrel~nind, data=restable[restable$nsnps==1600,], col=xc[6], pch=6)
points(rsqrel~nind, data=restable[restable$nsnps==3200,], col=xc[7], pch=7)
points(rsqrel~nind, data=restable[restable$nsnps==4235,], col=xc[8], pch=8)

lines(rsqrel~nind, data=restable[restable$nsnps==50,], col=xc[1], pch=1)
lines(rsqrel~nind, data=restable[restable$nsnps==100,], col=xc[2], pch=2)
lines(rsqrel~nind, data=restable[restable$nsnps==200,], col=xc[3], pch=3)
lines(rsqrel~nind, data=restable[restable$nsnps==400,], col=xc[4], pch=4)
lines(rsqrel~nind, data=restable[restable$nsnps==800,], col=xc[5], pch=5)
lines(rsqrel~nind, data=restable[restable$nsnps==1600,], col=xc[6], pch=6)
lines(rsqrel~nind, data=restable[restable$nsnps==3200,], col=xc[7], pch=7)
lines(rsqrel~nind, data=restable[restable$nsnps==4235,], col=xc[8], pch=8)


legend("topright", legend=sort(unique(restable$nsnps)), col=xc, pch=1:8, lwd=1,
       title="Number of SNPs", y.intersp = 0.9)

abline(h=restable$rsqrel[nrow(restable)], lty=2)



text(restable$nind, restable$proprel, labels=restable$nsnps)


pdata<-restable[,c(1,2,5)]

pdata<-reshape2::acast(pdata, nsnps~nind)

windows();barplot(pdata, beside=TRUE, col=xc)

windows();boxplot(pdata, beside=TRUE, col=xc) #organize for multiple
