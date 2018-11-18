options(stringsAsFactors = FALSE)

library(colorspace)

# rm(list=setdiff(ls(), "rels"))

estimator<-"MEst"

v5<-reals #or rels for simulation
v5$IBD<-v5[,estimator]
# v5$expected<-v5$biparental

catcol<-rainbow_hcl(12, alpha=0.5)[c(1,5,8,10,12)]
borcol<-rainbow_hcl(12, alpha=0.9)[c(1,5,8,10,12)]
lcol<-rainbow_hcl(12)[c(1,5,8,10,12)]

breks<-seq(min(v5$IBD), 0.6, 0.01)

#add negative values to breaks
#implement backward search

windows()
# pdf(file="rcthresholds.pdf")
{hist(v5$IBD[which(v5$expected==0.25)], probability = FALSE, 
      col=catcol[3],
      yaxt="n",
      xlim=c(min(v5$IBD),0.6), 
      ylim=c(0,65),
      xlab="Genomic Relatedness", 
      main=" ",
      border=borcol[3])
  hist(v5$IBD[which(v5$expected==0)], probability = FALSE, 
       col=catcol[1], add=TRUE, border=borcol[1], breaks=breks)
  hist(v5$IBD[which(v5$expected==0.5 & v5$RelCat=="PO")], probability = FALSE, 
       col=catcol[4], add=TRUE, border=borcol[4], breaks=breks)
  hist(v5$IBD[which(v5$expected==0.125)], probability = FALSE, 
       col=catcol[5], add=TRUE, border=borcol[5], breaks=breks)
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
x1<-density(x)
y1<-density(y)
z1<-density(z)
a1<-density(a)

lines(w1$x, w1$y, lwd=2, col=lcol[1])
lines(x1$x, x1$y, lwd=2, col=lcol[4])
lines(y1$x, y1$y, lwd=2, col=lcol[3])
lines(z1$x, z1$y, lwd=2, col=lcol[5])
lines(a1$x, a1$y, lwd=2, col=lcol[2])

lines(w1$x, w1$y, lty=2)
lines(x1$x, x1$y, lty=2)
lines(y1$x, y1$y, lty=2)
lines(z1$x, z1$y, lty=2)
lines(a1$x, a1$y, lty=2)

legend(x=0.25, y=50, legend=c("0", "0.0625", "0.125", "0.25", "0.5"), pch=22, pt.cex=2, title="Pedigree\nRelatedness", col=borcol, pt.bg = catcol, bty="n")

#try to search for a threshold 

v6<-v5[order(v5[,estimator]),]

thresh<-seq(0.0005, 0.6, 0.0005)

correct<-list()
incorrect<-list()

for (i in 1:length(thresh)) {
  
  v6$classification<-cut(v6$IBD, c(min(v6$IBD), thresh[i]), include.lowest = TRUE)
  
  tabl<-as.data.frame(table(v6$classification, v6$RelCat))
  
  correct[[i]]<-tabl[tabl$Var2=="UR", "Freq"]
  incorrect[[i]]<-sum(tabl[!tabl$Var2=="UR", "Freq"])
  
}

res<-as.data.frame(cbind(thresh, unlist(correct), unlist(incorrect)))

res$fp<-res[,3]/(res[,2]+res[,3])

co<-res[which(res$fp<=0.05),]

co_res<-res[which(res$thresh==min(co$thresh[co$V2==max(co$V2)])),]

co<-co_res$thresh

#between hc and ur
correct<-list()
incorrect<-list()

thresh<-thresh[which(thresh>co)]

for (i in 1:length(thresh)) {
  
  v6$classification<-cut(v6$IBD, c(co, thresh[i]), include.lowest = TRUE)
  
  tabl<-as.data.frame(table(v6$classification, v6$RelCat))
  
  correct[[i]]<-tabl[tabl$Var2=="HC", "Freq"]
  incorrect[[i]]<-sum(tabl[!tabl$Var2=="HC", "Freq"])
  
}

res<-as.data.frame(cbind(thresh, unlist(correct), unlist(incorrect)))

res$fp<-res[,3]/(res[,2]+res[,3])

co1<-res[which(res$fp<=0.05),]

co1_res<-res[which(res$thresh==min(co1$thresh[co1$V2==max(co1$V2)])),]

co1<-co1_res$thresh

#Avuncular
correct<-list()
incorrect<-list()

thresh<-thresh[which(thresh>co1)]

for (i in 1:length(thresh)) {
  
  v6$classification<-cut(v6$IBD, c(co1, thresh[i]), include.lowest = TRUE)
  
  tabl<-as.data.frame(table(v6$classification, v6$RelCat))
  
  correct[[i]]<-tabl[tabl$Var2=="FC", "Freq"]
  incorrect[[i]]<-sum(tabl[!tabl$Var2=="FC", "Freq"])
  
}

res<-as.data.frame(cbind(thresh, unlist(correct), unlist(incorrect)))

res$fp<-res[,3]/(res[,2]+res[,3])

co2<-res[which(res$fp<=0.05),]

co2_res<-res[which(res$thresh==min(co2$thresh[co2$V2==max(co2$V2)])),]

co2<-co2_res$thresh

#just take midpoint for co3 for now

correct<-list()
incorrect<-list()

thresh<-thresh[which(thresh>co2)]

for (i in 1:length(thresh)) {
  
  v6$classification<-cut(v6$IBD, c(co2, thresh[i]), include.lowest = TRUE)
  
  tabl<-as.data.frame(table(v6$classification, v6$RelCat))
  
  correct[[i]]<-tabl[tabl$Var2=="HS", "Freq"]
  incorrect[[i]]<-sum(tabl[!tabl$Var2=="HS", "Freq"])
  
}

res<-as.data.frame(cbind(thresh, unlist(correct), unlist(incorrect)))

res$fp<-res[,3]/(res[,2]+res[,3])

co3<-res[which(res$fp<=0.05),]

co3_res<-res[which(res$thresh==min(co3$thresh[co3$V2==max(co3$V2)])),]

co3<-co3_res$thresh

#for plotting
co3<-(max(v6$IBD[v6$RelCat=="HS"])+min(v6$IBD[v6$RelCat=="FS"]))/2

segments(x0=co, y0=0, x1=co, y1=60, lty=3)
segments(x0=co1, y0=0, x1=co1, y1=60, lty=3)
segments(x0=co2, y0=0, x1=co2, y1=60, lty=3)
segments(x0=co3, y0=0, x1=co3, y1=60, lty=3)

#set up thresholds backwards for other 

final<-do.call("rbind", list(co_res, co1_res, co2_res, co3_res))

final$type<-estimator

write.table(final, "class_res.csv", append = TRUE, sep = ",")

#No Solution for WEst, LREst, REst, and QGEst

#"QGEst" 

#Take the quantiles of the density curves 

