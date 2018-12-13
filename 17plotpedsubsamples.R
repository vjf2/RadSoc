#17 figure for kin classifications

options(stringsAsFactors = FALSE)

setwd("C:/Users/froug/Desktop/Real First Chapter")

library(colorspace)
catcol<-rainbow_hcl(12, alpha=0.5)[c(1,5,8,10,12)]
errcol<-"gray30"
allcol<-c(rbind(catcol, errcol))
borcol<-"darkgrey"

kin_class<-read.csv("kin_classifications.csv", row.names=1)

snps<-sort(unique(kin_class$nsnps))

plotpos<-function(i, b=-1){
cclass<-kin_class[kin_class$nsnps==snps[i],]
barplot(rbind(cclass$tp_prop,-cclass$fp_prop), beside=TRUE, col=allcol,
        ylab="", ylim=c(b,1), yaxt="n", xaxt="n", space=c(0,0))
text(1,-0.2, cclass$nsnps[1], cex=1.3)
abline(h=0.95, lty=3)
abline(h=0.80, lty=2)
}


windows()
#pdf(file="Figures/classification_rates.pdf")
par(mfrow=c(2,4), mar=c(0,0,1,0), oma=c(4,6,1,1), xpd=TRUE)

plotpos(1)
axis(2, at=c(-1, -0.5, 0, 0.5, 1), labels=c(1, 0.5, 0, 0.5, 1),las=1)
mtext("Classification Rate", 2, line=4, adj=-0.5, cex=1.1)
plotpos(2)
plotpos(3)
plotpos(4)
plotpos(5,-0.5)
axis(2, at=c(-1, -0.5, 0, 0.5, 1), labels=c(1, 0.5, 0, 0.5, 1),las=1)

plotpos(6, -0.5)

plotpos(7, -0.5)
plotpos(8, -0.5)
legend(x=-17, y=-0.4, legend=c("0", "0.125", "0.25", "0.5"), pch=22, cex=1.5 ,pt.cex=3, title="Pedigree Relatedness", col=borcol, pt.bg = catcol, bty="n", ncol=4, xpd =NA)

dev.off()


plotpos<-function(i){
  
  cv<-al[[pn[i]]]
  
  cv5<-merge_pairs(cv, v5[,c("ID1.x", "ID2.x", "expected")], 
                   "ID1.x", "ID2.x", all.x=FALSE, all.y=FALSE)
  
  cv5$classification<-ifelse(cv5$pi_hat>=ur[1] & cv5$pi_hat<ur[2], "ur", NA)
  cv5$classification<-ifelse(cv5$pi_hat>=fc[1] & cv5$pi_hat<fc[2], "fc", cv5$classification)
  cv5$classification<-ifelse(cv5$pi_hat>=hs[1] & cv5$pi_hat<hs[2], "hs", cv5$classification)
  cv5$classification<-ifelse(cv5$pi_hat>=po[1] & cv5$pi_hat<po[2], "po", cv5$classification)
  
  nms<-c("urE", "fcE", "hsE", "poE")
  
  tabl<-as.matrix(table(cv5$expected, cv5$classification))
  rownames(tabl)<-nms
  tabl<-tabl[,c("ur", "fc", "hs", "po")]
  
  #rates
  
  tp<-diag(tabl)/rowSums(tabl) #true positive
  fp<-(colSums(tabl)-diag(tabl))/colSums(tabl) #false positive
  fn<-(rowSums(tabl)-diag(tabl))/rowSums(tabl) #false negative
  
  nfn<-fn*-1
  nfp<-fp*-1
  
  barplot(rbind(tp,nfp, nfn), beside=TRUE, col=errcol,
          ylab="", ylim=c(-1,1), yaxt="n", xaxt="n", space=c(0,0))
  text(1,-0.8, cv$snps[pn[i]], cex=1.3)
  abline(h=0.95, lty=3)
  abline(h=0.80, lty=2)
}

