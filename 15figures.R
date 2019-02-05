#Figures

library(colorspace)
xc<-rainbow_hcl(8)

options(stringsAsFactors = FALSE)

load("Z:/RadSoc/restables.RData")

j<-1
q<-lapply(q, function(x) {x$iter<-j; x<-x[-nrow(x),];j<<-j+1;return(x)})

qdf<-do.call("rbind", q)

qdf<-qdf[order(qdf$nind),]

p<-aggregate(pvalue~nind+nsnps, data=qdf, function(x) {sum(x>=0.05)/length(x)})

p<-p[order(p$nind),]

max(p$pvalue)

# windows()

pdf(file="Figures/Figure_5.pdf", width=6.65354, height=4.1122, useDingbats=FALSE)
par(mar=c(0, 0, 0, 0), oma=c(4,4,0.5,0.1))
boxplot(rsqrel~nsnps+nind,  
        col=adjustcolor(xc, alpha.f=0.5)[c(2:8,1)], 
        at=c(1:8, 10:17, 19:26, 28:35, 37:44, 46:53, 55:62),
        data=qdf,
        xaxt="n",
        yaxt="n",
        ylab="",
        xlab="",
        ylim=c(min(qdf$rsqrel),0.72), 
        type="n",
        outcex=0.5, 
        whisklty=1,
        whiskcol="gray30",
        outcol="gray30") 

mtext("Number of Individuals", 1, line=2.5, cex=0.9)
mtext("Coefficient of Partial Determination - Relatedness", 2, line=2.5, cex=0.9)

axis(2, las=1, cex.axis=0.9)
axis(1, at=seq(4,62,9), labels=unique(qdf$nind), cex.axis=0.9)

# legend(55,0.70, legend=sort(unique(qdf$nsnps)), pt.bg=xc[c(2:8,1)], pch=22,
#        title="Number of SNPs", y.intersp = 0.9, bty="n")

legend("topright", legend=sort(unique(qdf$nsnps)), pt.bg=xc[c(2:8,1)],
       pch=22, title="Number of SNPs", bty="n", ncol=8,
       pt.cex=1, x.intersp = 0.75, cex = 0.8)

pvalueprop<-c(p$pvalue[1:8], NA,
              p$pvalue[9:16], NA, 
              p$pvalue[17:24], NA,  
              p$pvalue[25:32], NA,  
              p$pvalue[33:40], NA,  
              p$pvalue[41:48], NA,  
              p$pvalue[49:56])

points(c(1:8, 10:17, 19:26, 28:35, 37:44, 46:53, 55:62), p$pvalue, col="red", pch=18)
lines(1:62,pvalueprop, col="red", lwd=1)

abline(h=0.1287912, col="black", lty=3)
abline(h=0.05, col="red", lty=2)
abline(h=0.2, col="red", lty=2)

legend(34,0.65, pch=18, col="red", lty=1, lwd=1, legend="Type II error", bty="n", cex=0.8, pt.cex=1)

dev.off()



