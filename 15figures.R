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

windows()

# pdf(file="Figures/model_res100.pdf", width=14, height=6, useDingbats=FALSE)
par(mar=c(5, 4.1, 1, 1))
boxplot(rsqrel~nsnps+nind,  
        col=adjustcolor(xc, alpha.f=0.5)[c(2:8,1)], 
        at=c(1:8, 10:17, 19:26, 28:35, 37:44, 46:53, 55:62),
        data=qdf,
        xaxt="n",
        yaxt="n",
        ylab="Coefficient of Partial Determination - Relatedness",
        xlab="Number of Individuals",
        ylim=c(min(qdf$rsqrel),0.72), 
        type="n",
        outcex=0.75, 
        outcol="gray30") 

axis(2, las=1)
axis(1, at=seq(4,62,9), labels=unique(qdf$nind))

legend(55,0.70, legend=sort(unique(qdf$nsnps)), pt.bg=xc[c(2:8,1)], pch=22,
       title="Number of SNPs", y.intersp = 0.9, bty="n")

pvalueprop<-c(p$pvalue[1:8], NA,
              p$pvalue[9:16], NA, 
              p$pvalue[17:24], NA,  
              p$pvalue[25:32], NA,  
              p$pvalue[33:40], NA,  
              p$pvalue[41:48], NA,  
              p$pvalue[49:56])

points(c(1:8, 10:17, 19:26, 28:35, 37:44, 46:53, 55:62), p$pvalue, col="red", pch=18)
lines(1:62,pvalueprop, col="red", lwd=1)

abline(h=0.1287912, col="darkgrey", lty=2)
abline(h=0.05, col="black", lty=2)
abline(h=0.2, col="black", lty=3)

# boxplot(rsqrel~nsnps+nind,  
#         col=adjustcolor(xc, alpha.f=0.7)[c(2:8,1)], 
#         at=c(1:8, 10:17, 19:26, 28:35, 37:44, 46:53, 55:62),
#         data=qdf,
#         xaxt="n",
#         yaxt="n",
#         add=TRUE,
#         outcex=0.75) 

legend(54,0.38, pch=18, col="red", lty=1, lwd=1, legend="Type II error", bty="n")

dev.off()


#False Discovery Rate?
windows()

par(mar=c(5, 4.1, 1, 1))
boxplot(proprel~nsnps+nind,  
        col=adjustcolor(xc, alpha.f=0.5)[c(2:8,1)], 
        at=c(1:8, 10:17, 19:26, 28:35, 37:44, 46:53, 55:62),
        data=qdf,
        xaxt="n",
        yaxt="n",
        ylab="Coefficient of Partial Determination - Relatedness",
        xlab="Number of Individuals",
        ylim=c(0,0.7), 
        type="n") 

axis(2, las=1)
axis(1, at=seq(4,62,9), labels=unique(qdf$nind))

legend(55,0.70, legend=sort(unique(qdf$nsnps)), col=xc[c(2:8,1)], pch=15,
       title="Number of SNPs", y.intersp = 0.9, bty="n")





