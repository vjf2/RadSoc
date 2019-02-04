#Make pirate plots for subsampled genetic values vs pedigree estimates

#had to fork yarrr to fix horrible par choice

# devtools::install_github("vjf2/yarrr")
library(yarrr)
library(colorspace)

catcol<-rainbow_hcl(12, alpha=0.5)[c(1,3,5,8,10)]

#get ped_est list object from 16pedsubsamples.R

load("ped_est.RData")

pp<-function(df, xaxt=NULL, yaxt=NULL, nsnps=nsnps, ylim=ylim){
  pirateplot(formula = MEst ~ biparental,
             data = df[which(df$biparental %in% c(0, 0.0625, 0.125, 0.25, 0.5)),],
             xlab = NA,
             ylab = NA,
             xaxt = xaxt,
             yaxt = yaxt,
             ylim = ylim,
             xlim = c(0.6, 5.4),
             gl.y = NA,
             bean.b.col=NA,
             bean.f.col=catcol,
             bean.f.o=0.3,
             inf.f.o = 0, # Turn off inf fill
             inf.b.o = 0, # Turn off inf border
             avg.line.lwd = 1,
             cap.beans = TRUE,
             point.cex = 0.75, 
             theme=3)
  
  points(x=rep(1.5,5), y=df$MEst[df$biparental==0.03125], col=adjustcolor("black", alpha.f=0.2), pch=16, cex=0.75)
  
  points(x=rep(2.5,2), y=df$MEst[df$biparental==0.09375], col=adjustcolor("black", alpha.f=0.2), pch=16, cex=0.75)
  
  points(x=rep(3.5,5), y=df$MEst[df$biparental==0.1875], col=adjustcolor("black", alpha.f=0.2), pch=16, cex=0.75)
  
  points(x=rep(4.5,6), y=df$MEst[df$biparental==0.375], col=adjustcolor("black", alpha.f=0.2), pch=16, cex=0.75)
  
  text(1, max(ylim)-.05, nsnps, cex=1)
  
}

#will need to fix aspect ratios, etc. in final pdf plot

lvec<-c(rep(1,57), 
        rep(5, 43),
        rep(2, 57), 
        rep(6, 43), 
        rep(3, 57), 
        rep(7, 43), 
        rep(4, 57), 
        rep(8, 43))

lmat<-matrix(lvec, nrow=100, ncol=4)

windows()
#pdf(file="Figures/Figure_2.pdf", width=6.65354, height=4.132634)
par(mar=c(0,0,0,0), oma=c(4,3.5,0.5,0.1))

layout(lmat)

pp(ped_est[[7]], "n", "n", nsnps=50, ylim=c(-0.02, 0.81))
axis(2, las=1, at=c(0,0.1,0.2,0.3,0.4,0.5, 0.6, 0.7, 0.8), cex.axis=0.8, srt=45, xpd=TRUE)
pp(ped_est[[1]], "n", "n", nsnps=100, ylim=c(-0.02, 0.81))
pp(ped_est[[3]], "n", "n", nsnps=200, ylim=c(-0.02, 0.81))
pp(ped_est[[5]], "n", "n", nsnps=400, ylim=c(-0.02, 0.81))
pp(ped_est[[8]], "n", "n", nsnps=800, ylim=c(-0.02, 0.61))
axis(1, at=c(1,3,5), labels=c(0,0.125,0.5), cex.axis=0.8)
axis(1, at=c(2,4), labels=c(0.0625,0.25), cex.axis=0.8)
axis(2, las=1, at=c(0,0.1,0.2,0.3,0.4,0.5, 0.6), cex.axis=0.8)
pp(ped_est[[2]], "n", yaxt="n", nsnps=1600, ylim=c(-0.02, 0.61))
axis(1, at=c(1,3,5), labels=c(0,0.125,0.5), cex.axis=0.8)
axis(1, at=c(2,4), labels=c(0.0625,0.25), cex.axis=0.8)

pp(ped_est[[4]], "n", yaxt="n", nsnps=3200, ylim=c(-0.02, 0.61))
axis(1, at=c(1,3,5), labels=c(0,0.125,0.5), cex.axis=0.8)
axis(1, at=c(2,4), labels=c(0.0625,0.25), cex.axis=0.8)

pp(ped_est[[6]], "n", yaxt="n", nsnps=4235, ylim=c(-0.02, 0.61))
axis(1, at=c(1,3,5), labels=c(0,0.125,0.5), cex.axis=0.8)
axis(1, at=c(2,4), labels=c(0.0625,0.25), cex.axis=0.8)


mtext("Pedigree Relatedness",1, outer=TRUE, line=2.5, cex=0.75)
mtext("Genetic Relatedness",2, outer=TRUE, line=2.5, cex=0.75)

dev.off()




