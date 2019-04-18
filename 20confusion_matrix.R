#Plot confusion matrix 

#Based on StackOverflow confusion matrix answer by jbaums  

library(SocGen)

load("ped_est.RData")

#cutoff values
b1<-0.0362
b2<-0.0935
b3<-0.1923
b4<-0.3925

cats<-matrix(c(0, b1, 0, 
               b1, b2, 0.0625,
               b2, b3, 0.125,
               b3, b4, 0.25,
               b4, 1, 0.5), byrow = TRUE, 
             nrow=5)

colnames(cats)<-c("lb", "ub", "cat")

ord<-c(7,1,3,5,8,2,4,6)

# windows()

# pdf(file="Figures/Figure_3_R1.pdf", width=6.65354, height=4.1122)
pdf(file="Figures/Figure_3_R1.pdf", width=5.2, height=3.2138)
par(mfrow=c(2,4), mar=c(0,0,1.5,0.2), oma=c(4, 6, 0, 0))

for (i in ord) {

v5<-ped_est[[i]]

v5<-v5[v5$biparental %in% cats[,3],]

v5$class<-cut(v5$MEst, c(0,b1,b2,b3,b4,1), include.lowest = TRUE)

res<-as.matrix(table(v5$biparental, v5$class))

rels<-c("UR", "HC", "HA", "HS", "PO")
nrels<-c(0, 0.0625, 0.125, 0.25, 0.5)

dimnames(res)<-list(nrels, nrels)

#actual (expected) by predicted (observed)

x <- x.orig <- res #changed
x <- log(x + 0.5) * 2.33
x[x < 0] <- NA
x[x > 10] <- 10
diag(x) <- -diag(x)

image(1:ncol(x), 1:ncol(x),
      -(x[, nrow(x):1]), xlab='', ylab='',
      col=colorRampPalette(c(hsv(h = 0, s = 0.9, v = 0.9, alpha = 1), 
                             hsv(h = 0, s = 0, v = 0.9, alpha = 1), 
                             hsv(h = 0.6, s = 0.9, v = 0.9, alpha = 1)))(41), 
      xaxt='n', yaxt='n', zlim=c(-10, 10))

if(i %in% c(8,2,4,6)){
axis(1, at=1:ncol(x), labels=colnames(x), cex.axis=0.90)
axis(1, at=3, labels=0.125, cex.axis=0.90)}

if(i %in% c(7,8)){
axis(2, at=ncol(x):1, labels=colnames(x), las=1, cex.axis=0.90)
  }


abline(h = 0:ncol(x) + 0.5, col = 'gray')
abline(v = 0:ncol(x) + 0.5, col = 'gray')
text(1:5, rep(5:1, each=5), 
     labels = sub('^0$', '', round(c(x.orig), 0)), cex=1)
box(lwd=1)
text(3, 5.8, v5[1,"nsnps"], xpd=NA, cex=1)

}

mtext("Relationship Category Assignment from Pedigree",1, outer=TRUE, line=3, cex=0.70)
mtext("Relationship Category Assignment from Genetics",2, outer=TRUE, line=4, cex=0.70)

dev.off()