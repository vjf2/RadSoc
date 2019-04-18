
library(gridExtra)
library(latticeExtra)
library(reshape2)

ML<-read.csv("Mean_error_MEst.csv")
QG<-read.csv("Mean_error_QGEst.csv")

stats<-aggregate(rmse~nsnps+nind, data=ML, median) #try mean here too

s2<-with(stats, acast(stats, nsnps~nind))

#fill in real data number
s2[nrow(s2), ncol(s2)]<-0

errorcol="blue"

# The plot

MLplot<-cloud(s2, panel.3d.cloud = panel.3dbars, col="gray48",                     
      xbase = 1, ybase = 1, 
      # zlim = c(0, max(s2)+.01),  
      zlim = c(0, 0.25), 
      scales = list(arrows = FALSE, cex=0.4, distance=1.2, tck=1.3), 
      #shade=TRUE,
      lab.cex=0.5,
      cex.axis=0.5,
      xlab = "SNPs", ylab = "Individuals", zlab=list("Median RMSE",rot=90),
      col.facet = level.colors(s2, at = do.breaks(range(s2), 32),
                               col.regions = colorRampPalette(c("white", errorcol))(32),          
                               colors = TRUE),
      colorkey = list(col = colorRampPalette(c("white", errorcol))(32), at = do.breaks(range(s2), 32)),
      screen = list(z = -117, x = -55), 
      zoom=0.9 
      #,main="Method of moments mean error"
) + layer(panel.text(x=0.16, y=-0.19,labels="ML", cex=1))



#write out s2 matrices for all estimators
#need to compare rmse in simulated data with real for recommendations
stats<-aggregate(rmse~nsnps+nind, data=QG, median) #try mean here too

s2<-with(stats, acast(stats, nsnps~nind))

#fill in real data number
s2[nrow(s2), ncol(s2)]<-0

errorcol="darkgreen"

QGplot<-cloud(s2, panel.3d.cloud = panel.3dbars, col="gray48",                     
          xbase = 1, ybase = 1, 
          # zlim = c(0, max(s2)+.01),  
          zlim = c(0, 0.25), 
          scales = list(arrows = FALSE, cex=0.4, distance=1.2, tck=1.3), 
          #shade=TRUE,
          cex=0.5,
          xlab = "SNPs", ylab = "Individuals", zlab=list("Median RMSE",rot=90),
          col.facet = level.colors(s2, at = do.breaks(range(s2), 32),
                                   col.regions = colorRampPalette(c("white", errorcol))(32),          
                                   colors = TRUE),
          colorkey = list(col = colorRampPalette(c("white", errorcol))(32), at = do.breaks(range(s2), 32)),
          screen = list(z = -117, x = -55), 
          #perspective= 
          zoom=0.9 
) + layer(panel.text(x=0.16, y=-0.19,labels="QG", cex=1))


 CEX <- 0.7

 pset <- list(
             # par.sub.text     = list(cex = CEX),
             # par.main.text    = list(cex = CEX),
             par.zlab.text    = list(cex = CEX),
             par.ylab.text    = list(cex = CEX),
             par.xlab.text    = list(cex = CEX),
             # add.text         = list(cex = CEX),
             axis.text        = list(cex = CEX))

MLplot <- update(MLplot, par.settings = pset)
QGplot <- update(QGplot, par.settings = pset)

# MLplot<-update(MLplot, par.settings = list(fontsize = list(text = 8)))
# QGplot<-update(QGplot, par.settings = list(fontsize = list(text = 8)))

pdf(file="Figures/Figure_1_R1.pdf",width=6.65354, height=4.132634)
grid.arrange(MLplot, QGplot, ncol=2)
dev.off()
