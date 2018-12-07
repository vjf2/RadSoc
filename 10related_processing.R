options(stringsAsFactors = FALSE)

library(SocGen)

#load R object 

load(file="res.RData")

nb<-length(res)

fulldata<-read.table("C:/ZSL/Coancestry/emp4235/RelatednessEstimates.txt", sep=",", row.names=1)

names(fulldata)<-c("Ind1","Ind2","Pair","TrioEst","WEst","LLEst","LREst","REst","QGEst","MEst")

fulldata<-fulldata[,c(1:2,10)]
names(fulldata)[3]<-"full"

master_IBD<-list() #takes about 2 min

start<-Sys.time()

for (j in 1:nb) {
  
  dis<-res[[j]]
  names(dis)<-c("Ind1","Ind2","Pair","TrioEst","WEst","LLEst","LREst","REst","QGEst","MEst")
  dis<-dis[,-ncol(dis)]
  
  # names(res)[j]
  
  nsnps<-strsplit(names(res)[j], "_")[[1]][1]
  nsnps<-substr(nsnps, 2, nchar(nsnps))
  nind<-strsplit(names(res)[j], "_")[[1]][2]
  dis$nind<-nind
  dis$nsnps<-nsnps
  dis$iter<-j

  dis<-merge_pairs(dis, fulldata, "Ind1", "Ind2", all.x=TRUE, all.y=FALSE)
  
  dis$err<-abs(dis$MEst-dis$full)
  
  rmse<-sqrt((1/nrow(dis))*sum(dis$err^2))
  
  final<-c(dis$nsnps[1], dis$nind[1], rmse)
  
  master_IBD[[j]]<-final
}

end<-Sys.time()

fmaster_IBD<-do.call("rbind", master_IBD)
colnames(fmaster_IBD)<-c("nsnps", "nind", "rmse")
fmaster_IBD<-apply(fmaster_IBD, 2, as.numeric)

write.csv(fmaster_IBD, "Mean_error_MEst.csv", row.names = FALSE)

#try to make plot

library(latticeExtra)
library(reshape2)

stats<-aggregate(rmse~nsnps+nind, data=fmaster_IBD, median) #try mean here too

s2<-with(stats, acast(stats, nsnps~nind))

#fill in real data number
s2[nrow(s2), ncol(s2)]<-0

errorcol="darkgreen"

# The plot
windows()
#pdf(file="median_rmse_lynchli.pdf")
cloud(s2, panel.3d.cloud = panel.3dbars, col="gray48",                     
      xbase = 1, ybase = 1, 
      # zlim = c(0, max(s2)+.01),  
      zlim = c(0, 0.25), 
      scales = list(arrows = FALSE, just = "right"), 
      #shade=TRUE,
      #  lab.cex=0.5,
      xlab = "SNPs", ylab = "Individuals", zlab=list("Median Error",rot=90),
      col.facet = level.colors(s2, at = do.breaks(range(s2), 32),
                               col.regions = colorRampPalette(c("white", errorcol))(32),          
                               colors = TRUE),
      colorkey = list(col = colorRampPalette(c("white", errorcol))(32), at = do.breaks(range(s2), 32)),
      screen = list(z = -117, x = -55), 
      #perspective= 
      zoom=0.9 
      #,main="Method of moments mean error"
) + layer(panel.text(x=-0.17, y=0.16,labels="QG", cex=2))

dev.off()

#write out s2 matrices for all estimators
#need to compare rmse in simulated data with real for recommendations



