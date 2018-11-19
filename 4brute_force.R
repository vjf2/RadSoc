estimator<-"MEst"

v5<-rels
v5$IBD<-v5[,estimator]

ur<-v5$IBD[v5$expected==0]
hc<-v5$IBD[v5$expected==0.0625]
fc<-v5$IBD[v5$expected==0.125]

x1<-max(ur[which(ur<min(hc))])
x2<-min(hc[which(hc>max(ur))])

x3<-max(hc[which(hc<min(fc))])
x4<-min(fc[which(fc>max(hc))])

lt<-seq(x1, x2, 0.0005)
ut<-seq(x3, x4, 0.0005)

thresh_options<-expand.grid(lt, ut)

correct<-list()
incorrect<-list()


for (i in 1:nrow(thresh_options)) {
  
  v5$classification<-cut(v5$IBD, c(thresh_options[i,1], thresh_options[i,2]), include.lowest = TRUE)
  
  tabl<-as.data.frame(table(v5$classification, v5$RelCat))
  
  correct[[i]]<-tabl[tabl$Var2=="HC", "Freq"]
  incorrect[[i]]<-sum(tabl[!tabl$Var2=="HC", "Freq"])
  
}

res<-data.frame(thresh_options, tp=unlist(correct), fp=unlist(incorrect))

res$fpp<-res[,4]/(res[,3]+res[,4])

co1<-res[which(res$fpp<=0.05),]

co1<-co1[co1$tp==max(co1$tp),]

co1$width<-co1$Var2-co1$Var1

#add width column

ltfc<-seq(min(co1$Var2), x4, 0.0005)
utfc<-max(fc)


thresh_options<-expand.grid(ltfc, utfc)

correct<-list()
incorrect<-list()


for (i in 1:nrow(thresh_options)) {
  
  v5$classification<-cut(v5$IBD, c(thresh_options[i,1], thresh_options[i,2]), include.lowest = TRUE)
  
  tabl<-as.data.frame(table(v5$classification, v5$RelCat))
  
  correct[[i]]<-tabl[tabl$Var2=="FC", "Freq"]
  incorrect[[i]]<-sum(tabl[!tabl$Var2=="FC", "Freq"])
  
}

res<-data.frame(thresh_options, tp=unlist(correct), fp=unlist(incorrect))

res$fpp<-res[,4]/(res[,3]+res[,4])

co2<-res[which(res$fpp<=0.05),]

co2<-co2[co2$tp==max(co2$tp),]

co2$width<-co2$Var2-co2$Var1


ltfc<-min(ur)
utfc<-seq(x1, max(co1$Var1), 0.0005)


thresh_options<-expand.grid(ltfc, utfc)

correct<-list()
incorrect<-list()


for (i in 1:nrow(thresh_options)) {
  
  v5$classification<-cut(v5$IBD, c(thresh_options[i,1], thresh_options[i,2]), include.lowest = TRUE)
  
  tabl<-as.data.frame(table(v5$classification, v5$RelCat))
  
  correct[[i]]<-tabl[tabl$Var2=="UR", "Freq"]
  incorrect[[i]]<-sum(tabl[!tabl$Var2=="UR", "Freq"])
  
}

res<-data.frame(thresh_options, tp=unlist(correct), fp=unlist(incorrect))

res$fpp<-res[,4]/(res[,3]+res[,4])

co<-res[which(res$fpp<=0.05),]

co<-co[co$tp==max(co$tp),]

co$width<-co$Var2-co$Var1

#combine co, co1, and co2

# rect(xleft = 0.0357, ybottom = 0 , xright = 0.0372, ytop = 60 , col = adjustcolor("darkgrey", alpha=0.9), border=NA)
# 
# rect(xleft = 0.0901, ybottom = 0 , xright = 0.0976, ytop = 60 , col = adjustcolor("darkgrey", alpha=0.9), border=NA)

text(x = 0.021, y=40, "97%", cex=0.8)
text(x = 0.065, y=40, "85%", cex=0.8)
text(x = 0.13, y=40, "96%", cex=0.8)
text(x = 0.25, y=40, "100%", cex=0.8)
text(x = 0.45, y=40, "100%", cex=0.8)

text(x = 0.021, y=37, "5%", col="red", cex=0.8) #false positive
text(x = 0.065, y=37, "5%", col="red", cex=0.8) #false positive
text(x = 0.13, y=37, "5%", col="red", cex=0.8) #false positive
text(x = 0.25, y=37, "0%", col="red", cex=0.8) #false positive
text(x = 0.45, y=37, "0%", col="red", cex=0.8) #false positive

title(main="Simulated Data")

b1<-0.0362
b2<-0.0935
b3<-(max(v5$IBD[which(v5$expected==0.125)])+min(v5$IBD[which(v5$expected==0.25)]))/2
b4<-(max(v5$IBD[which(v5$expected==0.25)])+min(v5$IBD[which(v5$expected==0.5)]))/2

segments(x0=b1, y0=0, x1=b1, y1=60, lty=3)
segments(x0=b2, y0=0, x1=b2, y1=60, lty=3)
segments(x0=b3, y0=0, x1=b3, y1=60, lty=3)
segments(x0=b4, y0=0, x1=b4, y1=60, lty=3)

