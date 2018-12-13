v5<-ref
v5$IBD<-v5$relatedness.x
v5$expected<-v5$biparental

v5<-v5[which(v5$expected %in% c(0, 0.125, 0.25, 0.5)),]


ur<-v5$IBD[v5$expected==0]
fc<-v5$IBD[v5$expected==0.125]
hs<-v5$IBD[v5$expected==0.25]

x1<-max(ur[which(ur<min(fc))])
x2<-min(fc[which(fc>max(ur))])

x3<-max(fc[which(fc<min(hs))])
x4<-min(hs[which(hs>max(fc))])

lt<-seq(x1, x2, 0.0005)
ut<-seq(x3, x4, 0.0005)

thresh_options<-expand.grid(lt, ut)

correct<-list()
incorrect<-list()


for (i in 1:nrow(thresh_options)) {
  
  v5$classification<-cut(v5$IBD, c(thresh_options[i,1], thresh_options[i,2]), include.lowest = TRUE)
  
  tabl<-as.data.frame(table(v5$classification, v5$expected))
  
  correct[[i]]<-tabl[tabl$Var2==0.125, "Freq"]
  incorrect[[i]]<-sum(tabl[!tabl$Var2==0.125, "Freq"])
  
}

res<-data.frame(thresh_options, tp=unlist(correct), fp=unlist(incorrect))

res$fpp<-res[,4]/(res[,3]+res[,4])

co1<-res[which(res$fpp<=0.05),]

co1<-co1[co1$tp==max(co1$tp),]

co1$width<-co1$Var2-co1$Var1

#add width column

lths<-seq(min(co1$Var2), x4, 0.0005)
uths<-max(hs)


thresh_options<-expand.grid(lths, uths)

correct<-list()
incorrect<-list()


for (i in 1:nrow(thresh_options)) {
  
  v5$classification<-cut(v5$IBD, c(thresh_options[i,1], thresh_options[i,2]), include.lowest = TRUE)
  
  tabl<-as.data.frame(table(v5$classification, v5$expected))
  
  correct[[i]]<-tabl[tabl$Var2==0.25, "Freq"]
  incorrect[[i]]<-sum(tabl[!tabl$Var2==0.25, "Freq"])
  
}

res<-data.frame(thresh_options, tp=unlist(correct), fp=unlist(incorrect))

res$fpp<-res[,4]/(res[,3]+res[,4])

co2<-res[which(res$fpp<=0.05),]

co2<-co2[co2$tp==max(co2$tp),]

co2$width<-co2$Var2-co2$Var1


lths<-min(ur)
uths<-seq(x1, max(co1$Var1), 0.0005)


thresh_options<-expand.grid(lths, uths)

correct<-list()
incorrect<-list()


for (i in 1:nrow(thresh_options)) {
  
  v5$classification<-cut(v5$IBD, c(thresh_options[i,1], thresh_options[i,2]), include.lowest = TRUE)
  
  tabl<-as.data.frame(table(v5$classification, v5$expected))
  
  correct[[i]]<-tabl[tabl$Var2==0, "Freq"]
  incorrect[[i]]<-sum(tabl[!tabl$Var2==0, "Freq"])
  
}

res<-data.frame(thresh_options, tp=unlist(correct), fp=unlist(incorrect))

res$fpp<-res[,4]/(res[,3]+res[,4])

co<-res[which(res$fpp<=0.05),]

co<-co[co$tp==max(co$tp),]

co$width<-co$Var2-co$Var1

#combine co, co1, and co2

# rect(xleft = 0.0362, ybottom = 0 , xright = 0.0935, ytop = 60 , col = adjustcolor("darkgrey", alpha=0.9), border=NA)
# 
# rect(xleft = 0.0901, ybottom = 0 , xright = 0.0976, ytop = 80 , col = adjustcolor("darkgrey", alpha=0.9), border=NA)

tcex=1

text(x = 0.021, y=40, "100%", cex=1)
text(x = 0.13, y=40, "80%", cex=1)
text(x = 0.25, y=40, "98%", cex=1)
text(x = 0.45, y=40, "100%", cex=1)

text(x = 0.021, y=37, "0%", col="red", cex=1) #false positive
text(x = 0.13, y=37, "5%", col="red", cex=1) #false positive
text(x = 0.25, y=37, "2%", col="red", cex=1) #false positive
text(x = 0.45, y=37, "0%", col="red", cex=1) #false positive

# title(main="Real Data")

b1<-0.0362
b2<-0.0935
b3<-(max(v5$IBD[which(v5$expected==0.125)])+min(v5$IBD[which(v5$expected==0.25)]))/2
b4<-(max(v5$IBD[which(v5$expected==0.25)])+min(v5$IBD[which(v5$expected==0.5)]))/2

segments(x0=b1, y0=0, x1=b1, y1=80, lty=3)
segments(x0=b2, y0=0, x1=b2, y1=80, lty=3)
segments(x0=b3, y0=0, x1=b3, y1=80, lty=3)
segments(x0=b4, y0=0, x1=b4, y1=80, lty=3)

rect(xleft = b1, ybottom = 0 , xright = b2, ytop = 80 , col = adjustcolor("darkgrey", alpha=0.9), border=NA)


#what if we used categories from simulated data? 

# 0-0.362
#0.935-0.18
#0.18-0.4


correct<-list()
incorrect<-list()

i=1
thresh_options<-data.frame(0.0935, 0.175)

v5$classification<-cut(v5$IBD, c(thresh_options[i,1], thresh_options[i,2]), include.lowest = TRUE)

tabl<-as.data.frame(table(v5$classification, v5$expected))

correct[[i]]<-tabl[tabl$Var2==0.125, "Freq"]
incorrect[[i]]<-sum(tabl[!tabl$Var2==0.125, "Freq"])
res<-data.frame(thresh_options, tp=unlist(correct), fp=unlist(incorrect))



correct<-list()
incorrect<-list()

i=1
thresh_options<-data.frame(0.175, 0.4)

v5$classification<-cut(v5$IBD, c(thresh_options[i,1], thresh_options[i,2]), include.lowest = TRUE)

tabl<-as.data.frame(table(v5$classification, v5$expected))

correct[[i]]<-tabl[tabl$Var2==0.25, "Freq"]
incorrect[[i]]<-sum(tabl[!tabl$Var2==0.25, "Freq"])
res<-data.frame(thresh_options, tp=unlist(correct), fp=unlist(incorrect))


