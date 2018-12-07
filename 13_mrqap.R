#13 use res_social.Data in MRQAP models

library(SocGen)

options(stringsAsFactors = FALSE)

setwd("C:/Users/froug/Desktop/Real First Chapter")

load("C:/ZSL/Coancestry/res_social.RData") #in coancestry folder

#read in social data

smv<-read.csv("social_model_variabes.csv", row.names=1)

smv$together<-smv$sightings*smv$SRI

#make 0 and weight of 1 for those without overlap
smv$together[is.na(smv$together)]<-0
smv$SRI[is.na(smv$SRI)]<-0
smv$sightings[which(smv$sightings==0)]<-1
smv$mm<-ifelse(smv$sexpair=="MALEMALE", 1, 0)
smv$ff<-ifelse(smv$sexpair=="FEMALEFEMALE", 1, 0)

#add names to res

names(res) <- gsub(pattern = "/RelatednessEstimates.Txt", "", names(res))

#remove nind = 10

res<-res[grep("_10_1", x = names(res), invert = TRUE)]

#add fulldata to res

fulldata<-read.table("C:/ZSL/Coancestry/emp4235/RelatednessEstimates.txt", sep=",", row.names=1)

#need to subset full data for subjects
dolphins<-read.table("focal_dolphins.txt")[[1]]
fulldata<-fulldata[which(fulldata[,1] %in% dolphins & fulldata[,2] %in% dolphins),]

res<-c(res, list(fulldata))
names(res)[length(res)] <- "s4235_176_1"

rel_est <- list()

for (i in 1:length(res)) {
  df<-res[[i]][,c(1,2,10)]
  names(df)<-c("ID1", "ID2", "MEst")
  
  nsnps<-strsplit(names(res)[i], "_")[[1]][1]
  nsnps<-substr(nsnps, 2, nchar(nsnps))
  nind<-strsplit(names(res)[i], "_")[[1]][2]
  df$nind<-nind
  df$nsnps<-nsnps
  
  df<-merge_pairs(df, smv, "ID1", "ID2", all.x=TRUE, all.y=FALSE)
  
  rel_est[[i]]<-df}

#make matrices and run netlogitM

library(sna)

matricize<-function(x){
  g <- igraph::graph.data.frame(x, directed=FALSE)
  m <- igraph::get.adjacency(g, attr=names(x)[3], sparse=FALSE)
  mat <- m[sort(rownames(m)),sort(colnames(m))]
  diag(mat) <- 1
  return(mat)
}

source("Code/14netlogitM.R")

lqap <- list()
lqapH <- list()
llogmod <- list()
llogmodH <- list()
llogmodR <- list()
llogmodRH <- list()
llogmodtar <- list()
llogmodtarH <- list()

start<-Sys.time()

for (i in 1:length(rel_est)){
  
  cdata <- rel_est[[i]]
  
  SRIm <- matricize(cdata[, c("ID1", "ID2", "SRI")])
  relm <- matricize(cdata[, c("ID1", "ID2", "MEst")])
  HROm <- matricize(cdata[, c("ID1", "ID2", "HRO")])
  agem <- matricize(cdata[, c("ID1", "ID2", "agediff")])
  sexmm <- matricize(cdata[, c("ID1", "ID2", "mm")])
  sexff <- matricize(cdata[, c("ID1", "ID2", "ff")])
  weightsm <- matricize(cdata[, c("ID1", "ID2", "sightings")])
  weights <- gvectorize(weightsm, mode = "graph", diag = FALSE)
 
  lqap[[i]]<-netlogitM(y=SRIm,
                  x=list(relm, HROm, agem, sexmm, sexff),
                  intercept = TRUE,
                  mode = "graph" ,
                  nullhyp = "qapspp" ,
                  diag = FALSE ,
                  test.statistic = "z-value" ,
                  tol = 1e-07 ,
                  reps = 1000,
                  weights=weights) 
  
  llogmod[[i]]<-glm(cbind(together, (sightings-together))~MEst+HRO+agediff+mm+ff, data=cdata, family = binomial(link="logit"))
  
  llogmodR[[i]]<-glm(cbind(together, (sightings-together))~HRO+agediff+mm+ff, data=cdata, family = binomial(link="logit"))
  
  llogmodtar[[i]]<-glm(MEst~HRO+agediff+mm+ff, data=cdata, family = binomial(link="logit"))
  
  #without HRO
  
  lqapH[[i]]<-netlogitM(y=SRIm,
                 x=list(relm, agem, sexmm, sexff),
                 intercept = TRUE,
                 mode = "graph" ,
                 nullhyp = "qapspp" ,
                 diag = FALSE ,
                 test.statistic = "z-value" ,
                 tol = 1e-07 ,
                 reps = 1000,
                 weights=weights) 
  
  llogmodH[[i]]<-glm(cbind(together, (sightings-together))~MEst+agediff+mm+ff, data=cdata, family = binomial(link="logit"))
  
  llogmodRH[[i]]<-glm(cbind(together, (sightings-together))~agediff+mm+ff, data=cdata, family = binomial(link="logit"))
  
  llogmodtarH[[i]]<-glm(MEst~agediff+mm+ff, data=cdata, family = binomial(link="logit"))
  
  save.image(file=paste0("all_model_results",i,".RData"))
  print(i)
}

end<-Sys.time()





