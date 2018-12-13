#18 run full mrqap on both males and females 

library(SocGen)
library(rsq)
library(sna)

matricize<-function(x){
  g <- igraph::graph.data.frame(x, directed=FALSE)
  m <- igraph::get.adjacency(g, attr=names(x)[3], sparse=FALSE)
  mat <- m[sort(rownames(m)),sort(colnames(m))]
  diag(mat) <- 1
  return(mat)
}

source("Code/14netlogitM.R")


set.seed(1)

options(stringsAsFactors = FALSE)

setwd("C:/Users/froug/Desktop/Real First Chapter")

# load("C:/ZSL/Coancestry/res_social_females_ordered.RData") #in coancestry folder

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

#add fulldata to res

fulldata<-read.table("C:/ZSL/Coancestry/emp4235/RelatednessEstimates.txt", sep=",", row.names=1)

#need to subset full data for subjects
dolphins<-read.table("focal_dolphins.txt")[[1]]
fulldata<-fulldata[which(fulldata[,1] %in% dolphins & fulldata[,2] %in% dolphins),]

res<-list(fulldata)
names(res)[length(res)] <- "s4235_176_1"

rel_est <- list()

# for (i in 1:length(res)) {

i=1

  df<-res[[i]][,c(1,2,10)]
  names(df)<-c("ID1", "ID2", "MEst")
  
  nsnps<-strsplit(names(res)[i], "_")[[1]][1]
  nsnps<-substr(nsnps, 2, nchar(nsnps))
  nind<-strsplit(names(res)[i], "_")[[1]][2]
  df$nind<-nind
  df$nsnps<-nsnps
  
  df<-merge_pairs(df, smv, "ID1", "ID2", all.x=TRUE, all.y=FALSE)
  
  df<-df[which(df$sexpair=="MALEMALE"),]
  
  rel_est[[i]]<-df
  # }

#make matrices and run netlogitM

#think about nested list structure

# for (i in 1:length(rel_est)){

i=1

  cdata <- rel_est[[i]] 
  
  SRIm <- matricize(cdata[, c("ID1", "ID2", "SRI")])
  relm <- matricize(cdata[, c("ID1", "ID2", "MEst")])
  HROm <- matricize(cdata[, c("ID1", "ID2", "HRO")])
  agem <- matricize(cdata[, c("ID1", "ID2", "agediff")])
  # sexmm <- matricize(cdata[, c("ID1", "ID2", "mm")])
  # sexff <- matricize(cdata[, c("ID1", "ID2", "ff")])
  weightsm <- matricize(cdata[, c("ID1", "ID2", "sightings")])
  weights <- gvectorize(weightsm, mode = "graph", diag = FALSE)
  
  lqap<-netlogitM(y=SRIm,
                  x=list(relm, HROm, agem), #remove sex
                  intercept = TRUE,
                  mode = "graph" ,
                  nullhyp = "qapspp" ,
                  diag = FALSE ,
                  test.statistic = "z-value" ,
                  tol = 1e-07 ,
                  reps = 1000,
                  weights=weights) 
  
  rsqmod<-glm(SRI~MEst+HRO+agediff,
              weights = cdata$sightings,
              data=cdata,
              family = binomial(link="logit"))
  
  rsqs<-rsq.partial(rsqmod)
  
  rsqpartial<-rsqs$partial.rsq[1]
  
  proprel<-sum(cdata$MEst>0.0362)/length(cdata$MEst)
  
  
  #now re-run for males
  
  females_full_model<-list(lqap, rsqmod, rsqs, rsqpartial, proprel)
  
  males_full_model<-list(lqap, rsqmod, rsqs, rsqpartial, proprel)
  
  #have to include rsq in here too