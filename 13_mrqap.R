#13 use res_social.Data in MRQAP models

arg1<-commandArgs(trailingOnly = TRUE)
arg1<-as.numeric(arg1)

.libPaths("C:/Users/froug/OneDrive/Documents/R/win-library/3.5")

library(SocGen)
library(rsq)

options(stringsAsFactors = FALSE)

xfiles<-list.files(path="C:/ZSL/Coancestry", pattern="RelatednessEstimates", recursive = TRUE, full.names = TRUE)
xfiles<-xfiles[grep("social", xfiles)]
xfiles<-xfiles[grep(paste0("_", arg1, "/"), xfiles, fixed=TRUE)]

# #check and remove any other filenames

res<-lapply(xfiles, function(x) read.table(x, sep=",", row.names = 1))
names(res)<-xfiles

# save(res, file="res_social_females_ordered.RData")
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

names(res) <- gsub(pattern = "/RelatednessEstimates.Txt", "", names(res))

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
  
  df<-df[which(df$sexpair=="FEMALEFEMALE"),]
  
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
rsqpartial<-list()
proprel<-list()

set.seed(1)


library(foreach)
library(doParallel)

start<-Sys.time()

cl<-makeCluster(detectCores()-1, outfile="keep_track.txt")
clusterEvalQ(cl, library(sna))
clusterEvalQ(cl, library(rsq))
clusterExport(cl, c("rel_est", "matricize", "netlogitM", "lqap", "rsqpartial", "proprel"))
registerDoParallel(cl)

#think about nested list structure

# for (i in 1:length(rel_est)){
  
models<-foreach (i=1:length(rel_est), .errorhandling='pass') %dopar% {
  
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

  
  cat(paste0("network regression complete for number ", i, "\n"))
  
  return(list(lqap, rsqpartial, proprel))
  
} 

stopCluster(cl)

end<-Sys.time()

end-start

# save.image(file=paste0("all_female_model_results_seed",3,".RData"))

#make results table 

restable<-lapply(names(res), function(x) {
            nsnps<-strsplit(x, "_")[[1]][1]
            nsnps<-as.numeric(gsub("\\D", "", nsnps))
            nind<-strsplit(x, "_")[[1]][2]
            return(c(nind, nsnps))
          })
          
restable<-do.call("rbind", restable)
restable<-as.data.frame(apply(restable,2,as.numeric))

names(restable)<-c("nind", "nsnps")
restable[nrow(restable),1]<-92

lqap<-lapply(models, "[[", 1)

restable$pvalue<-sapply(lqap, function(x) x$pgreqabs[2]) 
restable$coeff<-sapply(lqap, function(x) coef(x)[2])

rsqpartial<-lapply(models, "[[", 2)
proprel<-lapply(models, "[[",3)

restable$rsqrel<-unlist(rsqpartial)
restable$proprel<-unlist(proprel)


write.csv(restable, file=paste0("restable", arg1, ".csv"))


