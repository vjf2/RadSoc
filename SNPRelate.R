#Calculate MLE for subsamples 

library(gdsfmt)
library(SNPRelate)
library(parallel)

#use --make-bed command in plink

setwd("C:/Users/froug/Desktop/Real First Chapter")

options(stringsAsFactors = FALSE)

bed.fn <- "original2018p.bed"
fam.fn <- "original2018p.fam"
bim.fn <- "original2018p.bim"

#snpgdsCreateGeno

snpgdsBED2GDS(bed.fn, fam.fn, bim.fn, "original2018p.gds")

snpgdsSummary("original2018p.gds")

genofile <- snpgdsOpen("original2018p.gds")

sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
snp.id <- read.gdsn(index.gdsn(genofile, "snp.id"))

#sample different numbers of individuals and snps

nmax<-length(sample.id)
nsnps<-length(snp.id)

snpset<-c(50,100,200,400,800,1600,3200,nsnps)
nind<-c(10,20,30,40,50,100,180)

#set up different loops for nind=200 and nind=nmax to force include behavioral IDs

#da megaloop

combns<-expand.grid(snpset, nind)
names(combns)<-c("nsnps", "nind")

combns<-combns[-nrow(combns),]

#get subjects from Social_model.R

ID_key<-read.csv("RawData/dolphin_sample_key.csv", colClasses = "character")

ID_key$Npoints<-as.numeric(ID_key$Npoints)

subjects<-ID_key$Dolphin_ID[which(ID_key$me_paper=="yes" & ID_key$Npoints>=35)]
subjects<-unique(subjects)
nsub<-length(subjects)

nt<-detectCores()-1

snp_order<-paste0("G", 1:4124)

subCalc<-function(num_snps, num_samp, genofile, nt){
  snp.ids<-sample(snp.id, size=num_snps)
  snp.ids<-snp.ids[order(match(snp.ids,snp_order))]
  sample.ids<-sample(subjects, size=num_samp)
  
  afreq <- snpgdsSNPRateFreq(genofile, snp.id=snp.ids, sample.id=sample.ids)$AlleleFreq
  
  ibd_MoM<-snpgdsIBDMoM(genofile, sample.id=sample.ids, snp.id=snp.ids, allele.freq=afreq, autosome.only=FALSE, kinship=TRUE, num.thread=nt)
  
  ibd_MLE<-snpgdsIBDMLE(genofile, sample.id=sample.ids, snp.id=snp.ids, allele.freq=afreq, method="EM", autosome.only=FALSE, kinship=TRUE, num.thread=nt)
  
  #get relatedness likelihood matrices
  
  fs_Lik<-snpgdsIBDMLELogLik(genofile, ibd_MLE, relatedness="fullsib")
  po_Lik<-snpgdsIBDMLELogLik(genofile, ibd_MLE, relatedness="offspring")
  hs_Lik<-snpgdsIBDMLELogLik(genofile, ibd_MLE, relatedness="halfsib")
  co_Lik<-snpgdsIBDMLELogLik(genofile, ibd_MLE, relatedness="cousin")
  un_Lik<-snpgdsIBDMLELogLik(genofile, ibd_MLE, relatedness="unrelated")
  
  dimnames(fs_Lik)<-list(sample.ids, sample.ids) 
  dimnames(po_Lik)<-list(sample.ids, sample.ids) 
  dimnames(hs_Lik)<-list(sample.ids, sample.ids) 
  dimnames(co_Lik)<-list(sample.ids, sample.ids) 
  dimnames(un_Lik)<-list(sample.ids, sample.ids) 
  
  #get coeff dataframes
  
  ibd_MoM.coeff <- snpgdsIBDSelection(ibd_MoM)
  ibd_MoM.coeff$pi_hat<-ibd_MoM.coeff$kinship*2
  
  ibd_MLE.coeff <- snpgdsIBDSelection(ibd_MLE)
  ibd_MLE.coeff$pi_hat<-ibd_MLE.coeff$kinship*2
  
  #bundle results

  results<-list(ibd_MoM.coeff, 
                ibd_MLE.coeff, 
                fs_Lik, 
                po_Lik,
                hs_Lik, 
                co_Lik,
                un_Lik)
  
  names(results)<-c(paste0(name_prefix,"ibd_MoM.coeff"),
                  paste0(name_prefix,"ibd_MLE.coeff"),
                  paste0(name_prefix,"fs_Lik"),
                  paste0(name_prefix,"po_Lik"),
                  paste0(name_prefix,"hs_Lik"),
                  paste0(name_prefix,"co_Lik"),
                  paste0(name_prefix,"un_Lik")) 

  return(results)
  
  }

 replicate_set<-list()

for (i in 1:nrow(combns)) {
  num_snps<-combns$nsnps[i]
  num_samp<-combns$nind[i]
  name_prefix<-paste0("S",num_snps, "N", num_samp, "_")
  
  finals<-replicate(10, subCalc(num_snps=num_snps, num_samp=num_samp, genofile=genofile, nt=nt), simplify = FALSE)
  
  names(finals)<-paste0(name_prefix, i)
  
  replicate_set[[i]]<-finals
  
  #save output
  
  save(finals, file=paste0(name_prefix, i, ".RData"))
  
  }

#write one for 200

diff.amount=200-nsub 
diff.poss=setdiff(sample.id, subjects)

comb200<-expand.grid(snpset,200)
names(comb200)<-c("nsnps", "nind")

replicate_set200<-list()

for (i in 1:nrow(comb200)) {
  num_snps<-combns$nsnps[i]
  num_samp<-combns$nind[i]
  name_prefix<-paste0("S",num_snps, "N", num_samp, "_")
  
  finals<-replicate(100, subCalc200(num_snps=num_snps, diff.amount=diff.amount, genofile=genofile, nt=nt), simplify = FALSE)
  
  names(finals)<-paste0(name_prefix, i)
  
  replicate_set200[[i]]<-finals
  
  #save output
  save(finals, file=paste0(name_prefix, i, ".RData"))
}

subCalc200<-function(num_snps, diff.amount, genofile, nt){
  snp.ids<-sample(snp.id, size=num_snps)
  sample.ids<-sample(diff.poss, size=diff.amount)
  
  sample.ids<-c(sample.ids,subjects)
  
  ibd_MoM<-snpgdsIBDMoM(genofile, sample.id=sample.ids, snp.id=snp.ids, autosome.only=FALSE, kinship=TRUE, num.thread=nt)
  
  ibd_MLE<-snpgdsIBDMLE(genofile, sample.id=sample.ids, snp.id=snp.ids, method="EM", autosome.only=FALSE, kinship=TRUE, num.thread=nt)
  
  #get relatedness likelihood matrices
  
  fs_Lik<-snpgdsIBDMLELogLik(genofile, ibd_MLE, relatedness="fullsib")
  po_Lik<-snpgdsIBDMLELogLik(genofile, ibd_MLE, relatedness="offspring")
  hs_Lik<-snpgdsIBDMLELogLik(genofile, ibd_MLE, relatedness="halfsib")
  co_Lik<-snpgdsIBDMLELogLik(genofile, ibd_MLE, relatedness="cousin")
  un_Lik<-snpgdsIBDMLELogLik(genofile, ibd_MLE, relatedness="unrelated")
  
  dimnames(fs_Lik)<-list(sample.ids, sample.ids) 
  dimnames(po_Lik)<-list(sample.ids, sample.ids) 
  dimnames(hs_Lik)<-list(sample.ids, sample.ids)  
  dimnames(co_Lik)<-list(sample.ids, sample.ids) 
  dimnames(un_Lik)<-list(sample.ids, sample.ids) 
  
  #get coeff dataframes
  
  ibd_MoM.coeff <- snpgdsIBDSelection(ibd_MoM)
  ibd_MoM.coeff$pi_hat<-ibd_MoM.coeff$kinship*2
  
  ibd_MLE.coeff <- snpgdsIBDSelection(ibd_MLE)
  ibd_MLE.coeff$pi_hat<-ibd_MLE.coeff$kinship*2
  
  #bundle results
  
  results<-list(ibd_MoM.coeff, 
                ibd_MLE.coeff, 
                fs_Lik, 
                po_Lik,
                hs_Lik, 
                co_Lik,
                un_Lik)
  
  names(results)<-c(paste0(name_prefix,"ibd_MoM.coeff"),
                    paste0(name_prefix,"ibd_MLE.coeff"),
                    paste0(name_prefix,"fs_Lik"),
                    paste0(name_prefix,"po_Lik"),
                    paste0(name_prefix,"hs_Lik"),
                    paste0(name_prefix,"co_Lik"),
                    paste0(name_prefix,"un_Lik")) 
  
  return(results)
  
}

# write one for 275

comb275<-expand.grid(snpset,275)
names(comb275)<-c("nsnps", "nind")
comb275<-comb275[-nrow(comb275),] #only need with both maxs 

replicate_set275<-list()

for (i in 1:nrow(comb200)) {
  num_snps<-combns$nsnps[i]
  num_samp<-combns$nind[i]
  name_prefix<-paste0("S",num_snps, "N", num_samp, "_")
  
  finals<-replicate(100, subCalc275(num_snps=num_snps, genofile=genofile, nt=nt), simplify = FALSE)
  
  names(finals)<-paste0(name_prefix, i)
  
  replicate_set200[[i]]<-finals
  
  #save output
  save(finals, file=paste0(name_prefix, i, ".RData"))
}

subCalc275<-function(num_snps, genofile, nt){
  snp.ids<-sample(snp.id, size=num_snps)
  
  ibd_MoM<-snpgdsIBDMoM(genofile, snp.id=snp.ids, autosome.only=FALSE, kinship=TRUE, num.thread=nt)
  
  ibd_MLE<-snpgdsIBDMLE(genofile, snp.id=snp.ids, method="EM", autosome.only=FALSE, kinship=TRUE, num.thread=nt)
  
  #get relatedness likelihood matrices
  
  fs_Lik<-snpgdsIBDMLELogLik(genofile, ibd_MLE, relatedness="fullsib")
  po_Lik<-snpgdsIBDMLELogLik(genofile, ibd_MLE, relatedness="offspring")
  hs_Lik<-snpgdsIBDMLELogLik(genofile, ibd_MLE, relatedness="halfsib")
  co_Lik<-snpgdsIBDMLELogLik(genofile, ibd_MLE, relatedness="cousin")
  un_Lik<-snpgdsIBDMLELogLik(genofile, ibd_MLE, relatedness="unrelated")
  
  dimnames(fs_Lik)<-list(sample.id, sample.id) 
  dimnames(po_Lik)<-list(sample.id, sample.id) 
  dimnames(hs_Lik)<-list(sample.id, sample.id) 
  dimnames(co_Lik)<-list(sample.id, sample.id) 
  dimnames(un_Lik)<-list(sample.id, sample.id) 
  
  #get coeff dataframes
  
  ibd_MoM.coeff <- snpgdsIBDSelection(ibd_MoM)
  ibd_MoM.coeff$pi_hat<-ibd_MoM.coeff$kinship*2
  
  ibd_MLE.coeff <- snpgdsIBDSelection(ibd_MLE)
  ibd_MLE.coeff$pi_hat<-ibd_MLE.coeff$kinship*2
  
  #bundle results
  
  results<-list(ibd_MoM.coeff, 
                ibd_MLE.coeff, 
                fs_Lik, 
                po_Lik,
                hs_Lik, 
                co_Lik,
                un_Lik)
  
  names(results)<-c(paste0(name_prefix,"ibd_MoM.coeff"),
                    paste0(name_prefix,"ibd_MLE.coeff"),
                    paste0(name_prefix,"fs_Lik"),
                    paste0(name_prefix,"po_Lik"),
                    paste0(name_prefix,"hs_Lik"),
                    paste0(name_prefix,"co_Lik"),
                    paste0(name_prefix,"un_Lik")) 
  
  return(results)
  
}

#Run once more with 275 and 4118 to get complete set

#for transfer

save(subjects, file="subjects.RData")

#ibd.coeff<-ibd.coeff[order(ibd.coeff$kinship, decreasing = TRUE),]

#convert to matrix

#values are right, but IDs are messed up?

g <- snpgdsGetGeno(genofile, sample.id=sample.ids, snp.id=snp.ids, with.id=TRUE)


#dimnames(g)<-list(sample.ids, snp.ids)

#create genofile

snpgdsCreateGeno("test.gds",g$genotype, sample.id=g$sample.id, snp.id=g$snp.id, snpfirstdim=FALSE)

#open genofile
genofile2 <- snpgdsOpen("test.gds")

snpgdsSummary(genofile2)

check <- read.gdsn(index.gdsn(genofile2, "sample.id"))

snpgdsClose(genofile2)
 
ibd_MLE<-snpgdsIBDMLE(genofile2, autosome.only=FALSE, kinship=TRUE, num.thread=nt)

ibd_MLE.coeff <- snpgdsIBDSelection(ibd_MLE)
ibd_MLE.coeff$pi_hat<-ibd_MLE.coeff$kinship*2

View(ibd_MLE.coeff[order(ibd_MLE.coeff$pi_hat, decreasing = TRUE),])
#also remember to close genofiles


#snpgdsClose(genofile)