#pedigree accuracy

#sample coa with different numbers of markers but same set of individuals 
#assume allele frequency is known 

#which individuals to select? #filter out of samples with full individuals

options(stringsAsFactors = FALSE)

full_co<-read.table("C:/ZSL/Coancestry/emp4235/GenotypeData.txt", row.names = 1, sep="")

dolphins<-rownames(full_co)

nmax<-length(dolphins)
nsnps<-ncol(full_co)/2

colIDs<-rep(seq(1,nsnps,1),each=2)

snpset<-c(50,100,200,400,800,1600,3200,nsnps)
nind<-nmax

combns<-expand.grid(snpset, nind)
names(combns)<-c("nsnps", "nind")

#combns<-combns[-nrow(combns),]

#freqs<-read.table("FrequencyData.txt")
#one line for alleles (4, 2)
#one line for frequencies

trio<-readLines("C:/ZSL/Coancestry/emp4235/TrioR.Dat")
#12 lines total 

#change line 1, output path
#change line 2, number of markers
#change line 3, length of number of alleles
#change line 9, number of individuals
#change line 10, length of error estimates
pathpre<-"'C:\\ZSL\\Coancestry\\"
pathsuf<-"\\OutPut.txt'"

#add iterations

#start loop 


for (i in 1:nrow(combns)) {
  
  j<-1
    
    snp<-combns[i,1]
    samp<-combns[i,2]
    
    newf<-paste0("ped", snp, "_", samp, "_", j)
    
    newpath<-paste0(pathpre, newf, pathsuf)
    
    trio[1]<-newpath
    
    trio[2]<-snp
    
    trio[3]<-paste(rep(2, snp), collapse = " ")
    
    trio[9]<-samp
    
    trio[10]<-paste(rep("0.01895", snp), collapse = " ")
    
    dir.create(paste0("C:\\ZSL\\Coancestry\\", newf))
    
    writeLines(trio, con=paste0("C:\\ZSL\\Coancestry\\", newf, "\\TrioR.Dat"))
    
    #may not need frequency data if unknown selected
    
    #generate genotype file
    
    new_snps<-sample(1:nsnps, snp)
    
    new_snps_col<-which(colIDs %in% new_snps)
    
    new_co<-full_co[, new_snps_col]
    
    write.table(new_co, file=paste0("C:\\ZSL\\Coancestry\\", newf, "\\GenotypeData.txt"), 
                row.names=TRUE, col.names = FALSE, sep=" ")

}


setwd("C:/ZSL/Coancestry/")

lf<-list.files()
lf<-lf[grep("ped[[:digit:]]", lf)]

lf<-paste0("cd C:\\ZSL\\Coancestry\\", lf)
runcmd<-rep("..\\trior10", length(lf))

script<-c(rbind(lf, runcmd))

writeLines(script, "batch_coancestry_ped.bat")

#run batch script
#move to new script

library(SocGen)

xfiles<-list.files(path="C:/ZSL/Coancestry", pattern="RelatednessEstimates", recursive = TRUE, full.names = TRUE)
xfiles<-xfiles[grep("ped", xfiles)]
# #check and remove any other filenames

res<-lapply(xfiles, function(x) read.table(x, sep=",", row.names = 1))
names(res)<-xfiles
names(res) <- gsub(pattern = "/RelatednessEstimates.Txt", "", names(res))
names(res) <- gsub(pattern = "C:/ZSL/Coancestry/ped", "", names(res))

ped_est <- list()

#replace smv with kinship data
known_kin<-read.csv("known_kin.csv")

for (i in 1:length(res)) {
  df<-res[[i]][,c(1,2,10)]
  names(df)<-c("ID1", "ID2", "MEst")
  
  nsnps<-strsplit(names(res)[i], "_")[[1]][1]
  # nsnps<-substr(nsnps, 2, nchar(nsnps))
  nind<-strsplit(names(res)[i], "_")[[1]][2]
  df$nind<-nind
  df$nsnps<-nsnps
  
  df<-merge_pairs(known_kin, df, "ID1", "ID2", all.x=TRUE, all.y=FALSE)
  
  ped_est[[i]]<-df}

#Classify false positives and false negatives

#cutoff values
b1<-0.0362
b2<-0.0935
b3<-0.175
b4<-0.4

cats<-matrix(c(0, b1, 0, 
               b1, b2, 0.0625,
               b2, b3, 0.125,
               b3, b4, 0.25,
               b4, 1, 0.5), byrow = TRUE, nrow=5)

#100 snps

cume_rates<-list()

for (i in 1:length(ped_est)){

v5<-ped_est[[i]]

class_rates<-as.data.frame(t(apply(cats, 1, function(x) {
  
  v5$classification<-cut(v5$MEst, x[1:2], include.lowest = TRUE)
  
  tabl<-as.data.frame(table(v5$classification, v5$biparental))
  
  correct<-tabl[tabl$Var2==x[3], "Freq"]
  incorrect<-sum(tabl[!tabl$Var2==x[3], "Freq"])
  
  res<-c(x, correct, incorrect, v5$nsnps[1])
  
  return(res)
})))

names(class_rates)<-c("lower_bound", "upper_bound", "expected", "true_positive",
                      "false_positive", "nsnps")


cume_rates[[i]]<-class_rates
}

cume_rates<-do.call("rbind", cume_rates)
cume_rates<-as.data.frame(apply(cume_rates, 2, as.numeric))

cume_rates$fp_prop<-cume_rates$false_positive/(cume_rates$true_positive+cume_rates$false_positive)

trues<-as.data.frame(table(v5$biparental))

cume_rates$total_trues<-trues$Freq[match(cume_rates$expected, trues$Var1)]

cume_rates$tp_prop<-cume_rates$true_positive/cume_rates$total_trues

windows();plot(tp_prop~nsnps, data=cume_rates, ylim=c(0,1))
points(cume_rates$nsnps, cume_rates$fp_prop, col="red")

write.csv(cume_rates,"kin_classifications.csv")


















