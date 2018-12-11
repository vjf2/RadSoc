#12 subsamples for social model

arg1 <- commandArgs(trailingOnly = TRUE)

arg1<-as.numeric(arg1)
print(arg1)
print(class(arg1))

.libPaths("C:/Users/froug/OneDrive/Documents/R/win-library/3.5")

setwd("C:/Users/froug/Desktop/Real First Chapter")

library(SocGen)

options(stringsAsFactors = FALSE)

# full_co<-read.table(file="coancestry_version.txt", row.names = 1)

dolphins<-read.table("focal_dolphins.txt")[[1]]

#add sexes

ID_key<-read.csv("RawData/dolphin_sample_key.csv")

ID_key$samplingdate<-as.Date(ID_key$samplingdate, format="%m/%d/%Y")
ID_key<-ID_key[order(ID_key$samplingdate),]

females<-unique(ID_key$Dolphin_ID[which(ID_key$Sex=="FEMALE")][which(ID_key$Dolphin_ID[which(ID_key$Sex=="FEMALE")] %in% dolphins)])

males<-unique(ID_key$Dolphin_ID[which(ID_key$Sex=="MALE")][which(ID_key$Dolphin_ID[which(ID_key$Sex=="MALE")] %in% dolphins)])

#scramble
set.seed(arg1)
females<-sample(females)
males<-sample(males)


full_co<-read.table("C:/ZSL/Coancestry/emp4235/GenotypeData.txt", row.names = 1, sep="")

nmax<-length(females) #or dolphins
nsnps<-ncol(full_co)/2

colIDs<-rep(seq(1,nsnps,1),each=2)

snpset<-c(50,100,200,400,800,1600,3200,nsnps)
nind<-c(20, 30, 40, 50, 60, 70, 80)

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

iter<-1

for (i in 1:nrow(combns)) {
  
  j<-arg1
  
  replicate(iter, {
    
    snp<-combns[i,1]
    samp<-combns[i,2]
    
    newf<-paste0("social", snp, "_", samp, "_", j)
    
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
    
    new_subjects<-c(females[1:samp], males[1:samp]) #changed to dolphins
    
    new_snps<-sample(1:nsnps, snp)
    
    new_snps_col<-which(colIDs %in% new_snps)
    
    new_co<-full_co[new_subjects, new_snps_col]
    
    write.table(new_co, file=paste0("C:\\ZSL\\Coancestry\\", newf, "\\GenotypeData.txt"), 
                row.names=TRUE, col.names = FALSE, sep=" ")
    
    j<<-j+1
    
  })
}

#Output genotype, freq, and trio into unique folder

#next time add output counter



#think about how to arrange IDs for social analysis and pedigree analysis

#just do pure random for graph, then do 


#176 dolphins in social analysis #redo social subsamples
#X in pedigree analysis #use allele freqs from subsamples with all indvs

#maybe don't both saving results if quick to read and overwrite?

#write batch script

setwd("C:/ZSL/Coancestry/")

lf<-list.files()
lf<-lf[grep("social[[:digit:]]", lf)]

lf<-paste0("cd C:\\ZSL\\Coancestry\\", lf)
runcmd<-rep("..\\trior10", length(lf))

script<-c(rbind(lf, runcmd))

writeLines(script, "batch_coancestry.bat")
