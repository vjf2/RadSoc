#Generate Coancestry subsamples

#setwd("..")

library(SocGen)

options(stringsAsFactors = FALSE)

# full_co<-read.table(file="coancestry_version.txt", row.names = 1)

full_co<-read.table("C:/ZSL/Coancestry/emp4235/GenotypeData.txt", row.names = 1, sep="")

nmax<-nrow(full_co)
nsnps<-ncol(full_co)/2

colIDs<-rep(seq(1,nsnps,1),each=2)

snpset<-c(50,100,200,400,800,1600,3200,nsnps)
nind<-c(10,20,30,40,50,100,150, nmax)

combns<-expand.grid(snpset, nind)
names(combns)<-c("nsnps", "nind")

combns<-combns[-nrow(combns),]

freqs<-read.table("FrequencyData.txt")
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

newf<-"test"

newpath<-paste0(pathpre, newf, pathsuf)

trio[1]<-newpath

trio[2]<-400

trio[3]<-substr(trio[3],1,(400*2))

trio[9]<-50

trio[10]<-paste(rep("0.01895 ", 400), collapse = " ")



# dir.create(paste0("C:\\ZSL\\Coancestry\\", newf))

writeLines(trio, con=paste0("C:\\ZSL\\Coancestry\\", newf, "\\TrioR.Dat"))

#may not need frequency data if unknown selected

#generate genotype file

new_subjects<-sample(rownames(full_co), 50)

new_snps<-sample(1:nsnps, 400)

new_snps_col<-which(colIDs %in% new_snps)

new_co<-full_co[new_subjects, new_snps_col]

write.table(new_co, file=paste0("C:\\ZSL\\Coancestry\\", newf, "\\GenotypeData.txt"), 
            row.names=TRUE, col.names = FALSE, sep=" ")


#Output genotype, freq, and trio into unique folder




#think about how to arrange IDs for social analysis and pedigree analysis

#just do pure random for graph, then do 


#176 dolphins in social analysis #redo social subsamples
#X in pedigree analysis #use allele freqs from subsamples with all indvs

#maybe don't both saving results if quick to read and overwrite?




