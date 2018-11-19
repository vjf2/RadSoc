#V4
#Vivienne Foroughirad
#Created June 13, 2017
#Modified November 3, 2017
#Modified September 5, 2018

#Read in raw data from SNPs and filter based on MAF, CallRate, Errors, Linkage, HWE, determine sexes
#Format for Plink

#Read in raw datasets for single and double file, and ID key
setwd("..")

options(stringsAsFactors = FALSE)

single_raw_old<-read.csv("RawData/OrderAppendix_3_DTur17-2726.csv", colClasses = "character", skip=6)

single_raw_new<-read.csv("RawData/Report_DTur17-3054_2_moreOrders_SNP_singlerow_2.csv", colClasses = "character", skip=6, na.strings="-")

double_raw<-read.csv("RawData/Report_DTur17-3054_2_moreOrders_SNP_2.csv", colClasses = "character", skip=6)

ID_key<-read.csv("RawData/dolphin_sample_key.csv", colClasses = "character")

alleleunion<-merge(single_raw_new, single_raw_old[,c(1,4:7)], by="AlleleID", all=FALSE)

single_raw<-alleleunion[,c(1:22, ncol(alleleunion):(ncol(alleleunion)-3),23:(ncol(alleleunion)-4))]

#Extract only samples needed for analysis

dolphins_to_include<-ID_key$Dolphin_ID[which(ID_key$me_paper=="yes" | ID_key$me_paper=="yes_v2")]

#need to match back to retain duplicates

samples_to_include<-ID_key$Sample_ID[which(ID_key$Dolphin_ID %in% dolphins_to_include)]

single_raw<-single_raw[,c(1:26, which(names(single_raw) %in% samples_to_include))]

double_raw<-double_raw[,c(1:22, which(names(double_raw) %in% samples_to_include))]

#set the min call rate

mincallrate<-0.95

#count numner of "-" per individual and per row

missing_per_sample<-apply(single_raw[,27:ncol(single_raw)],2, 
                          function(x) {dash<-sum(is.na(x))
                                        return(dash)})

percent_missing<-missing_per_sample/nrow(single_raw) #total snps

#call_rate

missing_per_snp<-apply(single_raw[,27:ncol(single_raw)],1, 
                          function(x) {dash<-sum(is.na(x))
                          return(dash)})

percent_missing_per_snp<-missing_per_snp/length(samples_to_include) #total samples (with dups)

#Raw stats - average number of reads

TotalAvg<-as.numeric(single_raw_new$AvgCountRef)+as.numeric(single_raw_new$AvgCountSnp)
summary(TotalAvg)

single_raw$TotalAvg<-as.numeric(single_raw$AvgCountRef)+as.numeric(single_raw$AvgCountSnp)
summary(single_raw$TotalAvg)

#Number of sequences that map to genome

1-length(single_raw$Chrom_Tursiops_v1_2016[which(single_raw$Chrom_Tursiops_v1_2016=="")])/nrow(single_raw) #plus filtered out sequences n=9928

#Max E-value

max(as.numeric(single_raw$AlnEvalue_Tursiops_v1_2016)[which(as.numeric(single_raw$AlnEvalue_Tursiops_v1_2016)<999)])

#Get sexes for all the individuals

#sex locuses

sex_snps<-c("15024211", "15025192", "15035552", "15035553", "21678180", "15035836")

#snp 15025192 gets filtered out because it doesn't map to genome

x<-unlist(sapply(sex_snps, function(x) grep(x, single_raw$AlleleID)))

sexes<-single_raw[x, 27:ncol(single_raw)]

sexes<-apply(sexes,2,as.numeric)

sexes<-t(rbind(sexes, colSums(sexes)))

# write.csv(single_raw[x,], file="sex_loci.csv")

#add sexes to ID key

ID_key<-merge(ID_key, sexes, by.x="Sample_ID", by.y=0, all.x=TRUE, all.y=FALSE)

#note - check for genotyping errors before calling sexes

#remove sex snps

single_raw<-single_raw[-x,]

#Check percent identity between duplicate samples, and remove snps from list that typed differently

#Identify duplicate samples
dups<-c("BLI", "SCA", "SOG", "FSB", "LAS", "STA", "JUP") #remove SHO and PHO bad samples
error_rate<-list()

#This calculates the error rates in duplicate samples

for (i in 1:length(dups)){

cdup<-ID_key$Sample_ID[ID_key$Dolphin_ID==dups[i]]
ctitle<-paste0("remove",i)

single_raw[,ctitle]<-ifelse(single_raw[,cdup[1]]==single_raw[,cdup[2]], "agree", "disagree")
error_rate[[i]]<-table(single_raw[,ctitle])[2]/sum(table(single_raw[,ctitle]))
}

ers<-unlist(error_rate)

mean(ers) #excluding missing data

#Sample error, remove duplicates

cdups<-ID_key$Sample_ID[ID_key$Dolphin_ID %in% dups]

#remove duplicate sample with fewer snps typed

#count snps typed

choice_table<-single_raw[,cdups]
nsnps_typed<-t(apply(choice_table, 2, table))

dupID<-ID_key[ID_key$Dolphin_ID %in% dups,1:2]

#match up

sample_choice<-merge(nsnps_typed, dupID, by.x=0, by.y="Sample_ID")

sample_choice$totmarkers<-sample_choice[,2]+sample_choice[,3]+sample_choice[,4]

sample_choice<-sample_choice[order(sample_choice$totmarkers),]

minus<-sample_choice$Row.names[!duplicated(sample_choice$Dolphin_ID)]

ID_key$Dolphin_ID[match(minus,ID_key$Sample_ID)]

single<-single_raw[,!names(single_raw) %in% minus]
double<-double_raw[,!names(double_raw) %in% minus]

#actually remove the disagreeing ones (could change the mismatch criterion to 2 in future)

single$dissum<-apply(single[,grep("remove",names(single))], 1,function(w) {ds<-sum(w=="disagree", na.rm=TRUE); return(ds)})

#calculate new error rate from filtered data at the end
#write.csv(single, file="for_error_rate.csv")

single<-single[which(single$dissum<=1),]
                 
#remove the removal columns and total avg
single<-single[,-(grep("remove",names(single)))]

#Note - skip HWE and MAF check here, do in plink instead

#Remove any samples with more than 10% missing (added to ID_key so don't have to redo)
# 
# g10<-percent_missing[which(percent_missing>0.10)]
# 
# single<-single[,setdiff(names(single), names(g10))]

#Filter on average number of reads

single$TotalAvg<-as.numeric(single$AvgCountRef)+as.numeric(single$AvgCountSnp)

single<-single[which(single$TotalAvg>=10),]

# dolphins_to_include<-setdiff(dolphins_to_include, ID_key$Dolphin_ID[which(ID_key$Sample_ID %in% names(g10))])

#Remove on CallRate and maximum heterozygosity filter

single$percent_missing_per_snp<-apply(single[,which(names(single) %in% samples_to_include)],1, 
                       function(x) {dash<-sum(is.na(x))
                       return(dash)})

single$percent_missing_per_snp<-single$percent_missing_per_snp/length(dolphins_to_include)

max_hetero=0.60

#redo min call rate

single<-subset(single, percent_missing_per_snp<=(1-mincallrate)); nrow(single)
single<-subset(single, FreqHets<=max_hetero); nrow(single)

#remove unmapped snps

single<-single[which(single$Chrom_Tursiops_v1_2016!=""),]; nrow(single)

#Organize SNPS per scaffold and export to plink for linkage calculation

#add scaffold name to map file

#format Allele ID column for matching by creating unique ID

AlleleNumber<-strsplit(single$AlleleID,"\\|")
AlleleNumber<-unlist(lapply(AlleleNumber, "[[",1))

single<-cbind(single, AlleleNumber)

AlleleNumberD<-strsplit(double$AlleleID,"\\|")
AlleleNumberD<-unlist(lapply(AlleleNumberD, "[[",1))

double<-cbind(double, AlleleNumberD)

#filter double

double<-double[which(double$AlleleNumberD %in% single$AlleleNumber),]

final_double<-subset(double, double$AlleleNumberD %in% single$AlleleNumber) #just to check

#still duplicates because of multiple snps in the same read
#give double a dummy ID

n<-dim(double)[1]/2

double$allele_unique<-rep(1:n, each=2)

final_ids<-subset(double$allele_unique, double$AlleleID %in% single$AlleleID)

final_double<-subset(double, double$allele_unique %in% final_ids)


#Reformat final double for plink

cervus_double<-t(final_double)

n2<-dim(cervus_double)[2]/2

al<-rep(c("a","b"), n2)

cervus_double<-as.data.frame(rbind(cervus_double, al))

tt<-paste0(cervus_double["AlleleNumberD",], cervus_double["al",])

names(cervus_double)<-tt

#plink is letter, franz is 1-4 based on letter

plink_double<-t(cervus_double)

#Make nuc1 and nuc0 columns

lt<-as.data.frame(plink_double[,c("SNP","allele_unique")], na.strings=c("", NA))

names(lt)<-c("SNP1", "allele_unique")

lt$SNP1[lt$SNP1==""] <- NA

lt<-lt[complete.cases(lt),]

plink_double<-merge(plink_double, lt, by="allele_unique")

AlleleLetter<-strsplit(plink_double$SNP1,":")
plink_double$AlleleLetter<-unlist(lapply(AlleleLetter, "[[",2))

RefLetter<-strsplit(plink_double$AlleleLetter, ">")
plink_double$RefLetter<-unlist(lapply(RefLetter, "[[",1))
plink_double$AltLetter<-unlist(lapply(RefLetter, "[[",2))

plink_double$Nuc0<-ifelse(plink_double$al=="b", plink_double$AltLetter, plink_double$RefLetter)
plink_double$Nuc1<-ifelse(plink_double$al=="a", plink_double$AltLetter, plink_double$RefLetter)

x<-plink_double

#get start and end of individuals
first<-which(names(x)=="RepAvg")+1
last<-which(names(x)=="AlleleNumberD")-1

for (i in first:last) {
  
  x[x[,i]=="0", names(x)[i]]<-x[x[,i]=="0", "Nuc0"]
  x[x[,i]=="1", names(x)[i]]<-x[x[,i]=="1", "Nuc1"]
  
}

x[,first:last]<-apply(x[,first:last], 2, function(x) gsub("-", 0, x))

plink1<-x

#add new scaffold and position for map file 

plink1$newscaffold<-single$Chrom_Tursiops_v1_2016[match(plink1$AlleleNumberD,single$AlleleNumber)]
plink1$newposition<-single$ChromPos_Tursiops_v1_2016[match(plink1$AlleleNumberD,single$AlleleNumber)]

#set up ped and map file with allele names that can be matched back to raw data (allele_unique)

plink<-plink1[,c(1, first:ncol(plink1))] 

prefix<-rep("G", nrow(plink))
suffix<-rep(c(".1",".2"), nrow(plink)/2)

map_loci<-paste0("G",unique(plink$allele_unique))
map_loci<-gsub(" ", "", map_loci, fixed=TRUE)

plink$name<-paste0(prefix,plink$allele_unique, suffix)
plink$name<-gsub(" ", "", plink$name, fixed=TRUE)

plink<-t(plink)

colnames(plink)<-plink[nrow(plink),]

x<-merge(plink, ID_key[,1:2], by.x="row.names", by.y="Sample_ID")

x[,1]<-x$Dolphin_ID
pl<-x[,-ncol(x)]

pl1<-data.frame(dummy_group=rep(1,nrow(pl)),pl)

m<-matrix(rep(0, 4*nrow(pl)), nrow=nrow(pl))

pl2<-data.frame(pl1[,1:2],m, pl1[,3:ncol(pl1)]) 

original.ped<-pl2

#set up ped and map file with allele names that can be matched back to raw data (allele_unique)

#fix map locus names 

nsnps<-length(map_loci)

original.map<-data.frame(rep(0,nsnps), map_loci, rep(0, nsnps), rep(0, nsnps))

original.map[,1]<-plink1$newscaffold[seq(1,nrow(plink1)-1,2)]
original.map[,4]<-plink1$newposition[seq(1,nrow(plink1)-1,2)]

original.map[,1]<-gsub(".1_scaffold","",original.map[,1])

write.table(original.ped, "original2018.ped", sep="\t", col.names=FALSE, row.names=FALSE,quote=FALSE)
write.table(original.map, "original2018.map", sep="\t", col.names=FALSE, row.names=FALSE,quote=FALSE)

#run plink_commands.bat

#calculate new error rate

fer<-read.csv("for_error_rate.csv")

filt<-read.table("original2018p.map", header=FALSE, sep="\t")
filt$V1<-paste0(filt$V1, ".1_scaffold")

finalerr<-merge(fer, filt, by.x=c("Chrom_Tursiops_v1_2016", "ChromPos_Tursiops_v1_2016"), 
            by.y=c("V1", "V4"), all.y=TRUE)

ctitle<-grep("remove", names(finalerr))

error_rate<-lapply(ctitle, function(x) table(finalerr[,x])[2]/sum(table(finalerr[,x])))

ers<-unlist(error_rate)
mean(ers)







