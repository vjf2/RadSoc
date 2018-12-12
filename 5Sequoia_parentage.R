#Parentage assignment with sequoia

setwd("C:/Users/froug/Desktop/Real First Chapter")

options(stringsAsFactors = FALSE)

library(sequoia)
library(SocGen)

#create life history input 
ID_key <- read.csv("RawData/dolphin_sample_key.csv")
lhl <- read.delim("RawData/LifeHistory20180717.txt",
             header = TRUE,
             sep = "\t")

lhl <- lhl[!duplicated(lhl$dolphin_id), ]

LH_dolphins <- lhl[, c("dolphin_id", "birth_date", "sex")]

LH_dolphins$sex[LH_dolphins$sex == "FEMALE"] <- 1
LH_dolphins$sex[LH_dolphins$sex == "MALE"] <- 2
LH_dolphins$sex[LH_dolphins$sex == "UNKNOWN"] <- NA
LH_dolphins$sex[LH_dolphins$sex == ""] <- NA

LH_dolphins$birthyear <- substr(LH_dolphins$birth_date, 1, 4)

LH_dolphins <- LH_dolphins[, c("dolphin_id", "sex", "birthyear")]

#manual updates

lhl$father_id[which(lhl$dolphin_id == "AGA")] <- "VEE"
lhl$father_id[which(lhl$dolphin_id == "WHP")] <- "TRN"
lhl$mother_id[which(lhl$dolphin_id == "SHO")] <- "SCU"
lhl$mother_id[which(lhl$dolphin_id == "KID")] <- ""

#manual updates pt2

lhl$father_id[which(lhl$father_id=="BOT")]<-""
lhl$father_id[which(lhl$dolphin_id=="PAS")]<-"HII"

#load in plink raw file
plink_raw <- "original2018p.raw"

maf4loci <- GenoConvert(plink_raw, "raw")

set.seed(1)

ped_with_sibs <- sequoia(GenoM = maf4loci,
                              LifeHistData = LH_dolphins,
                              Err=0.01895,
                              MaxMismatch = 30,
                              MaxSibIter = 10, 
                              MaxSibshipSize = 20)

summary(ped_with_sibs)

ped<-ped_with_sibs$Pedigree

# SNPd=c(ID_key$Dolphin_ID[which(ID_key$me_paper=="yes")], "JUP", "ABC", "ANE", "ZUL", "NET", "APO", "FCH", "ORY","ATC", "GRV", "MNK", "WON", "HEL", "DUD", "UVE") #adding in UVE causes compare to fail sometimes

SNPd=c(rownames(maf4loci))

compareOUT <- PedCompare(Ped1 = lhl[,c("dolphin_id", "mother_id", "father_id")], 
                         Ped2 = ped, na1="", SNPd=SNPd)

consensusped<-compareOUT$ConsensusPed
                         
consensusped$sex<-lhl$sex[match(consensusped$id, lhl$dolphin_id)] 

# newkin<-expected_kinship(mother_id = "dam",
#                          father_id = "sire", sex = "sex", data=consensusped)

                       
# write.csv(ped, "sequoia_pedigree_output.csv")

parented<-consensusped[which(consensusped$dam!="<NA>" & consensusped$sire!="<NA>" &
                               consensusped$id %in% SNPd),]


#missing mums
cc<-compareOUT$P1only[which(compareOUT$P1only$id!="<NA>"
                            & compareOUT$P1only$dam.1 %in% SNPd
                            & compareOUT$P1only$Parent=="dam"),]

nrow(cc)
nmoms<-ped[ped$id %in% SNPd & !is.na(ped$dam),]
1-nrow(cc)/nrow(nmoms[grep("F0", nmoms$dam, invert = TRUE),])

compareOUT$P2only



