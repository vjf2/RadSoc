#Parentage assignment with sequoia

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

kinship <- expected_kinship(id = "dolphin_id", data = lhl) 

lkin <- kinship[which(kinship$ID1 %in% ID_key$Dolphin_ID &
                        kinship$ID2 %in% ID_key$Dolphin_ID), ]

#load in plink raw file
plink_raw <- "sequoia4p.raw"

maf4loci <- GenoConvert(plink_raw, "raw")

ped_with_sibs <- sequoia(GenoM = maf4loci,
                              LifeHistData = LH_dolphins,
                              Err=0.0161,
                              MaxSibIter = 20, 
                              MaxSibshipSize = 20)

summary(ped_with_sibs)

ped<-ped_with_sibs$Pedigree

SNPd=c(ID_key$Dolphin_ID[which(ID_key$me_paper=="yes")], "JUP", "ABC", "ANE", "ZUL", "NET", "APO", "FCH", "ORY","ATC", "GRV", "MNK", "WON", "HEL", "DUD", "UVE") #adding in UVE causes compare to fail sometimes


compareOUT <- PedCompare(Ped1 = lhl[,c("dolphin_id", "mother_id", "father_id")], Ped2 = ped, na1="",
                         SNPd=SNPd)

consensusped<-compareOUT$ConsensusPed
                         
consensusped$sex<-lhl$sex[match(consensusped$id, lhl$dolphin_id)] 

newkin<-expected_kinship(mother_id = "dam",
                         father_id = "sire", sex = "sex", data=consensusped)

                       
# write.csv(ped, "sequoia_pedigree_output.csv")

parented<-consensusped[which(consensusped$dam!="<NA>" & consensusped$sire!="<NA>" &
                               consensusped$id %in% SNPd),]










