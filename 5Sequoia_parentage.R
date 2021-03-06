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

allloci <- GenoConvert(plink_raw, "raw")

set.seed(1)

ped_with_sibs <- sequoia(GenoM = allloci,
                              LifeHistData = LH_dolphins,
                              Err=0.01895,
                              Tassign = 0.5, # default 0.5
                              Tfilter = -2, #default -2
                              MaxSibIter = 10, 
                              MaxSibshipSize = 20)

summary(ped_with_sibs)

ped<-ped_with_sibs$Pedigree

SNPd=c(rownames(allloci))

compareOUT <- PedCompare(Ped1 = lhl[,c("dolphin_id", "mother_id", "father_id")], 
                         Ped2 = ped, na1="", SNPd=SNPd)

consensusped<-compareOUT$ConsensusPed
                         
consensusped$sex<-lhl$sex[match(consensusped$id, lhl$dolphin_id)] 

# write.csv(consensusped, "consensusped.csv")

#missing mums
cc<-compareOUT$P1only[which(compareOUT$P1only$id!="<NA>"
                            & compareOUT$P1only$dam.1 %in% SNPd
                            & compareOUT$P1only$Parent=="dam"),]

nrow(cc)
nmoms<-ped[ped$id %in% SNPd & !is.na(ped$dam),]
1-nrow(cc)/nrow(nmoms[grep("F0", nmoms$dam, invert = TRUE),])

compareOUT$P2only

ConfPr<-EstConf(Ped = ped, LifeHistData = LH_dolphins, Specs = ped_with_sibs$Specs, Full = TRUE,
        nSim = 20, ParMis = 0.6, args.sim = NULL, return.PC = TRUE,
        quiet = TRUE)

save(ConfPr, file="confrpr.RData")


#create list of IDs that don't have a parent in the sample

afreqid<-ped_with_sibs$Pedigree$id[which(is.na(ped_with_sibs$Pedigree$dam) & is.na(ped_with_sibs$Pedigree$sire))]

afreqid<-afreqid[1:161] #Remove BMT and dummies

#what about individuals with siblings in the data?

write.table(afreqid, sep="\t",file="afreqid.txt", row.names = FALSE, col.names = FALSE)


#All individuals with a known mom who isn't in the sample
compareOUT$MergedPed[which(is.na(compareOUT$MergedPed$dam.2) & 
                             !is.na(compareOUT$MergedPed$dam.1) &
                             !is.na(compareOUT$MergedPed$id)),]

mom_counts<-table(compareOUT$MergedPed$dam.1[which(is.na(compareOUT$MergedPed$dam.2) & 
                                   !is.na(compareOUT$MergedPed$dam.1) &
                                   !is.na(compareOUT$MergedPed$id))]
)

sibmoms<-names(mom_counts)[which(mom_counts>1)]

check_imputed<-compareOUT$MergedPed[which(compareOUT$MergedPed$dam.1 %in% sibmoms &
                             !is.na(compareOUT$MergedPed$id)),]

check_imputed[order(check_imputed$dam.1),]


