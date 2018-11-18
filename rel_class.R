#check out misclassificaitons in empirical data

options(stringsAsFactors = FALSE)

r<-read.csv("C:/Users/froug/Desktop/Chapter 1 - Relationships with Kin/RawData/max_likelihood_relatedness.csv")

r<-r[,c("xID1", "ID2", "relatedness")]

#add self-pairs into relatedness

r<-rbind(r, data.frame(xID1=consensusped$id, ID2=consensusped$id, relatedness=rep(1, length(consensusped$id))))

# bothpar<-consensusped[which(consensusped$dam.cat=="GG" & consensusped$sire.cat=="GG"
#                              | consensusped$id=="FCH" | consensusped$id=="CLV"),]

#CLV and FCH need to be added back in

bothpar<-consensusped[complete.cases(consensusped[,1:3]),]

#can add dummy parents in 


pedids<-bothpar$id

r1<-r[intersect(which(r$xID1 %in% pedids),  which(r$ID2 %in% pedids)),]

newkin<-expected_kinship(mother_id = "dam",
                         father_id = "sire", sex = "sex", data=consensusped)

re<-merge_pairs(r1, newkin[,c(1,2,5)], "xID1", "ID2", "ID1", "ID2", all.x=TRUE, all.y=FALSE)

re<-merge(re, bothpar[,1:3], by.x="xID1", by.y="id", all.x=TRUE, all.y=FALSE)

re<-merge(re, bothpar[,1:3], by.x="ID2", by.y="id", all.x=TRUE, all.y=FALSE)

#M1M2
re<-merge_pairs(re, r, "dam.x", "dam.y", "xID1", "ID2")
names(re)[ncol(re)]<-"M1M2"

#M1F2
re<-merge_pairs(re, r, "dam.x", "sire.y", "xID1", "ID2")
names(re)[ncol(re)]<-"M1F2"

#F1M2
re<-merge_pairs(re, r, "sire.x", "dam.y", "xID1", "ID2")
names(re)[ncol(re)]<-"F1M2"

#F1F2
re<-merge_pairs(re, r, "sire.x", "sire.y", "xID1", "ID2")
names(re)[ncol(re)]<-"F1F2"

write.csv(re, "../classifications.csv")

table(re$biparental)

#add in more known relatives

#add in all individs with two known parents, and all parent-offspring
#dummy parents

#Use 95/97% of simulated data to assign relatives 

#Get better error rate


