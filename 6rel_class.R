#check out misclassificaitons in empirical data

#run sequoia_parentage first

library(SocGen)

options(stringsAsFactors = FALSE)

reals<-read.table("C:/ZSL/Coancestry/emp4235/RelatednessEstimates.Txt", sep=",", row.names = 1)

names(reals)<-c("ID1", "ID2", "Pair", "TrioEst","WEst","LLEst","LREst","REst","QGEst","MEst")

r<-reals[,c("ID1", "ID2", "MEst")]

names(r)[3]<-"relatedness"

SNPd<-unique(c(reals$ID1, reals$ID2))

#add self-pairs into relatedness

consensusped<-read.csv("consensusped.csv", row.names=1)

r<-rbind(r, data.frame(ID1=consensusped$id, ID2=consensusped$id, relatedness=rep(1, length(consensusped$id))))

parented<-consensusped[which(consensusped$dam!="<NA>" & consensusped$sire!="<NA>" &
                               consensusped$id %in% SNPd),]

bothpar<-parented[,1:3]

#can add dummy parents in 

pedids<-bothpar$id

#try without crazy rules



#rm JMY, unsupported
# pedids<-setdiff(pedids, "JMY")

r1<-r[intersect(which(r$ID1 %in% pedids),  which(r$ID2 %in% pedids)),]

newkin<-expected_kinship(mother_id = "dam",
                         father_id = "sire", sex="sex", data=consensusped)

re<-merge_pairs(r1, newkin[,c(1,2,5)], "ID1", "ID2", "ID1", "ID2", all.x=TRUE, all.y=FALSE)

re<-merge(re, bothpar[,1:3], by.x="ID1", by.y="id", all.x=TRUE, all.y=FALSE)

re<-merge(re, bothpar[,1:3], by.x="ID2", by.y="id", all.x=TRUE, all.y=FALSE)

#M1M2
re<-merge_pairs(re, r, "dam.x", "dam.y", "ID1", "ID2")
names(re)[ncol(re)]<-"M1M2"

#M1F2
re<-merge_pairs(re, r, "dam.x", "sire.y", "ID1", "ID2")
names(re)[ncol(re)]<-"M1F2"

#F1M2
re<-merge_pairs(re, r, "sire.x", "dam.y", "ID1", "ID2")
names(re)[ncol(re)]<-"F1M2"

#F1F2
re<-merge_pairs(re, r, "sire.x", "sire.y", "ID1", "ID2")
names(re)[ncol(re)]<-"F1F2"

write.csv(re, "classifications.csv")

table(re$biparental)

aggregate(relatedness.x~biparental, data=re, mean)

#remove unknowns that don't have all parents typed
#remove ones were parents are related
#unsupported old data
#add parents back-in 

#rules for dyad inclusion
#unrelated pairs need to have all 4 parents genotyped and unrelated
#all r =0.5
#all r = 0.25 if all 4 parents in pedigree
# r = 0.0625 and 0.125 if all 4 parents known and at least 3 genotyped in sample

re$countna<-apply(re[,9:12], 1, function(x) sum(is.na(x)))

re$rmrel<-ifelse(re$biparental=="0" & re$countna>0, "yes", "")

re$relmax<-apply(re[,9:12], 1, max)

#second highest value
re$relmax2<-apply(re[,9:12], 1, function(x) sort(x, na.last=FALSE)[4-1])

re$rmrel<-ifelse(re$biparental=="0" & re$relmax>0.0362, "yes", re$rmrel)

re$rmrel<-ifelse(re$biparental %in% c(0.0625, 0.125) & re$countna>2, "yes", re$rmrel)


ref<-re[which(re$rmrel!="yes"),]

table(ref$biparental)

#add in known parents

re2<-merge_pairs(r, newkin[,c(1,2,5)], "ID1", "ID2", "ID1", "ID2", all.x=TRUE, all.y=FALSE)

repar<-re2[which(re2$biparental==0.5),]

names(repar)[3]<-"relatedness.x"

ref<-rbind(ref[,1:4], repar)

ref<-reduce_pairs(ref, "ID1", "ID2")

write.csv(ref, "known_kin.csv", row.names = FALSE)
