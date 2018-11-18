# setwd("..")

options(stringsAsFactors = FALSE)

raw<-read.delim("original2018p.ped", sep="", header=FALSE, colClasses = 
                  "character")

check<-as.matrix(raw[,c(7:ncol(raw))])

check[check=="A"]<-1
check[check=="T"]<-2
check[check=="G"]<-3
check[check=="C"]<-4

check<-apply(check, 2, as.numeric)

check<-cbind(raw[,2],as.data.frame(check))

write.table(check, sep="\t", file="coancestry_version.txt", row.names = FALSE, col.names = FALSE)

#use coancestry empirical to generate allele frequencies

#error file

tc<-t(check[,-1])

#zeros are missing data

cs<-colSums(check==0)

missing<-cs/280

odds<-seq(2,length(missing),2)

missing_rates<-missing[odds]

error_rates<-matrix(rbind(missing_rates, 
                    rep(0.0161, length(missing_rates)), 
                        rep(0, length(missing_rates))), nrow=3)

write.table(t(error_rates), sep=" ", file="coancestry_error_rates.txt", row.names = FALSE, col.names = FALSE)

write.table(t(error_rates[2,]), file="coancestry_just_error.txt", row.names = FALSE, col.names = FALSE)

#matrix of types of relationships

n=100 #number to simulate

#deltas 7-9

family_mat <- matrix(c(0, 1, 0, #parent-offspring
                .25, .5, .25, #full siblings
                0, .5, .5, #half-siblings, avuncular, grandparents
                0, .25, 0.75, #first cousins
                0, 0.125, 0.875, #half-first cousins
                0, 0, 1 #unrelated
                ), nrow=6, byrow=TRUE)

family_mat<-cbind(matrix(rep(0, 36), nrow=6),family_mat, rep(n, 6))

write.table(family_mat, sep="\t", file="true_rel2sim.txt", row.names = FALSE, col.names = FALSE)

#add into coancestry to simulate coefficients
#takes about 5 min to run

rels <- read.table("C:/ZSL/Coancestry/sim_maf05/RelatednessEstimates.txt", sep=",", row.names=1)

names(rels)<-c("Ind1","Ind2","RelCat", "TrioEst","WEst","LLEst","LREst","REst","QGEst","MEst")

rels<-rels[,c(1:10)]

#Rename relcat

rels$RelCat <- substr(rels$Ind1, 1, 4)

rels$RelCat[rels$RelCat=="R001"]<-"PO"
rels$RelCat[rels$RelCat=="R002"]<-"FS"
rels$RelCat[rels$RelCat=="R003"]<-"HS"
rels$RelCat[rels$RelCat=="R004"]<-"FC"
rels$RelCat[rels$RelCat=="R005"]<-"HC"
rels$RelCat[rels$RelCat=="R006"]<-"UR"

#melt

m<-reshape2::melt(rels, id.vars=c("Ind1", "Ind2", "RelCat"), na.rm=TRUE)

m$value<-as.numeric(m$value)

#check genotype data - do we really have mistypings and missing data


#don't forget to remove either wang or lynchli at end

################
#4 plots

windows()
# pdf(file="allcats2.pdf", width=7)
# par(mfrow=c(4,1),mar=c(2.6,10.8,2.1,2.1), mgp=c(4,1,0))

par(mfrow=c(6,1), oma = (c(6,6,0,0) + 0.1), mar = (c(0.5,0,1.5,1) + 0.1))

boxplot(value~variable, data=m[m$RelCat=="PO",], col=rainbow(7),
        ylab="", 
        xaxt="n", 
        yaxt="n", 
        cex.lab=1.3,
        boxwex=0.5,
        # ylim=c(0.33,0.6),
        main="Parent-Offspring")

axis(2,las=1, cex.axis=1.3)

abline(h=0.5, lty=2, col="darkgrey")

boxplot(value~variable, data=m[m$RelCat=="FS",], col=rainbow(7),
        ylab="", 
        xaxt="n", 
        yaxt="n", 
        cex.lab=1.3,
        boxwex=0.5,
        # ylim=c(0.33,0.6),
        main="Full Siblings", cex.main=1.3)

axis(2,las=1, cex.axis=1.3)

abline(h=0.5, lty=2, col="darkgrey")

boxplot(value~variable, data=m[m$RelCat=="HS",], col=rainbow(7),
        ylab="", 
        xaxt="n", 
        yaxt="n", 
        cex.lab=1.3,
        boxwex=0.5,
        main="Half Siblings")

axis(2,las=1, cex.axis=1.3)

abline(h=0.25, lty=2, col="darkgrey")

boxplot(value~variable, data=m[m$RelCat=="FC",], col=rainbow(7),
        ylab="", 
        xaxt="n", 
        yaxt="n", 
        cex.lab=1.3,
        boxwex=0.5,
        main="Avuncular")

axis(2,las=1, cex.axis=1.3)

abline(h=0.125, lty=2, col="darkgrey")

boxplot(value~variable, data=m[m$RelCat=="HC",], col=rainbow(7),
        ylab="", 
        xaxt="n", 
        yaxt="n", 
        cex.lab=1.3,
        boxwex=0.5,
        main="Half Cousins")

axis(2,las=1, cex.axis=1.3)

abline(h=0.0625, lty=2, col="darkgrey")

boxplot(value~variable, data=m[m$RelCat=="UR",], col=rainbow(7),
        ylab="", 
        xaxt="n", 
        yaxt="n", 
        cex.lab=1.3,
        boxwex=0.5,
        main="Unrelated")

axis(2,las=1, cex.axis=1.3)

axis(1, at=1:7, labels=c("TrioEst","WEst","LLEst","LREst","REst","QGEst","MEst"), 
     tick=TRUE, cex.axis=1.3, padj = 0.5, las=1)

abline(h=0, lty=2, col="darkgrey")

title(xlab = "Estimator",
      ylab = "Relatedness Coefficient",
      outer = TRUE, line = 4, cex.lab=1.5)


dev.off()

#Different (and worse?) results when account for mistypings is included
#Check with real data and simulated

#get cor and rmse 

rels$expected[rels$RelCat=="PO"]<-0.5
rels$expected[rels$RelCat=="FS"]<-0.5
rels$expected[rels$RelCat=="HS"]<-0.25
rels$expected[rels$RelCat=="FC"]<-0.125
rels$expected[rels$RelCat=="HC"]<-0.0625
rels$expected[rels$RelCat=="UR"]<-0

rels$expected<-as.numeric(rels$expected)

rmse_po<-apply(rels[rels$RelCat=="PO",c(4:10)], 2, function(x) {
  sqrt((1/nrow(rels[rels$RelCat=="PO",]))*sum(abs(x-rels[rels$RelCat=="PO","expected"])^2))})

rmse_fs<-apply(rels[rels$RelCat=="FS",c(4:10)], 2, function(x) {
  sqrt((1/nrow(rels[rels$RelCat=="FS",]))*sum(abs(x-rels[rels$RelCat=="FS","expected"])^2))})

rmse_hs<-apply(rels[rels$RelCat=="HS",c(4:10)], 2, function(x) {
  sqrt((1/nrow(rels[rels$RelCat=="HS",]))*sum(abs(x-rels[rels$RelCat=="HS","expected"])^2))})

rmse_fc<-apply(rels[rels$RelCat=="FC",c(4:10)], 2, function(x) {
  sqrt((1/nrow(rels[rels$RelCat=="FC",]))*sum(abs(x-rels[rels$RelCat=="FC","expected"])^2))})

rmse_hc<-apply(rels[rels$RelCat=="HC",c(4:10)], 2, function(x) {
  sqrt((1/nrow(rels[rels$RelCat=="HC",]))*sum(abs(x-rels[rels$RelCat=="HC","expected"])^2))})

rmse_ur<-apply(rels[rels$RelCat=="UR",c(4:10)], 2, function(x) {
  sqrt((1/nrow(rels[rels$RelCat=="UR",]))*sum(abs(x-rels[rels$RelCat=="UR","expected"])^2))})

rmse_all<-apply(rels[,c(4:10)], 2, function(x) {
  sqrt((1/nrow(rels))*sum(abs(x-rels$expected)^2))})


cor_all<-apply(rels[,c(4:10)], 2, function(x) {cor(x, rels$expected)})

rmse_comb<-rbind(rmse_po, rmse_fs, rmse_hs, rmse_fc, rmse_hc, rmse_ur, rmse_all, cor_all)

# write.csv(rmse_comb, "rmse_account_error.csv")
