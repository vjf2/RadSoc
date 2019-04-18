options(stringsAsFactors = FALSE)

library(SocGen)

library(adegenet)

library("hierfstat")

##adegenet uses two new S4 classes to store genetic data, "genind" for smaller (i.e. microsat) data and "genlight" for storing large, genome-wide SNP datasets

getClassDef("genind")
getClassDef("genlight")

dol<-read.PLINK("original2018p.raw", parallel=FALSE)


x.mat <- as.matrix(dol) # x is a genlight object
x.mat[x.mat == 0] <- "1/1" # homozygote reference
x.mat[x.mat == 1] <- "1/2" # heterozygote
x.mat[x.mat == 2] <- "2/2" # homozygote alternate
x.gid <- df2genind(x.mat, sep = "/", ploidy = 2)

x<-genind2hierfstat(x.gid,pop=rep(1,272))

sum_x <- basic.stats(x)

Ho = 0.2407383

He = 0.2433753

sd(sum_x$Hs[,1])
sd(sum_x$Ho[,1])


#Van Cise uses the strataG package


#install.packages("strataG")

