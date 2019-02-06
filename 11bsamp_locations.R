#Sampling Location and Centroid Map

options(stringsAsFactors = FALSE)

library(raster)
library(rgdal)
library(rgeos)

#load in a geotiff made in ArcMap (using composite bands tool)

lgt<-brick("C:/Users/froug/Desktop/Galezo Juvenile Paper/ls81200.tif")

#get centroids and get sampling locations

ID_key<-read.csv("RawData/dolphin_sample_key.csv") 

sl_dol<-ID_key[which((ID_key$me_paper=="yes" | ID_key$me_paper=="yes_v2") & ID_key$Npoints<35),] ###check to make sure corect

cent_dol<-ID_key[which(ID_key$me_paper=="yes" | ID_key$me_paper=="yes_v2" & ID_key$Npoints>=35),] ###check to make sure correct

sl_xydata<-as.matrix(sl_dol[,c("Sampling_location_east","Sampling_location_south")])

#get centroids from social attributes code

kernelcontours90 <- getverticeshr(re, percent=90, standardize = TRUE)

cents <- gCentroid(kernelcontours90,byid=TRUE) 

proj4string(cents)<-CRS("+proj=tmerc +lat_0=-25 +lon_0=113 +k=0.99999 +x_0=50000 +y_0=100000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")


cents<-spTransform(cents, CRS("+proj=utm +zone=49 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))

xydata<-as.matrix(as.data.frame(cents))

#generate good colors
pcol <- rgb(0, 20, 50, max = 255, alpha = 255, names = "background")
ccol <- rgb(255, 140, 0, max = 255, alpha = 255, names = "circle")
tcol <- rgb(173, 255, 47, max = 255, alpha = 255, names = "triangle")
bcol <- rgb(255, 255, 255, max = 255, alpha = 255, names = "border")
  
#test colors
windows();plot(rnorm(100), type="n")
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = pcol)
points(rnorm(100), pch=24, cex=0.8, col=tcol, bg=tcol)
points(rnorm(100), pch=21, cex=1, col=ccol, bg=ccol)  

#convert points to coordinate system of landsat raster
sl_xydata<-as.data.frame(project(sl_xydata, "+proj=utm +zone=49 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))

####PLOTTING

tiff(filename="check.tif", res=600, height=4200, width=4200, units="px")

pdf(file="Figures/Figure_S1.pdf")
#windows()
plotRGB(lgt, r=1, g=2, b=3, maxpixels=(4200^2), stretch="lin", ext=extent3)

points(sl_xydata, pch=24, cex=1, col="black", bg="red")
points(xydata, pch=21, cex=1, col="black", bg="yellow")

rect(724589.4, -2882507, 740445, -2874816, col="white")

sb2<-list(x=726860.1, y=-2879246)

scalebar(10000, xy = c(sb2$x, sb2$y) , type = "bar", divs=4, label=c(0, 5, 10), below="kilometers", font=2, col="black", adj=c(0.5,-1.2))

dev.off()

#inset map
library(oz)
pdf(file="wa_inset.pdf")
par(mar=c(0,0,0,0))
wa()
dev.off()

#select extent

ex<-locator(2)

extent3<-extent(c(ex$x[1], ex$x[2], ex$y[2], ex$y[1]))

# extent3
# class       : Extent 
# xmin        : 724589.4 
# xmax        : 801816.3 
# ymin        : -2882507 
# ymax        : -2809656 
