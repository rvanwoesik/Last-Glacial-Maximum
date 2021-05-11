## distance matrix between LGM and modern day reefs.
#>12 all spp? or sby each individual species. both? 
#>12 would give a reef comunity adaptation
#species specific is more accurate as far as gene flow is only by species. takes longer. 



library(raster)
library(foreign)
############# plots 
library(plotrix)
library(raster)
library(MASS)
library(audio)
library(sp)
library(foreign)
library(rgdal)
library(maptools)
library(rgeos)
library(doParallel)
library(rasterVis)
library(dismo)
library(plotKML)
library(SDMTools)
library(PBSmapping)
library(lme4)
library(blme)
library(mailR)
library(raster)
library(fields)

load("E:/Work for RA/Comps/load data comps 10-1-2014.RData")



list1<-c(1,1,2,2,2,2,2,1,2,2,2,2,1,2,2,2,2,2)                       #species to choose from
list2<-c(58,30,65,129,140,180,238,233,395,398,455,470,19,125,143,176,191,331)
# list1<-c(1,2,2,2,2,2)   
# list2<-c(19,125,143,176,191,331)

#vers<-4    # this is the species choosing number, 10 is p.lobata, 1 is a hyacinthus

for (vers in 1:18){ ############### WARNING TAKES 3 hours
  
  set<-list1[vers]
  bnm<-list2[vers]
  setwd("E:/Work for RA/IUCN shapes")
  corals2<-read.dbf(paste("CORALS",set,".dbf",sep=""))
  corals<-corals2[order(corals2$binomial),]
  Species<-paste(corals[bnm,2])
 # setwd("E:/Paleoreconstruction/plots jan3 2019")
  
  SPP<-Species
  #time<-"A2"  
  
  
  time<-"LGM"  #initial
  
  R<-raster(paste("E:/paleoreconstruction/paleoreconstruction/",SPP,"/",SPP," 1 ",time,' distribution.nc',sep=''))
  
  for (i in 2:25){ #2:25
    x<-raster(paste("E:/paleoreconstruction/paleoreconstruction/",SPP,"/",SPP," ",i," ",time,' distribution.nc',sep=''))
    R<-sum(R,x,na.rm=T)
    print(i)
  }
  
  time<-"initial"
  
  Z<-raster(paste("E:/paleoreconstruction/paleoreconstruction/",SPP,"/",SPP," 1 ",time,' dist.nc',sep=''))
  
  for (i in 2:25){ #2:25
    x<-raster(paste("E:/paleoreconstruction/paleoreconstruction/",SPP,"/",SPP," ",i," ",time,' dist.nc',sep=''))
    Z<-sum(Z,x,na.rm=T)
    print(i)
  }
  
  R12<-R>12
  Z12<-Z>12
  




  if(vers==1){ R2<-R12; Z2<-Z12}
  if(vers!=1){R2<-sum(R2,R12,na.rm=T);Z2<-sum(Z2,Z12,na.rm=T)}
}  






setwd("E:/Paleoreconstruction/Stats plots feb 2020/sums all greater than 12")
#writeRaster(R2,'LGM12.nc')
#writeRaster(Z2,'Modern12.nc')
R12<-raster("LGM12.nc" )
R12<-raster::extend(R12,c(-180,180,-45,45),value=NA)


a2R12<-aggregate(R12,3,fun=min,na.rm=T)  #to find dense LGM reefs. 1/4 degree cell size
cellStats(a2R12>0,sum)
a2R12[a2R12==0]<-NA


writeRaster(a2R12,'aggreagted LGM.nc')
# Z12<-raster('Modern12.nc')
# R12[R12==0]<-NA
# Z12[Z12==0]<-NA
# cellStats(R12,sum)
# 
# aR12<-aggregate(R12,4,fun=mean,na.rm=T) #2min instead of forever
# cellStats(aR12>0,sum)

a2R12<-raster('aggreagted LGM.nc')
a2R12b<-mask(a2R12,riversBufferLGM,inverse=T)
a2R12<-a2R12b
writeRaster(a2R12,'aggreagted LGM.nc')

# ncols<-200
# nrows<-1080
# ptm <- proc.time()
# r <- raster(ncol=ncols,nrow=nrows)
# r[] <- NA
# r[round(runif(1000,1,nrows*ncols))] <- 1
# dist <- raster::distance(r) 
# proc.time() - ptm

ptm <- proc.time()
dist2 <- raster::distance(a2R12)
proc.time() - ptm
setwd("E:/Paleoreconstruction/distance mat")
writeRaster(dist2,'18spp LGM reef distance agmin2.nc',overwrite=T)
setwd("E:/Paleoreconstruction/Stats plots feb 2020/sums all greater than 12")
writeRaster(dist2,'18dist from LGM reef.nc',overwrite=T)

setwd("E:/Paleoreconstruction/Stats plots feb 2020/sums all greater than 12")
Z12<-raster("Modern12.nc" )
# z <- raster(ncol=100,nrow=200)
# z[] <- NA
# z[round(runif(100,1,20000))] <- 1
zpt<-rasterToPoints(Z12,spatial=T)

#plot(dist)
dpts<-extract(dist2,zpt)
#points(zpt,col=dpts)
#zpt2<-zpt

zo2<-Z12
zo2[zpt]<-dpts

setwd("E:/Paleoreconstruction/Stats plots feb 2020/sums all greater than 12")
writeRaster(zo2,'dist from LGM reef.nc',overwrite=T)







## Hierarchical Clustering
setwd("E:/Paleoreconstruction/Stats plots feb 2020/sums all greater than 12")
a2R12<-raster('aggreagted LGM.nc')
R18<-rasterToPoints(a2R12)
pts<-R18[,c(1,2)]
library(maptools)
library(sp)
library(raster)
#data(wrld_simpl)
#pts<-spsample(a2R12,n=100,'random') #randomly select points in the world (this will be the major LGM reef points either all or by species)


setwd("E:/Paleoreconstruction/Stats plots feb 2020/sums all greater than 12")
dist2<-raster('18dist from LGM reef.nc')
plot(dist2)
points(pts,add=T,pch=16)

mydata<-pointDistance(pts,lonlat=T)
d <- dist(mydata, method = "euclidean") # distance matrix

#determine number of clusters
wssm<-matrix(nrow=15,ncol=100)
for(j in 1:100){
wss <- c()
for (i in 2:15) wss[i] <- sum(kmeans(d,centers=i)$withinss)
#plot(1:15, wss, type="b", xlab="Number of Clusters",ylab="Within groups sum of squares")
#looks like 6 is lowest number where most variance can be accounted for. 
wssm[,j]<-wss
}
plot(1:15, wssm[,1], type="b", xlab="Number of Clusters",ylab="Within groups sum of squares")
for(i in 2:100){
  lines(1:15, wssm[,i], type='b',xlab="Number of Clusters",ylab="Within groups sum of squares")
}
plot(1:15, rowMeans(wssm), type='b',xlab="Number of Clusters",ylab="Within groups sum of squares")


fit <- hclust(d, method="centroid")
plot(fit) # display dendogram
groups <- cutree(fit, k=6) # cut tree into 5 clusters
# draw dendogram with red borders around the 5 clusters
rect.hclust(fit, k=6, border="red")



plot(dist2) #plot geographically
points(pts,col=groups,pch=16)



setwd("E:/Paleoreconstruction/Stats plots feb 2020/sums all greater than 12")
#writeRaster(dist2,'18dist from LGM reef.nc',overwrite=T)






##################################### plot

setwd("E:/Paleoreconstruction/Stats plots feb 2020/sums all greater than 12")
#aggregated distance raster only of high density reefs (27*27km or 1/4 degree)
#plot(zo2)
setwd("E:/Paleoreconstruction/Stats plots feb 2020/sums all greater than 12")
dist2<-raster('18dist from LGM reef.nc')

library(raster)
library(foreign)
############# plots 
library(plotrix)
library(raster)
library(MASS)
library(audio)
library(sp)
library(foreign)
library(rgdal)
library(maptools)
library(rgeos)
library(doParallel)
library(rasterVis)
library(dismo)
library(plotKML)
library(SDMTools)
library(PBSmapping)
library(lme4)
library(blme)
library(mailR)
library(raster)
library(fields)

compassRose<-function(x,y,rot=0,cex=1,cex.dir=1,llwd=1) { 
  oldcex<-par(cex=cex) 
  mheight<-strheight("M") 
  xylim<-par("usr") 
  plotdim<-par("pin") 
  xmult<-(xylim[2]-xylim[1])/(xylim[4]-xylim[3])*plotdim[2]/plotdim[1] 
  point.angles<-seq(0,7*pi/4,by=pi/4)+pi*rot/180 
  crspans<-rep(c(mheight*3,mheight/2),4) 
  xpoints<-cos(point.angles)*crspans*xmult+x 
  ypoints<-sin(point.angles)*crspans+y 
  polygon(xpoints,ypoints,lwd=llwd) 
  txtxpoints<-cos(point.angles[c(1,3,5,7)])*1.33*crspans[1]*xmult+x 
  txtypoints<-sin(point.angles[c(1,3,5,7)])*1.33*crspans[1]+y 
  text(txtxpoints,txtypoints,c("E","N","W","S"),cex=cex.dir) 
  par(oldcex) 
} 

moll<-"+proj=moll"


plot.map<- function(database,center,transf=T,projectione=newproj,...){
  Obj <- map(database,...,plot=F)
  coord <- cbind(Obj[[1]],Obj[[2]])
  newproj <- "+proj=merc +lon_0=150 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs" #utm
  nextproj<-"+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0" #latlong
  moll<-"+proj=moll"
  # split up the coordinates
  id <- rle(!is.na(coord[,1]))
  id <- matrix(c(1,cumsum(id$lengths)),ncol=2,byrow=T)
  polygons <- apply(id,1,function(i){coord[i[1]:i[2],]})
  
  # split up polygons that differ too much
  polygons <- lapply(polygons,function(x){
    x[,1] <- x[,1] + center
    x[,1] <- ifelse(x[,1]>180,x[,1]-360,x[,1])
    if(sum(diff(x[,1])>300,na.rm=T) >0){
      id <- x[,1] < 0
      x <- rbind(x[id,],c(NA,NA),x[!id,])
    }
    x
  })
  # reconstruct the object
  polygons <- do.call(rbind,polygons)
  
  
  colnames(polygons)<-c("x",'y')
  polygons<-as.data.frame(polygons)
  z<-complete.cases(polygons)
  p<-z
  z<-cbind(z,z)
  polygons<-polygons[complete.cases(polygons),]
  coordinates(polygons)<-~x+y
  proj4string(polygons)<-CRS(nextproj)
  
  if(transf==T){ polygons<-spTransform(polygons,CRS(projectione))}
  
  z[p==F,]<-c(NA,NA)
  z[which(p==T),]<-coordinates(polygons)
  Obj[[1]] <- z[,1]
  Obj[[2]] <- z[,2]
  
  
  map(Obj,...)
}


source("E:/Work for RA/R functions/revrotate..R")
source("E:/Work for RA/R functions/rotater.R")



zo2<-rotater(zo2,center=210)
extent(zo2)<-c(-180,180,-37,37)

#LGMshelf<-raster('E:/paleoreconstruction/sea-level/LGMshelf.nc')
# extent(LGMshelf)<-c(c(-180,180,-70,70))
# tLGMshelf<-rotater(LGMshelf>0,center=210)
# extent(tLGMshelf)<-c(-180,180,-70,70)
# tLGMshelf<-crop(tLGMshelf,c(-180,180,-45,45))
# tLGMshelf[tLGMshelf==1]<-NA
# ################# FIGURES



bluered<-colorRampPalette(c("blue",'green','yellow','orange',"red3"),alpha=F)(100)



setwd("E:/Paleoreconstruction/distance mat")
zo2<-raster('dist from LGM reef.nc')   #aggregated distance raster only of high density reefs (27*27km or 1/4 degree)
#plot(zo2)


setwd("E:/Paleoreconstruction/Stats plots feb 2020/sums all greater than 12")
dist2<-raster('18dist from LGM reef.nc')
dist2<-rotater(dist2,center=210)
extent(dist2)<-c(-180,180,-45,45)
pgATL<-readOGR('E:/Paleoreconstruction/distance mat','pgATL')
dist2<-mask(dist2,pgATL,inverse=F)


# pts1<-pts
# pts1[,1]<-pts1[,1]+210-360
# x<-which(pts1[,1]<(-180))
# pts1[x,1]<-pts1[x,1]+360
# groups[which(groups==5)]<-7

setwd("E:/Paleoreconstruction/Stats plots feb 2020/sums all greater than 12")
tiff(paste("LGM reef dist mat18 nopts.png"),width=3900,height=1500, res = 300)
plot.map("world", transf=F, center=210 , col="burlywood",bg="white",ylim=c(-36.5,36.5),fill=TRUE,mar=c(2,5,2,0),add=F) #center is still 0
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col ="lightsteelblue1")
#image(tLGMshelf,add=T,col=rev(c('lightsteelblue1','khaki1')),legend=F)
image(dist2,add=T,maxpixels=3836160,col=rev(terrain.colors(100)),legend=F)  ### dist mat
plot.map("world", transf=F,center=210 , col="burlywood",bg="white",ylim=c(-36.5,36.5),fill=TRUE,mar=c(2,5,2,0),add=T,xlab='longitude',ylab='latitude',lwd=1.5) #center is still 0
box(lwd=1.5)
axis(1,at=c(-150,-60,30,120,210),labels=c("",'','','',''), tck = .025,mgp=c(3,.3,0),cex.axis=.975,lwd=1.5)
mtext("90", 1, line=-1.5, adj=0.335,cex=.975)
mtext('270',1,line=-1.5,adj=0.836,cex=.975)
#image(zo2,add=T,maxpixels=3836160,col=rev(terrain.colors(100)),legend=F)  ### dist mat
axis(2,at=c(-30,-15,0,15,30),labels=c('','','','',''), tck = .025,mgp=c(3,.405,0),cex.axis=.975,lwd=1.5)
mtext("-15", 2, line=-1.5, adj=0.315,cex=.975)
mtext('15',2,line=-1.5,adj=0.675,cex=.975)
color.legend(139,28,173,32,legend=round(seq(0,5000,length.out=3)), rect.col=rev(terrain.colors(100)),cex=.975)
#legend("topright",c('LGM shelf',"LGM coral","modern coral"),col=c('khaki1',"red3","blue"),pch=c(16,16,16),box.col=NA,bg="#00000000",cex=1.1)
mtext("kilometers", 1, line=-11.1, adj=0.965,cex=.975)
compassRose(-155,-23,cex=.75,cex.dir=1.2,llwd=1.5)
text(-72.31,-22.336,'Indian Ocean',cex=1.2)
text(87.3,3.4,'Pacific Ocean',cex=1.2)
#points(pts1,col=groups,pch=16,cex=2)
scalebar(d=4000,xy=c(66,-32.5),label=c(0,'',4000),cex=.9,type='bar',divs=4,below="kilometers",adj=c(0.5,-1.1))
dev.off()




#pgATL<-drawPoly(sp=TRUE, col='red', lwd=2)
#pgATL<-SpatialPolygonsDataFrame(pgATL,data='NA')
# pgATL<-as(pgATL, "SpatialPolygonsDataFrame")
# setwd("E:/Paleoreconstruction/distance mat")
# writePolyShape(pgATL,'pgATL.shp')
# writeOGR(pgATL,'pgATL',layer=1,driver="ESRI Shapefile")
# 










########################################


pointDist(x, y = NULL, distFunct = NULL, longLat = NULL, ...)
mydata <- matrix(rnorm(500), nrow = 10)
d <- dist(mydata, method = "euclidean") # distance matrix
fit <- hclust(d, method="ward")
plot(fit) # display dendogram
groups <- cutree(fit, k=5) # cut tree into 5 clusters
# draw dendogram with red borders around the 5 clusters
rect.hclust(fit, k=3, border="red")



