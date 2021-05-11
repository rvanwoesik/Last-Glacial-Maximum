#Conf plotting LGM

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



list1<-c(1,1,2,2,2,2,2,1,2,2,2,2)                       #species to choose from
list2<-c(58,30,65,129,140,180,238,233,395,398,455,470)
vers<-4    # this is the species choosing number, 10 is p.lobata, 1 is a hyacinthus

for (vers in 1:12){ ############### WARNING TAKES 3 hours
  
  set<-list1[vers]
  bnm<-list2[vers]
  setwd("E:/Work for RA/IUCN shapes")
  corals2<-read.dbf(paste("CORALS",set,".dbf",sep=""))
  corals<-corals2[order(corals2$binomial),]
  Species<-paste(corals[bnm,2])
  setwd("E:/Paleoreconstruction/plots Oct 11 2019")
  
  SPP<-Species
  #time<-"A2"  
  
  
  time<-"LGM"  
  R<-stack(paste("E:/paleoreconstruction/paleoreconstruction/",SPP,"/",SPP," LGM stack.nc",sep=''))
  R<-sum(R)
  
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
  
  
  
  
  
  
  
  
  
  setwd("E:/Paleoreconstruction/plots Oct 11 2019")
   # writeRaster(Z2,'Modern12.nc')
   # writeRaster(R2,'LGM12.nc')
   # 

 Z2<-raster('Modern12.nc')
 R2<-raster('LGM12.nc')
# 
#   writeRaster(R2,'LGM12.tif')
#   writeOGR(R2, "LGM12.kml", layer='R2', driver="KML") 
  
  riversBufferLGM<-readShapeSpatial('E:/paleoreconstruction/rivers buffer/riversBufferLGM.shp')
  proj4string(riversBufferLGM)<-'+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0'
  proj4string(riversBuffer)<-'+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0'
  
  R2<-mask(R2,riversBufferLGM,inverse=T)
  
  palette<-rev(colorRampPalette(c("red3","orange",'yellow','white'))(12))
  setwd("E:/Paleoreconstruction/plots Oct 11 2019/Earth/LGM")
  LGMKML<-R2
  LGMKML[LGMKML==0]<-NA
  plotKML(LGMKML,colour_scale=palette,z.lim=c(0,12))
  setwd("E:/Paleoreconstruction/plots Oct 11 2019/Earth/Modern")
  MODKML<-Z2
  MODKML[MODKML==0]<-NA
  plotKML(MODKML,colour_scale=palette,z.lim=c(0,12))
  
  setwd("E:/Paleoreconstruction/plots Oct 11 2019/Earth/RiversLGM")
    plotKML(riversBufferLGM)
  setwd("E:/Paleoreconstruction/plots Oct 11 2019/Earth/Rivers Modern")
  ( p.df <- data.frame( ID=1:length(riversBuffer)) ) 
  modRivKml<-SpatialPolygonsDataFrame(riversBuffer,data=p.df)
  plotKML(modRivKml)
  
  
  R3<-rotater(R2,center=210+32)
  extent(R3)<-c(-180,180,-37,37)
  Z3<-rotater(Z2,center=210+32)
  extent(Z3)<-c(-180,180,-37,37)
  

  
  LGMshelf<-raster('E:/paleoreconstruction/sea-level/LGMshelf.nc')
  extent(LGMshelf)<-c(c(-180,180,-70,70))
  tLGMshelf<-rotater(LGMshelf>0,center=210+32)
  extent(tLGMshelf)<-c(-180,180,-70,70)
  tLGMshelf<-crop(tLGMshelf,c(-180,180,-45,45))
  tLGMshelf[tLGMshelf==1]<-NA
  # ################# FIGURES
  
  
  makeTransparent<-function(someColor, alpha=100)
  {
    newColor<-col2rgb(someColor)
    apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
                                                blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
  }
  makeTransparent('red3',1)
  redcol<-colorRampPalette(c("#00000000","red3"),alpha=T)(12)
  bluecol<-colorRampPalette(c("#00000000","blue"),alpha=T)(12)
  
  
  
  
  
  
  
  
  
  
  #lgmbound<-boundaries(tLGMshelf,type='outer',directions=4)
  
  setwd("E:/Paleoreconstruction/Stats plots feb 2020")
  tiff(paste("overall conf zoom zoneA mar2020.tif"),width=3000,height=2000, res = 300)
  #plot.map("world", transf=F, center=210+32 , col="burlywood",bg="white",ylim=c(-25,25),fill=TRUE,mar=c(2,5,2,0),add=T) #center is still 0
  #,xlim=c(210-62,2210+16.5)
 
  
  image(tLGMshelf,add=F,col=rev(c('lightsteelblue1','khaki1')),xlim=c(-35,40),ylim=c(-25,25), axes=F,xlab='',ylab='')
  #rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col ="lightsteelblue1")
  image(tLGMshelf,add=F,col=rev(c('lightsteelblue1','khaki1')),xlim=c(-35,40),ylim=c(-25,25), axes=F,xlab='',ylab='')
  plot.map("world", transf=F, center=210+32 , col="burlywood",bg="white",ylim=c(-25,25),fill=TRUE,mar=c(2,5,2,0),add=T) #center is still 0
  
  box(lwd=1.5)
  compassRose(-28,-17,cex=.75,cex.dir=1.2,llwd=1.5)
  axis(1,at=c(-18,2,22),labels=c("",'',''), tck = .025,mgp=c(3,.3,0),cex.axis=.975,lwd=1.5)
  
  mtext("100", 1, line=-1.5, adj=0.22,cex=.975)
  mtext("120", 1, line=-1.5, adj=0.49,cex=.975)
  mtext('140',1,line=-1.5,adj=0.77,cex=.975)
  
  axis(2,at=c(-15,0,15),labels=c('','',''), tck = .025,mgp=c(3,.405,0),cex.axis=.975,lwd=1.5)
  
  mtext("15", 2, line=-1.5, adj=0.817,cex=.975)
  mtext("0", 2, line=-1.5, adj=0.5,cex=.975)
  #mtext('140',2,line=-1.5,adj=0.675,cex=.975)
  
  #image(get(paste('RTDone',scenario,sep='')),maxpixels=3836160,add=T,col = c(color.palette2,color.palette1),breaks=palette.breaks)
  image(R3,add=T,maxpixels=3836160,col=redcol,legend=F)  ### new habitat 
  image(Z3,maxpixels=3836160,add=T,col=bluecol,legend=F)  ### currently inhabited
  #image((losthab2>0),add=T,col=c("red"),maxpixels=3836160)
  
  legend("topright",c('LGM land',"LGM coral","Modern coral"),col=c('khaki1',"red3","blue"),pch=c(16,16,16),box.col=NA,bg="#00000000",cex=1.1)
  legend("topright",c('LGM land',"LGM coral","Modern coral"),col=c('black',"black","black"),pch=c(1,NA,NA),box.col=NA,bg="#00000000",cex=1.1)
  color.legend(9,19,20,21,legend=c(0,6,12,18), rect.col=redcol,cex=.7)
  color.legend(9,15,20,17,legend=c(0,6,12,18), rect.col=bluecol,cex=.7)
  
  
  #mtext("mm above sea-level", 1, line=-12.15, adj=0.99,cex=.975)
  
  
  text(-14,-14,'Indian Ocean',cex=1.2)
  text(20,3,'Pacific Ocean',cex=1.2)
  scalebar(d=500,xy=c(11,-22),label=c(0,'',500),cex=.9,type='bar',divs=4,below="kilometers",adj=c(0.5,-1.1))
  dev.off()
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  

##################### zoom 2 ###################
  
  #setwd("E:/Paleoreconstruction/plots Oct 11 2019")
  tiff(paste("overall conf zoom seychelles mar2020.tif"),width=2000,height=2000, res = 300)

  image(LGMshelf,add=F,col=rev(c('white','khaki1')),legend=F,xlim=c(35,90),ylim=c(-30,25), axes=F,xlab='',ylab='')
  plot.map("world", transf=F, center=0 , col="burlywood",bg="white",ylim=c(-25,25),fill=TRUE,mar=c(2,5,2,0),add=T) #center is still 0
  
  box(lwd=1.5)
  compassRose(83,-3,cex=.75,cex.dir=1.2,llwd=1.5)
  axis(1,at=c(50,65,80),labels=c("",'',''), tck = .025,mgp=c(3,.3,0),cex.axis=.975,lwd=1.5)
  
  mtext("50", 1, line=-1.5, adj=0.26,cex=.975)
  mtext("65", 1, line=-1.5, adj=0.545,cex=.975)
  mtext('80',1,line=-1.5,adj=0.83,cex=.975)
  
  axis(2,at=c(-15,0,15),labels=c('','',''), tck = .025,mgp=c(3,.405,0),cex.axis=.975,lwd=1.5)
  
  mtext("15", 2, line=-1.5, adj=0.827,cex=.975)
  mtext("0", 2, line=-1.5, adj=0.55,cex=.975)
  mtext('-15',2,line=-1.5,adj=0.255,cex=.975)
  
  image(R2,add=T,maxpixels=3836160,col=redcol,legend=F)  ### new habitat 
  image(Z2,maxpixels=3836160,add=T,col=bluecol,legend=F)  ### currently inhabited
  #image((losthab2>0),add=T,col=c("red"),maxpixels=3836160)
  legend("topright",c('LGM land',"LGM coral","Modern coral"),col=c('khaki1',"red3","blue"),pch=c(16,16,16),box.col=NA,bg="#00000000",cex=1.1)
  legend("topright",c('LGM land',"LGM coral","Modern coral"),col=c('black',"black","black"),pch=c(1,NA,NA),box.col=NA,bg="#00000000",cex=1.1)
  color.legend(75,-19,85,-17,legend=c(0,6,12,18), rect.col=redcol,cex=.7)
  color.legend(75,-23,85,-21,legend=c(0,6,12,18), rect.col=bluecol,cex=.7)

  text(75,-12,'Indian Ocean',cex=1.2)
  #text(20,3,'Pacific Ocean',cex=1.2)
  scalebar(d=500,xy=c(83,-27),label=c(0,'',500),cex=.9,type='bar',divs=4,below="kilometers",adj=c(0.5,-1.1))
  dev.off()
  


  
  
  
  
  
  
  
  
  R4<-R2
  extent(R4)<-c(0,360,-37,37)
  Z4<-(Z2)
  extent(Z4)<-c(0,360,-37,37)
  
  LGMshelf<-raster('E:/paleoreconstruction/sea-level/LGMshelf.nc')
  extent(LGMshelf)<-c(c(-180,180,-70,70))
  tLGMshelf2<-LGMshelf
  extent(tLGMshelf2)<-c(0,360,-70,70)
  #tLGMshelf2<-crop(tLGMshelf,c(0,360,-45,45))
  #tLGMshelf2[tLGMshelf2==1]<-NA
  
  
  
  
  ############################### zoom 3 ##############
 # setwd("E:/Paleoreconstruction/plots Oct 11 2019")
  tiff(paste("overall conf zoom central eastern pacific mar2020.tif"),width=3200,height=2000, res = 300)
  
  image(tLGMshelf2,add=F,col=rev(c('white','khaki1')),legend=F,xlim=c(0,110),ylim=c(-30,25), axes=F,xlab='',ylab='')
  plot.map("world", transf=F, center=180 , col="burlywood",bg="white",ylim=c(-25,25),fill=TRUE,mar=c(2,5,2,0),add=T) #center is still 0
  
  box(lwd=1.5)
  compassRose(83,-8,cex=.75,cex.dir=1.2,llwd=1.5)
  axis(1,at=c(20,50,80),labels=c("",'',''), tck = .025,mgp=c(3,.3,0),cex.axis=.975,lwd=1.5)
  
  mtext("-160", 1, line=-1.5, adj=0.17,cex=.975)
  mtext("-130", 1, line=-1.5, adj=0.45,cex=.975)
  mtext('-100',1,line=-1.5,adj=0.73,cex=.975)
  
  axis(2,at=c(-15,0,15),labels=c('','',''), tck = .025,mgp=c(3,.405,0),cex.axis=.975,lwd=1.5)
  
  mtext("15", 2, line=-1.5, adj=0.827,cex=.975)
  mtext("0", 2, line=-1.5, adj=0.55,cex=.975)
  mtext('-15',2,line=-1.5,adj=0.255,cex=.975)
  
  image(R4,add=T,maxpixels=3836160,col=redcol,legend=F)  ### new habitat 
  image(Z4,maxpixels=3836160,add=T,col=bluecol,legend=F)  ### currently inhabited
  #image((losthab2>0),add=T,col=c("red"),maxpixels=3836160)
  legend("topleft",c('LGM land',"LGM coral","Modern coral"),col=c('khaki1',"red3","blue"),pch=c(16,16,16),box.col=NA,bg="#00000000",cex=1.1)
  legend("topleft",c('LGM land',"LGM coral","Modern coral"),col=c('black',"black","black"),pch=c(1,NA,NA),box.col=NA,bg="#00000000",cex=1.1)
  color.legend(75,-19,90,-17,legend=c(0,6,12,18), rect.col=redcol,cex=.7)
  color.legend(75,-23,90,-21,legend=c(0,6,12,18), rect.col=bluecol,cex=.7)
  
  #text(75,-12,'Indian Ocean',cex=1.2)
  text(180,3,'Pacific Ocean',cex=1.2)
  scalebar(d=500,xy=c(100,-27),label=c(0,'',500),cex=.9,type='bar',divs=4,below="kilometers",adj=c(0.5,-1.1))
  dev.off()