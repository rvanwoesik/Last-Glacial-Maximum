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
vers<-1    # this is the species choosing number, 10 is p.lobata, 1 is a hyacinthus

for (vers in 1:12){ ############### WARNING TAKES 3 hours
  
  set<-list1[vers]
  bnm<-list2[vers]
  setwd("E:/Work for RA/IUCN shapes")
  corals2<-read.dbf(paste("CORALS",set,".dbf",sep=""))
  corals<-corals2[order(corals2$binomial),]
  Species<-paste(corals[bnm,2])
  setwd("E:/Paleoreconstruction/plots May 16 2019")
  
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



R2<-rotater(R,center=210)
extent(R2)<-c(-180,180,-37,37)
Z2<-rotater(Z,center=210)
extent(Z2)<-c(-180,180,-37,37)

LGMshelf<-raster('E:/paleoreconstruction/sea-level/LGMshelf.nc')
extent(LGMshelf)<-c(c(-180,180,-70,70))
tLGMshelf<-rotater(LGMshelf>0,center=210)
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
redcol<-colorRampPalette(c("#00000000","red3"),alpha=T)(25)
bluecol<-colorRampPalette(c("#00000000","blue"),alpha=T)(25)







setwd("E:/Paleoreconstruction/plots May 16 2019")
tiff(paste(Species,"conf.png"),width=3900,height=1500, res = 300)
plot.map("world", transf=F, center=210 , col="burlywood",bg="white",ylim=c(-45,45),fill=TRUE,mar=c(2,5,2,0),add=F) #center is still 0
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col ="lightsteelblue1")

image(tLGMshelf,add=T,col=rev(c('lightsteelblue1','khaki1')),legend=F)

plot.map("world", transf=F,center=210 , col="burlywood",bg="white",ylim=c(-45,45),fill=TRUE,mar=c(2,5,2,0),add=T,xlab='longitude',ylab='latitude',lwd=1.5) #center is still 0

box(lwd=1.5)
compassRose(-155,-30,cex=.75,cex.dir=1.2,llwd=1.5)
axis(1,at=c(-150,-60,30,120,210),labels=c("",'','','',''), tck = .025,mgp=c(3,.3,0),cex.axis=.975,lwd=1.5)

mtext("90", 1, line=-1.5, adj=0.335,cex=.975)
mtext('270',1,line=-1.5,adj=0.836,cex=.975)

axis(2,at=c(-30,-15,0,15,30),labels=c('','','','',''), tck = .025,mgp=c(3,.405,0),cex.axis=.975,lwd=1.5)

mtext("-15", 2, line=-1.5, adj=0.315,cex=.975)
mtext('15',2,line=-1.5,adj=0.675,cex=.975)

#image(get(paste('RTDone',scenario,sep='')),maxpixels=3836160,add=T,col = c(color.palette2,color.palette1),breaks=palette.breaks)
image(R2,add=T,maxpixels=3836160,col=redcol,legend=F)  ### new habitat 
image(Z2,maxpixels=3836160,add=T,col=bluecol,legend=F)  ### currently inhabited
#image((losthab2>0),add=T,col=c("red"),maxpixels=3836160)

#color.legend(140,28,175,32,legend=round(seq(-845,336,length.out=3)), rect.col=c(color.palette4,color.palette3),cex=.975)
legend("topright",c('LGM shelf',"LGM coral","modern coral"),col=c('khaki1',"red3","blue"),pch=c(16,16,16),box.col=NA,bg="#00000000",cex=1.1)

#mtext("mm above sea-level", 1, line=-12.15, adj=0.99,cex=.975)


text(-72.31,-25.336,'Indian Ocean',cex=1.2)
text(87.3,3.4,'Pacific Ocean',cex=1.2)
scalebar(d=36,xy=c(66,-39),label=c(0,'',4000),cex=.9,type='bar',divs=4,below="kilometers",adj=c(0.5,-1.1))
dev.off()


}



################################### PLotting only the LGM for highlighting centers... maybe do the aggregate...
###################################
###############################
#####################################
#############################3
#











setwd("E:/Paleoreconstruction/plots jan3 2019")
tiff(paste(Species,"conf.png"),width=3900,height=1500, res = 300)
plot.map("world", transf=F, center=210 , col="burlywood",bg="white",ylim=c(-45,45),fill=TRUE,mar=c(2,5,2,0),add=F) #center is still 0
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col ="lightsteelblue1")

image(tLGMshelf,add=T,col=rev(c('lightsteelblue1','khaki1')),legend=F)

plot.map("world", transf=F,center=210 , col="burlywood",bg="white",ylim=c(-45,45),fill=TRUE,mar=c(2,5,2,0),add=T,xlab='longitude',ylab='latitude',lwd=1.5) #center is still 0

box(lwd=1.5)
compassRose(-155,-30,cex=.75,cex.dir=1.2,llwd=1.5)
axis(1,at=c(-150,-60,30,120,210),labels=c("",'','','',''), tck = .025,mgp=c(3,.3,0),cex.axis=.975,lwd=1.5)

mtext("90", 1, line=-1.5, adj=0.335,cex=.975)
mtext('270',1,line=-1.5,adj=0.836,cex=.975)

axis(2,at=c(-30,-15,0,15,30),labels=c('','','','',''), tck = .025,mgp=c(3,.405,0),cex.axis=.975,lwd=1.5)

mtext("-15", 2, line=-1.5, adj=0.315,cex=.975)
mtext('15',2,line=-1.5,adj=0.675,cex=.975)

#image(get(paste('RTDone',scenario,sep='')),maxpixels=3836160,add=T,col = c(color.palette2,color.palette1),breaks=palette.breaks)
image(R2,add=T,maxpixels=3836160,col=redcol,legend=F)  ### new habitat 
image(Z2,maxpixels=3836160,add=T,col=bluecol,legend=F)  ### currently inhabited
#image((losthab2>0),add=T,col=c("red"),maxpixels=3836160)

#color.legend(140,28,175,32,legend=round(seq(-845,336,length.out=3)), rect.col=c(color.palette4,color.palette3),cex=.975)
legend("topright",c('LGM shelf',"LGM coral","modern coral"),col=c('khaki1',"red3","blue"),pch=c(16,16,16),box.col=NA,bg="#00000000",cex=1.1)

#mtext("mm above sea-level", 1, line=-12.15, adj=0.99,cex=.975)


text(-72.31,-25.336,'Indian Ocean',cex=1.2)
text(87.3,3.4,'Pacific Ocean',cex=1.2)
scalebar(d=36,xy=c(66,-39),label=c(0,'',4000),cex=.9,type='bar',divs=4,below="kilometers",adj=c(0.5,-1.1))
dev.off()


}