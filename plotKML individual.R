#each spp KML plot
library(plotKML)
library(raster)

list1<-c(1,2,2,2,2,2)
list2<-c(19,125,143,176,191,331)
#read IUCN
list1<-c(1,1,2,2,2,2,2,1,2,2,2,2,list1)                       #species to choose from    ### important to change to the 6 most valuable species
list2<-c(58,30,65,129,140,180,238,233,395,398,455,470,list2)

vers<-4    # this is the species choosing number, 10 is p.lobata, 1 is a hyacinthus

riversBufferLGM<-readShapeSpatial('E:/paleoreconstruction/rivers buffer/riversBufferLGM.shp')
proj4string(riversBufferLGM)<-'+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0'
palette<-rev(colorRampPalette(c("red3","orange",'yellow','white'))(18))
setwd("E:/Paleoreconstruction/Stats plots feb 2020/Earth/LGM")

for (vers in 1:18){

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

  R<-mask(R,riversBufferLGM,inverse=T)
  
  
  R[R==0]<-NA
  Z[Z==0]<-NA
  
  LGM<-R
  Modern<-Z
  
setwd("E:/Paleoreconstruction/Stats plots feb 2020/Earth")
dir.create(paste(getwd(),'/',SPP,sep=''))
 
setwd(paste("E:/Paleoreconstruction/Stats plots feb 2020/Earth/",SPP,sep=''))
dir.create(paste(getwd(),'/','Modern',sep='')) 
setwd(paste("E:/Paleoreconstruction/Stats plots feb 2020/Earth/",SPP,'/','Modern',sep=''))
plotKML(Modern,colour_scale=palette,z.lim=c(0,25)) 
 
setwd(paste("E:/Paleoreconstruction/Stats plots feb 2020/Earth/",SPP,sep=''))
dir.create(paste(getwd(),'/','LGM',sep='')) 
setwd(paste("E:/Paleoreconstruction/Stats plots feb 2020/Earth/",SPP,'/','LGM',sep=''))
plotKML(LGM,colour_scale=palette,z.lim=c(0,25)) 

}





palette<-rev(colorRampPalette(c("red3","orange",'yellow','white'))(18))
setwd("E:/Paleoreconstruction/Stats plots feb 2020/Earth/LGM")
LGMKML<-R2

LGMKML[LGMKML==0]<-NA
plotKML(LGMKML,colour_scale=palette,z.lim=c(0,18))
