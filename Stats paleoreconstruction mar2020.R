#Stats for Paleo runs
#allStats need to 

#summarize all spp.? maybe average each spp individually. 


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
kd<-raster('E:/Work for RA/Light Attenuation kd/mean kd 1998-2007 9km.nc')

list1<-c(1,1,2,2,2,2,2,1,2,2,2,2)                       #species to choose from
list2<-c(58,30,65,129,140,180,238,233,395,398,455,470)
vers<-1    # this is the species choosing number, 10 is p.lobata, 1 is a hyacinthus

list1<-c(1,2,2,2,2,2)   
list2<-c(19,125,143,176,191,331)
vers<-1


#masks
#shelf size
mask_deep_bathyLGM<-raster('E:/paleoreconstruction/sea-level/sealevelLGM130.nc')
mask_deep_bathy<-raster('E:/paleoreconstruction/sea-level/modern 0-30m.nc')
extent(mask_deep_bathy)<-c(-180,180,-37,37)
mask_deep_bathyLGM<-resample(mask_deep_bathyLGM,pPARmeanFINAL)

a<-raster::area(mask_deep_bathyLGM)
LGMmsk1<-mask(a,mask_deep_bathyLGM,maskval=0)
LGMmsk<-(mask(LGMmsk1,ATL,maskval=1,inverse=T))

MODmsk1<-mask(a,mask_deep_bathy,maskval=NA,inverse=F)
MODmsk<-(mask(MODmsk1,ATL,maskval=1,inverse=T))

shallowAreaMltpl<-cellStats(MODmsk,sum)/cellStats(LGMmsk,sum)
shallowAreaMltpl
#~ 3.17 times more shallow water habitat available now compared to LGM over the whole planet





for(vers in c(1:6)){ 
  set<-list1[vers]
  bnm<-list2[vers]
  setwd("E:/Work for RA/IUCN shapes")
  corals2<-read.dbf(paste("CORALS",set,".dbf",sep=""))
  corals<-corals2[order(corals2$binomial),]
  Species<-paste(corals[bnm,2])
  SPP<-Species
  
  time<-"LGM"  
  LGM<-raster(paste("E:/paleoreconstruction/paleoreconstruction/",SPP,"/",SPP," 1 ",time,' distribution.nc',sep=''))
  for (i in 2:25){ #2:25
    x<-raster(paste("E:/paleoreconstruction/paleoreconstruction/",SPP,"/",SPP," ",i," ",time,' distribution.nc',sep=''))
    LGM<-stack(LGM,x,na.rm=T)
    if(vers==1){ print(i)}
  }
  #LGM2<-mask(LGM,riversBufferLGM,inverse=T)  
  setwd(paste("E:/paleoreconstruction/paleoreconstruction/",SPP,sep=''))
  writeRaster(LGM,paste(Species,'LGM stack.nc'),overwrite=T)
  print(SPP)
}










library(doParallel)
library(doSNOW)
cores<-4
cl <- makeCluster(cores)

#progress bar
registerDoSNOW(cl)
iterations <- 6
pb <- txtProgressBar(max = iterations, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)


#foreach (vers = 1:12,.packages=c('raster','rgdal','sp','foreign'), .options.snow = opts)%dopar%{ ########
 
for(vers in c(1:6)){
  set<-list1[vers]
  bnm<-list2[vers]
  setwd("E:/Work for RA/IUCN shapes")
  corals2<-read.dbf(paste("CORALS",set,".dbf",sep=""))
  corals<-corals2[order(corals2$binomial),]
  Species<-paste(corals[bnm,2])
  SPP<-Species

  time<-"LGM"  
LGM<-stack(paste("E:/paleoreconstruction/paleoreconstruction/",SPP,"/",SPP," LGM stack.nc",sep=''))
  
  time<-"initial"
  NOW<-raster(paste("E:/paleoreconstruction/paleoreconstruction/",SPP,"/",SPP," 1 ",time,' dist.nc',sep=''))
  for (i in 2:25){ #2:25
    x<-raster(paste("E:/paleoreconstruction/paleoreconstruction/",SPP,"/",SPP," ",i," ",time,' dist.nc',sep=''))
    NOW<-stack(NOW,x,na.rm=T)
    print(i)
  }
  
 func<-function(x,y){
   sum(x,-1*y,na.rm=T)
   }
 
 
 
 
 
 ################################# CHANGE MASKVAL TO NOT INCLUDE NA!!!!!!!!!!!!!!!!
 a<-raster::area(NOW)
 NOWa<-mask(a,NOW,maskval=NA)
 NOWa<-mask(NOWa,NOW,maskval=1,inverse=T)
 
 LGMa<-mask(a,LGM,maskval=c(NA))
 LGMa<-mask(LGMa,LGM,maskval=1,inverse=T)
 
 NOWas<-cellStats(NOWa,sum)
 LGMas<-cellStats(LGMa,sum)
 
 slick<-cbind(NOWas,LGMas)  #area of reefs km2
 
 LGMovrNOW<-LGMas/NOWas
 PpnOfNOWLost<-1-LGMas/NOWas
 slick<-cbind(slick,cbind(LGMovrNOW,PpnOfNOWLost))
 
 
difs<-(NOW+(-1*LGM))  #NOW-LGM overalpping area (negative means only occur in LGM, positive only MODERN, and 0 is no change)
  a<-raster::area(difs)
  difsa<-mask(a,difs)
  #difsas<-cellStats(difsa[difs<1],sum) 

    
overlapNOWonly<-c()  
  for(i in 1:25){
  x<-sum(difsa[[i]][difs[[i]]<1 & difs[[i]]>.1])#how much OVERLAPPING area is inhabited NOW, but not during the LGM. 
  overlapNOWonly[i]<-x
  print(x)
  }

overlapLGMonly<-c()  
  for(i in 1:25){
    x<-sum(difsa[[i]][difs[[i]]<0]) #how much OVErLAPPING area was inhabited during LGM, but is no longer inhabited. 
    overlapLGMonly[i]<-x 
    print(x)
  }
  
slick<-cbind(slick,cbind(overlapLGMonly,overlapNOWonly))
 
 
  
#how many modern reefs of this spp. were restricted by SST min during the LGM...  
  #SST min
  SST_min_LGM<-raster('E:/paleoreconstruction/CMIP5 LGM tos/SSTminLGM_CMIP5.nc') 
  #SST_min_current
  
  #NOWall<-sum(NOW,na.rm=T)
  #NOWall[NOWall==0]<-NA
  minTolerance<-18
  SST_min_LGMr<-resample(SST_min_LGM,NOW)
  NOW_sst18<-mask(NOW,SST_min_LGMr>minTolerance,maskvalue=0)
  
  #NOWa
  #NOWas
  NOW_sst18a<-mask(a,NOW_sst18,maskval=NA)
  NOW_sst18a<-mask(NOW_sst18a,NOW_sst18,maskval=0)
  NOW_sst18as<-cellStats(NOW_sst18a,sum)
  
  SSTminLGM_lost_km2<-NOWas-NOW_sst18as #area in km2 of modern reefs of this spp.  restricted by SST min during the LGM... 
  SSTminLGM_lost_ppn<-SSTminLGM_lost_km2/NOWas #proportion of reefs lost from SST min
  
  slick<-cbind(slick,cbind(SSTminLGM_lost_km2,SSTminLGM_lost_ppn))
  
  
  
  
  
  
  
  #rivers, should be the same size, may overalp more but whatever. 
#  riversBufferLGM<-readShapeSpatial('E:/paleoreconstruction/rivers buffer/riversBufferLGM.shp')
#  riversBuffer
  
  
  
  
  
  #PAR
  PAR.limit=21.6
  
  #PAR.mean<-raster("C:/Users/Chris/Desktop/van Woesik/PAR global data/SWFMO_PAR.CR.timeAverage1998-2007.nc") #used to determine effects of light
  #pPARmeanFINAL<-crop(PAR.mean,extent(PLD))*PLD 
#prediction raster

  pastPARmean<-raster('E:/paleoreconstruction/insolation/LGM PARmean.nc')  #PAR mean for LGM 
  change.kd_to_PLD<-function(x,Zmax=3){ #for rasters
    x2<-.6677*x^.6763 #kdPAR~kd490 pierson et al. 2007
    Y2<-exp((-1*x2)*Zmax) #percent light at depth zmax
    return(Y2)
  }
  PLD<-change.kd_to_PLD(kd) #assume same turbidity, resampled to coarser res maintains regional patterns
  pastPLD_mean<-pastPARmean*PLD
  #pastPLD_mean<PAR.limit
#resample  
 
  
  NOW_LGMPARlim<-mask(NOW,pastPLD_mean<PAR.limit,maskvalue=0)
  #NOWa
  #NOWas
  NOW_LGMPARlima<-mask(a,NOW_LGMPARlim)
  NOW_LGMPARlimas<-cellStats(NOW_LGMPARlima,sum)
  
  
  
  PARminLGM_lost_km2<-NOWas-(NOWas-NOW_LGMPARlimas) #area in km2 of modern reefs of this spp.  restricted by SST min during the LGM... 
  PARminLGM_lost_ppn<-1-PARminLGM_lost_km2/NOWas #proportion of reefs lost from SST min
  
  slick<-cbind(slick,cbind(PARminLGM_lost_km2,PARminLGM_lost_ppn))
  
  
  
  
  #change in habitat adjusting for changes in bathymetry
  LGMas/(NOWas/shallowAreaMltpl)
  adjustedBathy<-LGMas/NOWas*shallowAreaMltpl
  
  slick<-cbind(slick,adjustedBathy)
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  conf95<-function(x){
    1.96*sd(x)/sqrt(length(x))
  }
  means<-apply(slick,2,mean)
  confs<-apply(slick,2,conf95)
  
  saveme<-rbind(slick,rbind(means,confs))
  
  setwd("E:/Paleoreconstruction/stats Oct 11 2019")
  write.csv(saveme,paste('paleostats',SPP,'.csv'))
  print(paste(SPP,'DONE'))
  
  
  
  
  
  
  
  
}

close(pb)
stopCluster(cl)
