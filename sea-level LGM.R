# marmap LGM

library(marmap)

#depth<-getNOAA.bathy(-180,180,-37,37,resolution=4) #need new directory I think

setwd("E:/Paleoreconstruction/sea-level")

depth<-raster("high 5min resolution global.nc")


# sea-level was 134 meters lower in the LGM (Lambeck et al. 2014)
plot(depth<(-(134)))
plot(depth<(-(164)))

e<-c(-180,180,-37,37)

LGM<-((depth<(-(130)))-(depth<(-(160))));LGM<-crop(LGM,e)
NOW<-((depth<(-(0)))-(depth<(-(30))));NOW<-crop(NOW,e)


a<-raster::area(LGM)

LGM<-mask(a,LGM,maskval=0)
NOW<-mask(a,NOW,maskval=0)

NOWa<-cellStats(NOW,sum)
LGMa<-cellStats(LGM,sum)

#Contemporary shorelines have 9x more area within 30m of sea-level.
sealevelLGM130<-((depth<(-(130)))-(depth<(-(160))))

sealevelLGM130<-crop(sealevelLGM130,e)
extent(sealevelLGM130)<-e

plot(rotater(sealevelLGM130,center=210))
setwd("E:/Paleoreconstruction/sea-level")
#writeRaster(sealevelLGM130,'sealevelLGM130.nc',overwrite=T)



LGMshelf<-((depth<(-(130))))
plot(LGMshelf,col=c('orangered','lightblue'),legend = FALSE)
plot(wrld_simpl,add=TRUE,col='white')
plot(riversBuffer,add=T,col='lawngreen')
plot(riversBufferLGM,col='blue',add=T)

writeRaster(LGMshelf,'LGMshelf.nc')


