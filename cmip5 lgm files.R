# CMIP5 LGM gcm tos
#https://esgf-node.llnl.gov/search/cmip5/
#LGM | Oclim,Omon | lgm | tos | CMIP5 | #find tos

setwd("E:/Paleoreconstruction/CMIP5 LGM tos")
tos1<-stack("tos_Omon_GISS-E2-R_lgm_r1i1p151_300001-302412.nc")
tos2<-stack("tos_Omon_GISS-E2-R_lgm_r1i1p151_302501-304912.nc")
tos3<-stack("tos_Omon_GISS-E2-R_lgm_r1i1p151_305001-307412.nc")
tos4<-stack("tos_Omon_GISS-E2-R_lgm_r1i1p151_307501-309912.nc")

LGM<-stack(tos1,tos2,tos3,tos4)
LGM<-LGM-273.15
#extent(LGM)<-c(-180,180,-90,90)
e<-c(-180,180,-37,37)

#loop through each year and take range and min, then average among the 100 years.
for (i in 1:100){
  print(paste(i,'/100'))
  #setTxtProgressBar(pb, i) 
  year<-subset(LGM,seq(1,1200,12)[i]:seq(12,1200,12)[i])
  rangeyr<-range(year)
  
  if(i==1){
    minLGM<-min(year)
    rangeLGM<-rangeyr$range_max-rangeyr$range_min
  }
  
  if(i!=1){
    minLGM<-stack(minLGM,min(year))
    rangeLGM<-stack(rangeLGM,rangeyr$range_max-rangeyr$range_min)
  }
}
  
LGMminSST<-mean(minLGM)
LGMrangeSST<-mean(rangeLGM)

LGMminSST<-rotate(LGMminSST)
LGMrangeSST<-rotate(LGMrangeSST)

# setwd("E:/Paleoreconstruction/CMIP5 LGM tos")
# writeRaster(LGMminSST, 'SSTminLGM_CMIP5.nc')
# writeRaster(LGMrangeSST, 'SSTrangeLGM_CMIP5.nc')

LGMminSST<-crop(LGMminSST,e)

################comparison

setwd("E:/Paleoreconstruction")
minLGMMARGO<-raster('sstminLGMcc.nc')

setwd("E:/Paleoreconstruction/climap")
dataug<-raster("E:/Paleoreconstruction/climap/data aug.nc")
datfeb<-raster("E:/Paleoreconstruction/climap/data feb.nc")
plot(min(datfeb,dataug))
ClimapLGMmin<-min(datfeb,dataug)

plot(ClimapLGMmin<18) #CLIMAP only
plot((LGMminSST)<18) #CMIP5
plot(minLGMMARGO<18) #MARGO CLIMAP hybrid



# margjfm<-raster('Margo09.nc',varname='SSTJFM')
# margjas<-raster('Margo09.nc',varname='SSTJAS')
# plot(rotate(margjas-margjfm))
# MargoLGMrange<-rotate(margjas-margjfm)