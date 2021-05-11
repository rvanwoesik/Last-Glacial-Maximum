#solar insolation for LGM

dat<-read.csv("E:/Paleoreconstruction/insolation/watts m2 by lat.csv", header=FALSE)
#change in w/m2 during LGM (Bush and Philander 1999)
plot(dat[,1],dat[,2])
colnames(dat)<-c('lat','change')

#can calculate the change in mean insolation. But what about summer, winter changes / range changes...

library(palinsol)

#t=time years after 1950
#solution One of ber78, ber90 or la04
orbit<- astro(t=-21000,solution=la04) #LGM using the Laskar 2004 solution
#eps is obliquity, ecc is eccentricty, and varpi is solar longitude of perhelion (precession)
  

#S0 1365 the total solar irradiance (could convert to PAR here)
# True solar longitude is measured in radians:
# pi/2 for June solstice
# pi for September equinox
# 3 * pi/2 for December solstice
# 0 for Spring equinox

colset<-c()
for(lata in seq(-89.5,89.5,1)){
for(longa in c(pi/2,pi,3*pi/2,0)){
x<-(Insol (orbit,long=longa, lat=lata*pi/180,S0=1365))
colset<-c(colset,x)
}}
colset<-matrix(colset,nrow=180,byrow=T)
plot(seq(-89.5, 89.5, 1),rowMeans(colset))



orbit2<- astro(t=0,solution=la04) #LGM using the Laskar 2004 solution
colset2<-c()
for(lata in seq(-89.5,89.5,1)){
  for(longa in c(pi/2,pi,3*pi/2,0)){
    x<-(Insol (orbit2,long=longa, lat=lata*pi/180,S0=1365))
    colset2<-c(colset2,x)
  }}
colset2<-matrix(colset2,nrow=180,byrow=T)
plot(seq(-89.5, 89.5, 1),rowMeans(colset2))


#difference between LGM and now. 
plot(seq(-89.5, 89.5, 1),rowMeans(colset)-rowMeans(colset2))

#as a percentage to multiply PAR by lat. 
plot(seq(-89.5, 89.5, 1),rowMeans(colset)/rowMeans(colset2))
#poles have lower irradiance and equator has higher irradiance. FOR MEAN PAR AND MIN avg PAR
pcntLatDif<-rowMeans(colset)/rowMeans(colset2)



par(mfrow=c(3,1))
par(oma = c(4, 4, 0, 0)) # make room (i.e. the 4's) for the overall x and y axis titles
par(mar = c(2, 2, 1, 1)) 
plot(seq(-89.5, 89.5, 1),(colset[,1])/(colset2[,1]),type='l',ylim=c(.98,1.08),lwd=2,xlab='',xaxt='n',ylab='',cex.axis=1.4)
abline(v=c(37,-37))
plot(seq(-89.5, 89.5, 1),(colset[,3])/(colset2[,3]),type='l',ylim=c(.98,1.08),lwd=2,xlab='',ylab='',cex.lab=1.5,xaxt='n',cex.axis=1.4)
abline(v=c(37,-37))
plot(seq(-89.5, 89.5, 1),rowMeans(colset)/rowMeans(colset2),type='l',lwd=2,xlab='',ylab='',cex.axis=1.4)
abline(v=c(37,-37))
mtext('Latitude', side = 1, outer = TRUE, line = 2)
mtext('Proportion of modern insolation', side = 2, outer = TRUE, line = 2)

#top is june, mid december, bottom yearly avg. 


lgmrange<-apply(colset,1,range)[2,]-apply(colset,1,range)[1,]
modrange<-apply(colset2,1,range)[2,]-apply(colset2,1,range)[1,]

modrange-lgmrange
###################################################################################
#fit to raster ranges
orbit<- astro(t=-21000,solution=la04) #LGM using the Laskar 2004 solution
sequ<-seq(-37,37,length.out=888)

colset<-c()
for(lata in sequ){
  for(longa in c(pi/2,pi,3*pi/2,0)){
    x<-(Insol (orbit,long=longa, lat=lata*pi/180,S0=1365))
    colset<-c(colset,x)
  }}
colset<-matrix(colset,nrow=888,byrow=T)
plot(sequ,rowMeans(colset))



orbit2<- astro(t=0,solution=la04) #LGM using the Laskar 2004 solution
colset2<-c()
for(lata in sequ){
  for(longa in c(pi/2,pi,3*pi/2,0)){
    x<-(Insol (orbit2,long=longa, lat=lata*pi/180,S0=1365))
    colset2<-c(colset2,x)
  }}
colset2<-matrix(colset2,nrow=888,byrow=T)
plot(sequ,rowMeans(colset2))


#difference between LGM and now. 
plot(sequ,rowMeans(colset)-rowMeans(colset2))

#as a percentage to multiply PAR by lat. 
plot(sequ,rowMeans(colset)/rowMeans(colset2))
#poles have lower irradiance and equator has higher irradiance. FOR MEAN PAR AND MIN avg PAR
pcntLatDif<-rowMeans(colset)/rowMeans(colset2)



#for range in PAR
load("E:/work for RA/Comps/load data comps 10-1-2014.RData") #load initial data
plot(Future$PAR_range)

range1<-(abs(apply(colset,1,range)[1,]-apply(colset,1,range)[2,]))
range2<-(abs(apply(colset2,1,range)[1,]-apply(colset2,1,range)[2,]))

plot(range1-range2) #less range in PAR at higher latitudes in past. 


plot(range1/range2)   #multiply PAR range by this lat difference to gain LGM PAR range. 
pcntRangeLatDif<-(range1/range2)

#Future$PAR_range[is.na(Future$PAR_range)]<-0
temp<-Future$PAR_range;values(temp)<-NA
values(temp)<-(sweep(as.matrix(Future$PAR_range),MARGIN=1,pcntRangeLatDif,FUN = "*"))
temp[temp==0]<-NA


plot(temp)
setwd("E:/Paleoreconstruction/insolation")
#writeRaster(temp,'LGM PAR_range.nc')



#################################### mean + turbidity
kd<-raster("E:/work for RA/Light Attenuation kd/mean kd 1998-2007 9km.nc ")#average light attenuation 490 coefficiten to 1998 to 1007 mean of each layer NA values removed. resolution = .0833333 (9km)
PAR.mean<-crop(PAR.mean,extent(kd))

temp2<-PAR.mean
values(temp2)<-NA
values(temp2)<-(sweep(as.matrix(PAR.mean),MARGIN=1,pcntLatDif,FUN = "*"))

plot((temp2))
setwd("E:/Paleoreconstruction/insolation")
#writeRaster(temp2,'LGM PARmean.nc')



