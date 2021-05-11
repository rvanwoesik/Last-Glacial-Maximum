#stats overall

library(readr)
library(rgdal)
library(foreign)


list1<-c(1,2,2,2,2,2)   
list2<-c(19,125,143,176,191,331)
#read IUCN 
list1<-c(1,1,2,2,2,2,2,1,2,2,2,2,list1)                       #species to choose from    ### important to change to the 6 most valuable species
list2<-c(58,30,65,129,140,180,238,233,395,398,455,470,list2)


saveme1<-matrix(ncol=12,nrow=36)
splist<-c()
for(vers in 1:18){
#vers<-1  saveme1<-as.data.frame(saveme1)
set<-list1[vers]
bnm<-list2[vers]
setwd("E:/Work for RA/IUCN shapes")
corals2<-read.dbf(paste("CORALS",set,".dbf",sep=""))
corals<-corals2[order(corals2$binomial),]
Species<-paste(corals[bnm,2])
SPP<-Species

setwd("E:/Paleoreconstruction/stats Oct 11 2019")
x<-read.csv(paste('paleostats',Species,".csv" ))
x1<-x[26:27,]
saveme1[vers*2-1,]<-as.numeric(x1[1,])
saveme1[vers*2,]<-as.numeric(x1[2,])

colnames(saveme1)<-colnames(x1)
splist[vers]<-SPP
}
rownames(saveme1)<-rep(splist,each=2)

colMeans(saveme1)

saveme1<-as.data.frame(saveme1)

#attach(saveme)

plot(saveme1$NOWas[seq(1,36,2)],saveme1$LGMas[seq(1,36,2)]/saveme1$NOWas[seq(1,36,2)])
plot(saveme1$adjustedBathy[seq(1,36,2)],saveme1$NOWas[seq(1,36,2)])
plot(saveme1$LGMas[seq(1,36,2)],saveme1$adjustedBathy[seq(1,36,2)])



par(mar=c(5,5,2,1))
plot(saveme1$adjustedBathy[seq(1,24,2)],saveme1$NOWas[seq(1,24,2)],xlab='LGM relative abundance',ylab=expression(paste('Modern distribution ( ',km^2,')',sep='')))

NOWavg<-saveme1$NOWas[seq(1,24,2)]
ABavg<-saveme1$adjustedBathy[seq(1,24,2)]
reg1<-lm(log(NOWavg)~ABavg)

xnew<-data.frame(ABavg=seq(.2,2.7,0.01))
preda<-predict(reg1,newdata=xnew,interval='confidence')
lines(xnew$ABavg,exp(preda[,'fit']),lwd=2)
lines(xnew$ABavg,exp(preda[,'lwr']),lwd=1,lty=2)
lines(xnew$ABavg,exp(preda[,'upr']),lwd=1,lty=2)










########################################## ALL POINTS

#read IUCN 
# list1<-c(1,1,2,2,2,2,2,1,2,2,2,2)                       #species to choose from    ### important to change to the 6 most valuable species
# list2<-c(58,30,65,129,140,180,238,233,395,398,455,470)

saveme<-data.frame()
splist<-c()
for(vers in 1:18){
  #vers<-1  
  set<-list1[vers]
  bnm<-list2[vers]
  setwd("E:/Work for RA/IUCN shapes")
  corals2<-read.dbf(paste("CORALS",set,".dbf",sep=""))
  corals<-corals2[order(corals2$binomial),]
  Species<-paste(corals[bnm,2])
  SPP<-Species
  
  setwd("E:/Paleoreconstruction/stats Oct 11 2019")
  x<-read.csv(paste('paleostats',Species,".csv" ))
  x1<-x[c(-26,-27),]
  saveme<-rbind(saveme,x1)
  colnames(saveme)<-colnames(x1)
  splist[vers]<-SPP
}

#rownames(saveme)<-rep(splist,25)
#saveme<-as.data.frame(saveme)
saveme<-saveme[,-1]
colMeans(saveme)

#skips<-which(saveme$LGMas==0)
#saveme<-saveme[saveme$LGMas!=0,]




plot(saveme$adjustedBathy,saveme$NOWas)
#plot(saveme$LGMas,saveme$NOWas)
plot(saveme$LGMas,saveme$adjustedBathy)
plot(saveme$NOWas,saveme$LGMas)

cols<-(rainbow(18))
cols<-rep(cols,each=25)
#cols<-cols[-skips]



setwd("E:/Paleoreconstruction/Stats plots feb 2020")
tiff(paste("relative vs modern mar20.png"),width=2300,height=1900, res = 300)
par(mar=c(5,5,2,1))
plot(saveme$adjustedBathy,saveme$NOWas,xlab='LGM relative abundance',ylab=expression(paste('Modern distribution ( ',km^2,')',sep='')),col=cols,pch=16)
points(saveme1$adjustedBathy[seq(1,36,2)],saveme1$NOWas[seq(1,36,2)],pch=16,cex=2)
points(saveme1$adjustedBathy[seq(1,36,2)],saveme1$NOWas[seq(1,36,2)],pch=16,cex=.5,col=rainbow(18))

NOWavg<-saveme$NOWas
ABavg<-saveme$adjustedBathy
reg1<-lm(log(NOWavg)~ABavg)

xnew<-data.frame(ABavg=seq(.2,2.7,0.01))
preda<-predict(reg1,newdata=xnew,interval='confidence')
lines(xnew$ABavg,exp(preda[,'fit']),lwd=2)
lines(xnew$ABavg,exp(preda[,'lwr']),lwd=1,lty=2)
lines(xnew$ABavg,exp(preda[,'upr']),lwd=1,lty=2)

legend('topright',pch=16,col=rainbow(18),legend = splist,cex=.9)
dev.off()










SST_min_LGM



















##############################################




# #read IUCN 
# list1<-c(1,1,2,2,2,2,2,1,2,2,2,2)                       #species to choose from    ### important to change to the 6 most valuable species
# list2<-c(58,30,65,129,140,180,238,233,395,398,455,470)

saveme<-data.frame()
splist<-c()
for(vers in 1:18){
  #vers<-1  
  set<-list1[vers]
  bnm<-list2[vers]
  setwd("E:/Work for RA/IUCN shapes")
  corals2<-read.dbf(paste("CORALS",set,".dbf",sep=""))
  corals<-corals2[order(corals2$binomial),]
  Species<-paste(corals[bnm,2])
  SPP<-Species
  
  setwd("E:/Paleoreconstruction/stats Oct 11 2019")
  x<-read.csv(paste('paleostats',Species,".csv" ))
  x1<-x[c(-26,-27),]
  saveme<-rbind(saveme,x1)
  colnames(saveme)<-colnames(x1)
  splist[vers]<-SPP
}

#rownames(saveme)<-rep(splist,25)
#saveme<-as.data.frame(saveme)
saveme<-saveme[,-1]
colMeans(saveme)

skips<-which(saveme$LGMas==0)
saveme<-saveme[saveme$LGMas!=0,]



plot(saveme$adjustedBathy,saveme$NOWas)
plot(saveme$LGMas,saveme$adjustedBathy)
plot(saveme$NOWas,saveme$LGMas)
#plot(saveme$LGMas,saveme$NOWas)



cols<-(rainbow(18))
cols<-rep(cols,each=25)
#cols<-cols[-skips]




LGMavg<-saveme$LGMas/500000
NOWavg<-saveme$NOWas/500000


#summary(lm(log(LGMavg)~log(NOWavg)))

#a * (1.0 - exp(bx))
#reg1<-nls((LGMavg)~a*(1.0 - exp(b*NOWavg)),start=list(a=-0.17639 ,b=0.36856 ))

#Asym/(1 + exp(b + K*Age))
reg1<-nls((LGMavg)~Asym/(1.0 + exp(b+K*NOWavg)),start=list(K=-0.17639 ,b=0.36856,Asym=1.5))
#summary(reg1)

#reg1<-lm((LGMavg)~log(NOWavg))
xnew<-data.frame(NOWavg=seq(400000,3000000,100))
preda<-predict(reg1,newdata=xnew/500000,interval='confidence')


df <- 
  structure(list(x = c(NOWavg), 
                 y = c(LGMavg)), 
            .Names = c("x", "y"), row.names = c(NA, -448L),class = "data.frame")

m0<-reg1
df$pred <- predict(m0)
se = summary(m0)$sigma
ci = outer(df$pred, c(outer(se, c(-1,1), '*'))*1.96, '+')
ii = order(df$x)



splist[splist=="Favia speciosa"]<-"Dipsastraea speciosa"
splist[splist=='Favia pallida']<-'Dipsastraea pallida'

setwd("E:/Paleoreconstruction/Stats plots feb 2020")
tiff(paste("km2 modern vs LGM mar20.tif"),width=2300,height=1900, res = 300)

par(mar=c(5,5,2,1))
plot(saveme$NOWas,saveme$LGMas,xlab=expression(paste('Contemporary distribution (',km^2,')',sep='')),ylab=expression(paste('LGM distribution (',km^2,')',sep='')),col=cols,pch=16,ylim=range(ci*500000))
points(saveme1$NOWas[seq(1,36,2)],saveme1$LGMas[seq(1,36,2)],pch=16,cex=2)
points(saveme1$NOWas[seq(1,36,2)],saveme1$LGMas[seq(1,36,2)],pch=16,cex=.5,col=rainbow(18))
with(df[ii,], lines(x*500000, pred*500000, ylim=range(ci*500000), type='l'))
matlines(df[ii,'x']*500000, ci[ii,]*500000, lty=2, col=1)
legend('bottomright',pch=16,col=rainbow(18),legend = splist,cex=.8)
dev.off()


#lines((xnew$NOWavg),(preda*500000),lwd=2)
#lines((xnew$NOWavg),(preda[,'fit']),lwd=2)
#lines(xnew$NOWavg,exp(preda[,'lwr']),lwd=1,lty=2)
#lines(xnew$NOWavg,exp(preda[,'upr']),lwd=1,lty=2)

#legend('topright',pch=16,col=rainbow(12),legend = splist,cex=.6)


##############################################################
