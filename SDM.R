#Paleoreconstruction SDM

#changes:
#dispersal mask is impossible to run in reverse, especially 21000 years, assume dispersal was possible anywhere, except Caribbean/Pacific.
#turbidity would be different because of tidal amplitude on new shelves.... hard to fix. will run with resampled version of contemporary for now (just changes mask and resolution)

#can't mask rivers the same, sea-level was lower and masks would be in the wrong spot, also river volume was potentially much different.- #masking out rivers for contemporary, but unknown location for the past?, I can move based on google earth outflow...
#fixed, new locations,assume rivers have same flow rate. 

#assume no adaptation over 21000 years...

#there is >3x less habitat for LGM reefs (<30m) because of steepness of continental shelf, so analyses will only be good comparing %change in total habitat..

# Check Milankovitch cycles for changes in irradiance!!! PAR min, maybe wont be able to get PAR range. PARmin is finction of turb and PAR mean. DONE

# Save AUC, gam1e/slick

# NEED slick habitat before, after and corrected for habitat availability

#TO DO


#-will take ~4 days to run, put on other pc. still ties up external drive

library(blme)
library(beepr)
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
#Sys.setenv(JAVA_HOME='C:/Program Files (x86)/Java/jre1.8.0_201');loadNamespace('rJava')
library(mailR)
library(raster)
library(foreign)

load("E:/work for RA/Comps/load data comps 10-1-2014.RData") #load initial data

kd<-raster("E:/work for RA/Light Attenuation kd/mean kd 1998-2007 9km.nc ")#average light attenuation 490 coefficiten to 1998 to 1007 mean of each layer NA values removed. resolution = .0833333 (9km)

####################
####################
####################
####################
#################### SET Mechanistic values
######################################## these assumptions are more uncertain into the past, because of unknown adaptation rates


minTolerance<-18.0         #degrees C of minimum tolerance
dispersalDist<-10         #how far the species can disperse each year (KM)
LongitudinalBuffer<-5     #the initial longitudinal buffer for dispersal uncertainty values(0:10) numeric
ProjTime<-21000             #years until projection (LGM)
AdaptationRange<-0       #how much range will the species be able to tolerate in 2100 (can change to time relationship) - 1 degree gain in range tolerance
PAR.limit<-21.6 #value from Kleypas 1997, 250 microE/m^2/s = 21.6 E/m^2/day
Zmax<-3  #set  depth at which light calculation happens

ncores<-5 #set the number of cores for PC to run with
ATLvPAC<-1


########################################
####################
####################
####################

rasters #rasters has training data
rasters2 #rasters2 has predicitng contemporary data


#switch out future sst range for paleoSST range
#pastsstrange2<-raster('E:/paleoreconstruction/sstrangeLGMcc.nc')  #MARGO and Climap data on LGM 21000 bp
pastsstrange<-raster('E:/paleoreconstruction/CMIP5 LGM tos/SSTrangeLGM_CMIP5.nc')
Past<-Future
Past<-resample(Past,pastsstrange)
Past$SST_range_current<-pastsstrange

pastParRANGE<-raster('E:/paleoreconstruction/insolation/LGM PAR_range.nc')  #PAR range for LGM
pastParRANGE<-resample(pastParRANGE,pastsstrange)
Past$PAR_range<-pastParRANGE


pastPARmean<-raster('E:/paleoreconstruction/insolation/LGM PARmean.nc')  #PAR mean for LGM 

change.kd_to_PLD<-function(x,Zmax=3){ #for rasters
  x2<-.6677*x^.6763 #kdPAR~kd490 pierson et al. 2007
  Y2<-exp((-1*x2)*Zmax) #percent light at depth zmax
  return(Y2)
}

PLD<-change.kd_to_PLD(kd) #assume same turbidity...


pastPLD_mean<-pastPARmean*PLD
pastPLD_mean<-resample(pastPLD_mean,pastsstrange)
Past$PLD_mean<-pastPLD_mean


PAR.mean<-raster("E:/Work for RA/github compress/load data RA/SWFMO_PAR.CR.timeAverage1998-2007.nc") #used to determine effects of light
pPARmeanFINAL<-crop(PAR.mean,extent(PLD))*PLD #prediction raster

Past #has predicting LGM data

#set the prediction data to the past

#SST_min_LGM<-raster('E:/paleoreconstruction/sstminLGMcc.nc') 
SST_min_LGM<-raster('E:/paleoreconstruction/CMIP5 LGM tos/SSTminLGM_CMIP5.nc')
mask_deep_bathyLGM<-raster('E:/paleoreconstruction/sea-level/sealevelLGM130.nc')
riversBufferLGM<-readShapeSpatial('E:/paleoreconstruction/rivers buffer/riversBufferLGM.shp')
mask_deep_bathy<-raster('E:/paleoreconstruction/sea-level/modern 0-30m.nc')







#acropora hyacinthus = set 1 i<-58
#plobata = set<-2;i<-398

# cl <- makeCluster(ncores) # 
# registerDoParallel(cl)
# registerDoParallel(cl)

#read IUCN 
# list1<-c(1,1,2,2,2,2,2,1,2,2,2,2)                       #species to choose from    ### important to change to the 6 most valuable species
# list2<-c(58,30,65,129,140,180,238,233,395,398,455,470)


#add 6 intermediate species

#Dipsastaea pallida
#Galaxea fascicularis
#Platygyra sinensis
#Goniastrea retiformis
#Favites pentagona
#Acropora cerealis

list1<-c(1,2,2,2,2,2)   
list2<-c(19,125,143,176,191,331)


for (vers in c(1)){ #set to 1:12 for all spp
  #vers<-1    # this is the species choosing number, 10 is p.lobata, 1 is a hyacinthus
  set<-list1[vers]
  i<-list2[vers]
  
  setwd("E:/Work for RA/IUCN shapes")
  corals2<-read.dbf(paste("CORALS",set,".dbf",sep=""))
  corals<-corals2[order(corals2$binomial),]
  Species<-paste(corals[i,2])
 
   
  #convert to CRS
  setwd(paste("E:/Work for RA/IUCN shapes/Coral shapes/Corals",set,sep=""))
  IUCN_dist2<-readOGR(".",paste("CORALS",set,"_binomial__",Species,sep=''))
  IUCN_dist<-spTransform(IUCN_dist2,CRS(newproj))
  
  sparea<-raster::area(IUCN_dist)
  #find ecregions the species is in
  IUCNinECO2<-sort(as.numeric(as.character(over(IUCN_dist,ecoregions,returnList = T)[[1]]$SP_ID)))+1 
  IUCNinECO<-replace(IUCNinECO2, IUCNinECO2==142, 141)
  
  print(paste("IUCN data converted to ecoregions",Species))
  
  
  
  
  
  
  
  
  
  
  
  
  
  ########## probability ##################################################  Used in the first manuscript
  regions<-ereg1[IUCNinECO,2];nregions<-regions #regions the spp is in  
  one2141<-c(1:141)
  absent<- one2141[!one2141 %in% regions] #find the regions the species is not in
  
  
  adjlist<-data.frame() #list of all adjacent regions to the species (includes regions species is in)
  for (j in nregions){
    AdjRegions<-AdjacencyMat[j,1:10]
    adjlist<-rbind(adjlist,AdjRegions)}
  noNAadjlist<-adjlist[!is.na(adjlist)]
  No_Present_Or_NA_adjlist<-noNAadjlist[!noNAadjlist%in%(nregions)]######################### No_Present_Or_NA_adjlist is list of adjacent regions that do not have the species detected, not unique because we sample through the number of iterations
  
  ########## give each initial absent region a probability of 30% to find species on first try.
  problist<-data.frame()
  for (i in 1:141){
    if (i %in% nregions){ # if it is in a region P(detect)=1
      x<-1}
    else{x<-.3} # if not in region P(detect)= .3 
    problist<-rbind(problist,x)                          
  }
  ##################################### initial values set
  
  
  problist1<-c(do.call("cbind",problist))  ################ change .6 as decrease of adjacent regions
  for (i in No_Present_Or_NA_adjlist){
    problist1[i]<-problist1[i]*.6}  ########### everytime a species is detected adjacent to a region the species is absent from, probablility of detection is decreased by 40% because it is more likely the samplers missed it.
  
  ############### incorporate size into P
  
  Plist<-problist1
  for (i in absent){
    size<-sizes[i]
    ratioArea<-(minsize/size)^4 ##for regions the species is absent in take the log area of the region and divide log min area by it
    Plist[i]<-Plist[i]*ratioArea}
  problist<-Plist
  
  #atlantic and pacific vicariance - excluded in paleoreconstruction  
  # if (any(!regions %in% c(123:137))==T){ ## they are not in atlantic
  #   problist[c(123:137)]<-1 
  #   ATLvPAC<-1}#caribbean <- known absent
  # 
  # if (any(regions %in% c(123:137))==T){  ## they are in atlantic
  #   (problist[c(1:122,139:141)]<-1)
  #   ATLvPAC<-0} # pacific <- known absent
  # 
  #ATLvPAC is 1 if it is a pacific species
  
  probmat<-problist
  
  cl <- makeCluster(ncores) # 
  registerDoParallel(cl)
  
  #determine which regions the species is likely to be absent
  
  print("p2")
  variance<-foreach (i = 1:141,.combine=rbind,.verbose=F) %dopar% {  ##########################loop it  for variance and Prob of Geometric
    p<-probmat[i]
    nfails<-SEF[i]
    
    maxSE<-max(SEF) #find the maximum value of sampling effort
    
    #geometric distribution for SE out of maxSE
    pgm<-pgeom(nfails,p)#probablility of finding this spp after x trials 
    mx<-pgeom(maxSE,p)#prob of finding after max attempts
    PGM<-pgm /mx  #difference in probablility of finding after x trials - max attempts
    
    #determine variance
    VAR<-(1-p)/p^2  # get the variance
    EXP<-1/p    #expected value
    c(VAR,EXP,nfails,PGM) #combine variance, expected, sampling effort, probability of finding after SE trials
  }
  print("done var")
  
  n<-1
  se<-(sqrt(variance[,2]))/sqrt(n)
  p<-probmat
  #confidence of 95%
  c95<-se*1.96 # determine the 95% confidence of regions
  s95<-((log(1-.95)/log(1-p))-1) # determine the SE at which we are 95% sure we have found it.
  plotsp<-cbind(variance,c95,s95,p)
  rownames(plotsp)<-scip[,1]
  colnames(plotsp)<-(c("var","exp","nfails","ProbFound","c95","s95","p"))
  oplotsp<-plotsp[order(plotsp[,6]),]
  
  
  low_confidence_regions<-foreach (i = 1:141,.combine=rbind)%dopar%{
    if(plotsp[i,4]<.95){i} ####################################### only take ones greater than 95% (could include the confidence of that confidence?)
  }
  lowregions<-scip[low_confidence_regions,1]#if you need the region names
  
  
  stopCluster(cl)
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  ############# shapes of classified regions ############
  #uncertain
  n.regions<-as.numeric(low_confidence_regions[,1]) 
  n.shpregions<-ereg[n.regions,604]
  lowconf<-ecoregions[c(n.shpregions),]   
  
  #present
  p_regions<-regions
  shpregions<-ereg[regions,604]
  this_spp<-ecoregions[c(shpregions),]   
  this_spp<-gBuffer(this_spp,width=0,byid=T)
  
  #absent
  pres_and_uncert<-c(n.regions,p_regions)
  all<-c(1:141)
  no<-setdiff(all,pres_and_uncert)
  certainGone<-ereg[no,604]
  none<-ecoregions[c(certainGone),]  
  
  
  
  
  #sampling numbers
  if(ATLvPAC==1){absentSampleNumber<-(length(none)-15)}
  if(ATLvPAC==0){absentSampleNumber<-(length(none)-125)}
  absentSampleNumber<-abs(absentSampleNumber)
  absentSampleNumber/(length(this_spp)+absentSampleNumber)->ratio
  asamp<-round(3000*ratio);psamp<-round(3000-asamp)
  
  
  print(paste("probability is determined", length(regions),"present",length(no),"absent,","presence samples=",psamp,"absence samples=",asamp))
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  #run 25 times for confidence in projection
  for (CONF in 25:25){ #set to 1:25 for all runs
    
    
    
    print('start points')
    ############################## extract points
    if(ATLvPAC==1){Atldepth<-AtldepthP}
    if(ATLvPAC==0){Atldepth<-AtldepthC}
    
    cl <- makeCluster(ncores) # 
    registerDoParallel(cl)
    
    
    #absence points
    ecoA<-none
    ecoA<-spTransform(ecoA,CRS(proj4string(Atldepth)),force_ring=TRUE)
    ecoA<-gBuffer(ecoA,width=0,byid=F)
    Atldepth<-gBuffer(Atldepth,width=0,byid=F)
    Absintersection<-gIntersection(ecoA,Atldepth)
    
    cuts<-round((asamp+300)/25)
    cutsamp<-round(seq(1,asamp,length.out=cuts))
    
    
    abspts<-foreach(i=1:(cuts-1),.packages="sp",.combine=rbind) %dopar%{
      abspt<-spsample(Absintersection,n=(cutsamp[i+1]-cutsamp[i]), type='random',iter=100)
      coordinates(abspt)}
    
    abspts<-SpatialPoints(abspts)
    abspts.<-SpatialPointsDataFrame(abspts,data=data.frame(1:length(abspts)))
    
    #presence points
    ecoP<-this_spp
    ecoP<-spTransform(ecoP,CRS(proj4string(reefs)),force_ring=TRUE)
    
    
    
    
    
    
    
    
    ####################################################################
    list<-over(ecoP,reefs,returnList =T)
    short<-c(unique(unlist(list)))
    #short<-(list[!is.na(list)])
    inreefs<-reefs[short]   ####### REEFS INSIDE ECOREGIONS ################
    
    blure<-paste("determine presence points");print(blure)
    
    #   stopCluster(cl)
    #   
    #   cl <- makeCluster(ncores) # 
    #   registerDoParallel(cl)
    
    cuts<-round(psamp/10)
    cutsamp<-round(seq(1,psamp+200,length.out=cuts))
    scuts<-round(seq(1,length(inreefs),length.out=cuts))
    inreefst1<-foreach(i=2:cuts,.packages=c("sp"),.combine=rbind) %dopar%{
      lrfs<-(inreefs[scuts[i-1]:scuts[i]])
      rpoints<-spsample(lrfs,n=(cutsamp[i]-cutsamp[i-1]),"random",iter=100)
      coordinates(rpoints)
    }
    
    inreefst<-SpatialPoints(inreefst1)
    Prespoints1<-SpatialPointsDataFrame(inreefst,data=data.frame(1:length(inreefst)))
    proj4string(Prespoints1)<-proj4string(ecoP)
    testpoint<-(over(Prespoints1,as(ecoP, "SpatialPolygons")))
    tp2<-cbind(testpoint,1:length(testpoint))
    shlst<-tp2[complete.cases(tp2),2] # take only complete cases
    snum<-sample(shlst,psamp) #randomly sample the needed number of values
    Prespoints<-Prespoints1[snum,]
    
    try(stopCluster(cl))
    
    
    print(paste("done sampling points"))
    
    ########################
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    ################# Regression
    #################  regression
    #################  regression
    #################  regression
    
    proj4string(abspts.)<-newproj
    tabs_pts<-abspts.
    tpts<-Prespoints
    
    pts<-spTransform(tpts,CRS=CRS(nextproj))
    abs_pts<-spTransform(tabs_pts,CRS=CRS(nextproj))
    
    
    ############# CONNECTIVITY LAYER - mask outside 5 degrees for contemporary #############
    
    longs<-(coordinates(pts)[,1]) #lon of pts
    lats<-(coordinates(pts)[,2])
    rlongs<-round(longs)
    rlats<-round(lats)
    urlon<-unique(rlongs)
    urlat<-unique(rlats)
    turlon<-180+urlon
    
    
    today<-data.frame()
    for (i in 1:LongitudinalBuffer) {
      today<-c(today,(unique(c(turlon,turlon+i,turlon-i))))  
    }
    
    todayL<-data.frame()
    for (i in 1:LongitudinalBuffer) {
      todayL<-c(todayL,(unique(c(urlat,urlat+i,urlat-i))))  
    }
    
    
    unturlat<-as.numeric(unique(todayL))
    unturlon<-as.numeric(unique(today))
    xc1<-sub(361,1,unturlon);xc2<-sub(362,2,xc1);xc3<-sub(363,3,xc2);xc4<-sub(364,4,xc3);xc5<-sub(365,5,xc4);xc6<-sub(366,6,xc5);xc7<-sub(367,7,xc6);xc8<-sub(368,8,xc7);xc9<-sub(369,9,xc8);xc10<-sub(370,10,xc9)
    xc<-sub( "^0",360,xc10)
    xc12<-sub( "^-1$" ,359,xc);xc13<-sub( -2 ,358,xc12);xc14<-sub( -3 ,357,xc13);xc15<-sub( -4 ,356,xc14);xc16<-sub( -5 ,355,xc15);xc17<-sub( -6 ,354,xc16);xc18<-sub( -7 ,353,xc17);xc19<-sub( -8 ,352,xc18);xc20<-sub( -9 ,351,xc19);xc21<-sub( -10 ,350,xc20)
    
    xct<-as.numeric(xc21)
    uxct<-unique(xct)
    #min(uxct);max(uxct)
    
    mlat<-matrix(ncol=360,nrow=180)
    mlon<-matrix(ncol=360,nrow=180)
    mlon[,uxct]<-1
    mlat[,]<-1
    mlat[90-unturlat,]<-0
    q<-raster(mlon)
    q2<-raster(mlat)
    extent(q)<-c(xmin=-180,xmax=180,ymin=-90,ymax=90)
    q<-crop(q,extent(-180,180,-37,37))
    proj4string(q)<-"xy"
    proj4string(q)<-nextproj
    extent(q2)<-c(xmin=-180,xmax=180,ymin=-90,ymax=90)
    q2<-crop(q2,extent(-180,180,-37,37))
    proj4string(q2)<-"xy"
    proj4string(q2)<-nextproj
    
    
    q<-mask(q,q2,maskvalue=1)
    # plot(q)
    # plot(q2,add=T)
    # points(pts)
    
    #combine lat and lon masks
    q<-resample(q,SST_min_current)
    
    print(paste("Connectivity layer done"))
    
    ##############################################
    
    
    
    
    
    
    
    
    
    
    ####################### KFOLD for HIGHEST AUC ############################ here is first run of the model
    
    cl <- makeCluster(ncores) # 4 cores
    registerDoParallel(cl)
    
    
    
    perMutations<-100
    pscoreholder <- foreach(i = 1:perMutations, .packages = c("dismo","raster","stats","blme")) %dopar% {
      kfold(pts, 5) ##subset into 5 groups
    }
    ascoreholder <- foreach(i = 1:perMutations, .packages = c("dismo","raster","stats")) %dopar% {
      kfold(abs_pts, 5) ##subset into 5 groups 
    }
    AUC <- foreach(i = 1:perMutations, .packages = c("dismo","raster","stats","lme4"),.combine='c') %dopar% {
      ##subset into 5 groups
      
      pres_train <- pts[c(pscoreholder[[i]]) != 1, ]
      pres_test <- pts[c(pscoreholder[[i]]) == 1, ]
      
      abs_train <- abs_pts[c(ascoreholder[[i]]) != 1, ]
      abs_test <- abs_pts[c(ascoreholder[[i]]) == 1, ]
      
      #extract training values
      presvals1 <- extract(rasters2, pres_train,method="bilinear") #rasters has training data
      absvals1 <- extract(rasters2, abs_train,method="bilinear")  #rasters2 has predicitng contemporary data
      presvals<-round(presvals1,1)
      absvals<-round(absvals1,1)
      ##format data
      pb <- c(rep(1, nrow(presvals)), rep(0, nrow(absvals)))
      sdmdata <- data.frame(cbind(pb, rbind(presvals, absvals)))
      
      sdmdata<-sdmdata[complete.cases(sdmdata),]
      
      ####################### MODEL CODE!!!
      ###################
      ###############
      ##########
      ######
      try(gam1 <- bglmer(pb ~SST_range_current*PAR_range+(1|PLD_mean),data=sdmdata,cov.prior=gamma,verbose=T ,family = binomial(link = "logit"))) 
      ### AUC = true positives/false positives
      ######
      #########
      ##############
      ##################
      ########################
      
      
      
      ## set up evaluation data
      testpres <- data.frame( extract(rasters, pres_test) )
      testpres<-testpres[complete.cases(testpres),]
      testbackg <- data.frame( extract(rasters, abs_test) )
      testbackg<-testbackg[complete.cases(testbackg),]
      ##
      
      ### evaluate and set a threshold
      
      gam1e <- evaluate(testpres, testbackg, gam1,allow.new.levels=T)
      gam1e@auc
    }
    
    maxAUC<-which.max(AUC)
    
    pres_train <- pts[c(pscoreholder[[maxAUC]]) != 1, ]
    pres_test <- pts[c(pscoreholder[[maxAUC]]) == 1, ]
    
    abs_train <- abs_pts[c(ascoreholder[[maxAUC]]) != 1, ]
    abs_test <- abs_pts[c(ascoreholder[[maxAUC]]) == 1, ]
    
    #extract training values
    presvals <- extract(rasters2, pres_train,method="bilinear") 
    absvals <- extract(rasters2, abs_train,method="bilinear")
    
    ##format data
    pb <- c(rep(1, nrow(presvals)), rep(0, nrow(absvals)))
    sdmdata <- data.frame(cbind(pb, rbind(presvals, absvals)))
    sdmdata[,4]<-round((sdmdata[,4])*2,-1)/2   ####  sdmdata[,4]<-round(sdmdata[,4],1)  sdmdata[,4]<-round((sdmdata[,4]*2),0)/2  
    
    
    print(paste("Kfold finished"))
    
    
    
   
    
    
    
    ####################### RUN initail model with  highest AUC   #######
    #################
    ############
    #######
    ####
 gam1<-foreach(i = 1:1, .packages = c("dismo","raster","stats","lme4"),.combine='c') %dopar% {
    gam1 <- bglmer(pb ~SST_range_current*PAR_range+(1|PLD_mean),data=sdmdata, family = binomial(link = "logit"),verbose=T,cov.prior=gamma)#,control=glmerControl(optCtrl=list(maxfun=100000) ))
    gam1
    }
    
    stopCluster(cl)
    
    ####
    #######
    ###########
    ###########
    ################
    ####################
    ### AUC = true positives/false positives
    ## set up evaluation data
    testpres <- data.frame( extract(rasters, pres_test) )
    testpres<-testpres[complete.cases(testpres),]
    testbackg <- data.frame( extract(rasters, abs_test) )
    testbackg<-testbackg[complete.cases(testbackg),]
    enutriTPRES<-data.frame( extract(nutrients, pres_test) )
    enutriTPRES<-enutriTPRES[complete.cases(testpres),]
    enutriTBACKG<-data.frame( extract(nutrients, abs_test) ) 
    enutriTBACKG<-enutriTBACKG[complete.cases(testpres),]
    ##
    testinga<-rbind(testpres,testbackg)
    ### evaluate and set a threshold
    testing1<-raster(ncol=1,nrow=sum(nrow(testpres),nrow(testbackg)));values(testing1)<-testinga[,1]
    testing2<-raster(ncol=1,nrow=sum(nrow(testpres),nrow(testbackg)));values(testing2)<-testinga[,2]
    testing3<-raster(ncol=1,nrow=sum(nrow(testpres),nrow(testbackg)));values(testing3)<-testinga[,3]
    testing<-stack(testing1,testing2,testing3)
    names(testing)<-names(testinga)
    
    tstd.a <- predict(testing,gam1,type="response",na.action=na.omit,re.form=NULL,allow.new.levels=T)
    #tstd<-tstd.a-(.0432*(enutriTPRES*(2.86297283611898e+13/sparea))* 1.73819020691208e-07)
    evald<-evaluate(values(tstd.a)[1:nrow(testpres)],values(tstd.a)[nrow(testpres):length(values(tstd.a))])
    
    
    tr<-threshold(evald, 'spec_sens')
    
    gam1e <- evaluate(testpres, testbackg, gam1,re.form=NULL,allow.new.levels=T)
    gam1e@auc
    
    print('predict')
    pgam1.a <- predict(rasters2,gam1,type="response",na.action=na.omit,re.form=NULL,allow.new.levels=T)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    # 
    # 
    # ################### FUTURE MODEL RUN #######################
    pgam2.a <- predict(Past, gam1,type="response",na.action=na.omit,re.form=NULL,allow.new.levels=T)
    # pgam2B1 <- predict(FutureB1, gam1,type="response",na.action=na.omit,re.form=NULL,allow.new.levels=T)
    # pgam2A1B <- predict(FutureA1B, gam1,type="response",na.action=na.omit,re.form=NULL,allow.new.levels=T)
    # 
    # 
    # #################### ADAPTATION MODEL RUN ###################
    # FTA2<-TempRange2100x-AdaptationRange
    # FTAB1<-TempRangeB1x-AdaptationRange
    # FTAA1B<-TempRangeA1Bx-AdaptationRange
    # 
    # adapRasA2<-stack(PAR_rangex,FTA2,pPARmeanFINAL)
    # adapRasB1<-stack(PAR_rangex,FTAB1,pPARmeanFINAL)
    # adapRasA1B<-stack(PAR_rangex,FTAA1B,pPARmeanFINAL)
    # 
    # 
    # 
    # names(adapRasA2)<-c('PAR_range','SST_range_current','PLD_mean')
    # adapRasA2<-crop(adapRasA2,extent(-180,180,-37,37))
    # pgam3 <- predict(adapRasA2, gam1,type="response",na.action=na.omit,re.form=NULL,allow.new.levels=T)
    # 
    # names(adapRasB1)<-c('PAR_range','SST_range_current','PLD_mean')
    # adapRasB1<-crop(adapRasB1,extent(-180,180,-37,37))
    # pgam3B1 <- predict(adapRasB1, gam1,type="response",na.action=na.omit,re.form=NULL,allow.new.levels=T)
    # 
    # names(adapRasA1B)<-c('PAR_range','SST_range_current','PLD_mean')
    # adapRasA1B<-crop(adapRasA1B,extent(-180,180,-37,37))
    # pgam3A1B <- predict(adapRasA1B, gam1,type="response",na.action=na.omit,re.form=NULL,allow.new.levels=T)
    # 
    # print(paste("Contemporary, Future, and adaptation models run"))
    # 
    # 
    # 
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    #####################
    ####################  NUTRIENT EFFECT, bases on endinger Plobata study inference
    
    #pgam1<-pgam1.a-(.0432*(nutrients*(2.86297283611898e+13/sparea))* 1.73819020691208e-07)
    #pgam2<-pgam2.a-(.0432*(nutrients*(2.86297283611898e+13/sparea))* 1.73819020691208e-07)
    
    
    ######## Nutrients have no effect, to discard them from this analysis use
    
    pgam1<-pgam1.a
    pgam2<-pgam2.a
    
    
    
    ######################
    ##################
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    ############### First Masking ################
    
    ## mask out reef, temp<18, and longitude
    extent(mask_deep_bathy)<-c(-180,180,-37,37)
    print('masking 1')
    maskreef<-mask(pgam1,mask_deep_bathy,maskvalue=NA) #mask less than 20 meters
    mask_temp_reef<-mask(maskreef,SST_min_current>minTolerance,maskvalue=0) #mask less than the minimum tolerance degrees C
    mask_temp_reef_lon<-mask(mask_temp_reef,q,maskvalue=NA) #mask outside distribution
    mask_temp_reef_lon<-mask(mask_temp_reef_lon,riversBuffer,inverse=T) #mask rivers
    mask_temp_reef_lon<-mask(mask_temp_reef_lon,pPARmeanFINAL<PAR.limit,maskvalue=1)#mask too low light
    mask_temp_reef_lon_ATL<-mask(mask_temp_reef_lon,ATL,inverse=T) #mask the Atlantic
    
    
    
    
    
    
    
    
    ##################### DISPERSAL BETWEEN TIME PERIODS #################### - not in paleoreconstruction
    # print('start dispersal')
    # initdistr<-(mask_temp_reef_lon_ATL > tr)
    # initdistr[initdistr==0] <- NA
    # dispdistKM<-dispersalDist*ProjTime
    # 
    # 
    # lx<-5330-dispersalDist*333 ## random points to be taken from the distribution, should cover all # using this formula, 1km samples 5000 pts, and 10 samples 2000, ranging linearly between
    # prand<-randomPoints(initdistr,lx)  #take random background points to cover extent #transform CRS for projection
    # 
    # sprand<-SpatialPoints(prand)
    # proj4string(sprand) <- CRS(nextproj)
    # lox<-spTransform(sprand,CRS(newproj))
    # 
    # IDK<-gBuffer(lox,width=(dispdistKM*1000),byid=F) ## buffer by the distance of dispersal
    # 
    # #polygon to raster
    # sdxprand<-SpatialPolygonsDataFrame(IDK,data=data.frame(1),match.ID=F)
    # idk<-vect2rast(sdxprand,cell.size=70000)
    # rasBuf<-(raster(idk))
    # reprojRB2<-projectRaster(rasBuf,crs=CRS(nextproj))
    # reprojRB3<-reprojRB2>0
    # reprojRB4<-crop(reprojRB3,extent(pgam1))
    # reprojRB<-resample(reprojRB4,pgam1)
    # 
    # print(paste("DISPERSAL BETWEEN TIME PERIODS done"))
    
    
    
    
    
    
    
    
    
    
    
    
    
    ################# MODEL ACCURACY #####################
    presconfusion<-data.frame(extract(pgam1,pts));names(presconfusion)<-'x'
    absconfusion<-data.frame(extract(pgam1,abs_pts));names(absconfusion)<-"x"
    conmat<-rbind(presconfusion,absconfusion)
    pi.hat <- exp(conmat)/(1 + exp(conmat))
    pball <- c(rep(1, nrow(pts)), rep(0, nrow(abs_pts)))
    flop<-cbind(pball,pi.hat)
    flopx<-flop[complete.cases(flop),]
    Lost<-confusion.matrix(flopx$pball,flopx$x,threshold=.5)
    ModelAccuracy<-(Lost[1,1]+Lost[2,2])/(Lost[1,1]+Lost[1,2]+Lost[2,1]+Lost[2,2])
    OM<-omission(Lost)
    Sens<-sensitivity(Lost)
    Specifi<-specificity(Lost)
    PropC<-prop.correct(Lost)
    h1<-paste(ModelAccuracy,"Model Accuracy")
    h2<-paste(OM,"omission")
    h3<-paste(Sens,"sensitivity")
    h4<-paste(Specifi,"specificity")
    h5<-paste(PropC,"prop.correct")
    
    
    
    print(paste("Model accuracy done"))
    
    
    
    # ## mask out reef, temp<18, and longitude
    # maskreef<-mask(pgam1,mask_deep_bathy,maskvalue=0) #mask less than 20 meters
    # mask_temp_reef<-mask(maskreef,SST_min_current>minTolerance,maskvalue=0) #mask less than the minimum tolerance degrees C
    # mask_temp_reef_lon<-mask(mask_temp_reef,q,maskvalue=NA) #mask outside distribution
    # mask_temp_reef_lon<-mask(mask_temp_reef_lon,riversBuffer,inverse=T) #mask rivers
    # mask_temp_reef_lon_ATL<-mask(mask_temp_reef_lon,ATL,inverse=T) #mask the Atlantic
    
    
    ################# MASKING for PLOTTING2 ####################
    #masking for  future scenarios and adaptation
    
    # cl <- makeCluster(ncores) # 
    # registerDoParallel(cl)
    # namess<-c("2","2B1","2A1B","3","3B1","3A1B")
    # 
    # listedras<-pgam2
    # for(str in 2:6){
    #   listedras<-stack(listedras,get(paste("pgam",namess[str],sep='')))}
    # names(listedras)<-
    #   
    #   rastermask<-foreach (prep=1:6,.packages=c('raster','sp')) %dopar%{ 
    #     maskreef2A1B<-mask(subset(listedras,prep),mask_deep_bathy,maskvalue=0)
    #     #### mask out all future places where min is less than 18 degrees
    #     mask_temp_reef2A1B<-mask(maskreef2A1B,SST_min_future>minTolerance,maskvalue=0) #mask too low of temperatures
    #     mask_temp_reef_lon2A1B<-mask(mask_temp_reef2A1B,reprojRB>0,maskvalue=NA) #mask distribution
    #     mask_temp_reef_lon2A1B<-mask(mask_temp_reef_lon2A1B,riversBuffer,inverse=T) #mamsk rivers
    #     mask(mask_temp_reef_lon2A1B,ATL,inverse=T) #mask atlantic ocean
    #   }
    
    # 
    # ### assign the names of the masked rasters
    # for (fg in 1:6){
    #   assign(paste("mask_temp_reef_lon_ATL",namess[fg],sep=''),rastermask[[fg]])
    # }
    
    
    #dissaggregate the output to mask
    
    
     
  
    #### mask out all HISTORIC places where min is less than 18 degrees
    mask_temp_reef_lon2<-mask(pgam2,SST_min_LGM>minTolerance,maskvalue=0) #mask too low of temperatures in LGM
    
    #mask_temp_reef_lon2<-mask(mask_temp_reef2,reprojRB>0,maskvalue=NA) #mask future distribution, no way?
    
    mask_temp_reef_lon2<-mask(mask_temp_reef_lon2,riversBufferLGM,inverse=T) #mask LGM rivers
    
    mask_deep_bathyLGM<-resample(mask_deep_bathyLGM,pPARmeanFINAL)
    mask_temp_reef_lon2<-resample(mask_temp_reef_lon2,mask_deep_bathyLGM)
    pastPLD_mean<-resample(pastPLD_mean,mask_temp_reef_lon2)
    
    mask_temp_reef_lon2<-mask(mask_temp_reef_lon2,pastPLD_mean<PAR.limit,maskvalue=1) #mask too low light in LGM
    
    mask_temp_reef_lon_ATL2<-mask(mask_temp_reef_lon2,mask_deep_bathyLGM,maskvalue=0) #mask the depth #must be 130m lower- 160
    
    mask_temp_reef_lon_ATL2<- mask(mask_temp_reef_lon_ATL2,ATL,inverse=T) #mask atlantic ocean
    
    print(paste("masking 2 done"))
    # proj4string(pgam2)<-CRS("+proj=longlat")
    # proj4string(pgam2B1)<-CRS("+proj=longlat")
    # proj4string(pgam2A1B)<-CRS("+proj=longlat")
    # proj4string(pgam3)<-CRS("+proj=longlat")
    # proj4string(pgam3B1)<-CRS("+proj=longlat")
    # proj4string(pgam3A1B)<-CRS("+proj=longlat")
    proj4string(mask_temp_reef_lon_ATL)<-CRS("+proj=longlat")
    
    
    #### RasterVis
    #of all area
    
    # mask_temp1<-mask(pgam1>tr,SST_min_future>minTolerance,maskvalue=0)
    # mask_temp_lon1<-mask(mask_temp1,q,maskvalue=NA)
    # mask_temp2<-mask(pgam2>tr,SST_min_future>minTolerance,maskvalue=0)
    # mask_temp_lon2<-mask(mask_temp2,q,maskvalue=NA)
    # dif<-mask_temp_lon1-mask_temp_lon2
    # loss_of_hab<-mask(dif,dif,maskvalue=-1)
    # loss_of_hab<-mask(loss_of_hab,ATL,inverse=T)
    # loss_of_hab<-mask(loss_of_hab,riversBuffer,inverse=T)
    # #of reefs
    # lostras<-((mask_temp_reef_lon_ATL>tr)-(mask_temp_reef_lon_ATL2>tr))>0
    # lostrass<-mask(lostras,lostras,maskvalue=0)
    # 
    
    
    ####################### AREA STATISTICS #######################
    #### find area of distribution #deepmask needs to be second
    
    # slick<-matrix(ncol=1,nrow=5,data=NA,dimnames=list(c("current","lost",'future','new','change'),c("2")))
    # strings<-c("2")#,"3","3B1","3A1B") #2100 A2,B1,A1B,  Adaptation A2,Adap B1,Adap A1B
    # for (p in 1:1){
    #   i<-strings[p]
    # 
    #   a<-raster::area(mask_temp_reef_lon_ATL)## the area of cells in meters^2? MAKE SURE MASS IS DETACHED
    #   #deepmask<-mask(a,mask_temp_reef_lon_ATL2)
    #   #assign(paste0("deepmask",i),mask(a,get(paste("mask_temp_reef_lon_ATL",i,sep=''))>tr))
    # 
    #   #current area
    #   curhab1<-mask(a,mask_temp_reef_lon_ATL>tr)
    #   curhab<-mask(curhab1,mask_temp_reef_lon_ATL>tr,maskvalue=0)
    #   curhabval<-(cellStats(curhab,sum)) #-66120
    #   slick[1,p]<-curhabval
    # 
    # 
    #   #lost
    #   pizza<-sum((get(paste("mask_temp_reef_lon_ATL",i,sep=''))>tr),(-2*(mask_temp_reef_lon_ATL>tr)),na.rm=F)
    #   zz2<-mask(a,pizza==(-2),maskvalue=0)
    #   assign(paste0("losthab",i),mask(zz2,pizza==(-2),maskvalue=NA))
    #   assign(paste0("losthabval",i),cellStats(get(paste("losthab",i,sep='')),sum))
    #   slick[2,p]<-get(paste("losthabval",i,sep=''))
    # 
    #   #future area
    #   zz4<-mask(a,(get(paste("mask_temp_reef_lon_ATL",i,sep=''))>tr),maskvalue=0)
    #   assign(paste0("futarea",i),mask(zz4,(get(paste("mask_temp_reef_lon_ATL",i,sep=''))>tr),maskvalue=NA))
    #   assign(paste0("futareaval",i),cellStats(get(paste("futarea",i,sep="")),sum))
    #   slick[3,p]<-get(paste("futareaval",i,sep=''))
    # 
    # 
    #   #new habitat gained
    #   assign(paste0("gain",i),mask((get(paste("futarea",i,sep=''))),curhab,inverse=T))
    #   assign(paste0("gainval",i),cellStats(get(paste("gain",i,sep='')),sum))
    #   slick[4,p]<-get(paste("gainval",i,sep=''))
    # 
    # 
    #   #percent changed
    #   assign(paste0("habchange",i),-(1-(get(paste("futareaval",i,sep=''))/curhabval))) # change in habitat
    #   slick[5,p]<-get(paste("habchange",i,sep=''))
    # 
    #   print(i)
    # }
    # 
    # slick<-rbind(slick,c(gam1e@auc));rownames(slick)[6]<-"auc"
    # slick<-rbind(slick,c(Species));rownames(slick)[7]<-"species"
    # slick<-rbind(slick,c(paste(Sys.Date())));rownames(slick)[8]<-"time"
    # 
    # 
    
    slick<-matrix(ncol=1,nrow=5,data=NA,dimnames=list(c("current area","gained since LGM",'LGM area','LGM only hab','pptn change')))
    for (i in 2){
      p<-1
      
      a<-raster::area(mask_temp_reef_lon_ATL)## the area of cells in meters^2? MAKE SURE MASS IS DETACHED
      #deepmask<-mask(a,mask_temp_reef_lon_ATL2)
      #assign(paste0("deepmask",i),mask(a,get(paste("mask_temp_reef_lon_ATL",i,sep=''))>tr))
      
      #current area
      curhab1<-mask(a,mask_temp_reef_lon_ATL>tr)
      curhab<-mask(curhab1,mask_temp_reef_lon_ATL>tr,maskvalue=0)
      curhabval<-(cellStats(curhab,sum)) #-66120
      slick[1,p]<-curhabval
      
      
      #lost
      pizza<-sum((get(paste("mask_temp_reef_lon_ATL",i,sep=''))>tr),(-2*(mask_temp_reef_lon_ATL>tr)),na.rm=F)
      zz2<-mask(a,pizza==(-2),maskvalue=0)
      assign(paste0("losthab",i),mask(zz2,pizza==(-2),maskvalue=NA)) 
      assign(paste0("losthabval",i),cellStats(get(paste("losthab",i,sep='')),sum))
      slick[2,p]<-get(paste("losthabval",i,sep=''))
      
      #future area
      zz4<-mask(a,(get(paste("mask_temp_reef_lon_ATL",i,sep=''))>tr),maskvalue=0)
      assign(paste0("futarea",i),mask(zz4,(get(paste("mask_temp_reef_lon_ATL",i,sep=''))>tr),maskvalue=NA))  
      assign(paste0("futareaval",i),cellStats(get(paste("futarea",i,sep="")),sum))
      slick[3,p]<-get(paste("futareaval",i,sep=''))
      
      
      #new habitat gained
      assign(paste0("gain",i),mask((get(paste("futarea",i,sep=''))),curhab,inverse=T))
      assign(paste0("gainval",i),cellStats(get(paste("gain",i,sep='')),sum))
      slick[4,p]<-get(paste("gainval",i,sep=''))
      
      
      #proportion changed
      assign(paste0("habchange",i),-(1-(get(paste("futareaval",i,sep=''))/curhabval))) # change in habitat
      slick[5,p]<-get(paste("habchange",i,sep=''))
      
      print(i)
    }
    
    
    slick<-rbind(slick,c(gam1e@auc));rownames(slick)[6]<-"auc"
    slick<-rbind(slick,c(Species));rownames(slick)[7]<-"species"
    slick<-rbind(slick,c(paste(Sys.Date())));rownames(slick)[8]<-"time"
    #slick<-rbind(slick,c(paste(print(VarCorr(gam1),comp=c("Variance")))));rownames(slick)[9]<-"Variance"
    slick<-rbind(slick,c(paste(sqrt(diag(VarCorr(gam1)$PLD_mean)))));rownames(slick)[9]<-"Std.Dev."
    
    print('slick')
    
    
    
    
    
    
    
    
    
    if(CONF==1){
      
      valeus<-(zonal(losthab2>0, init(losthab2, v='row'), fun='sum'))
      pagain<-gain2>0
      gains<-(zonal(pagain,init(pagain, v='row'), fun='sum'))
      
      lonplot<-40 #what lattitude you want the map to show
      btwnlons<-lonplot-37
      
      kmlost2<-valeus[,2]*9.3^2 #9.3*9.3 per raster cell
      kmlost<-c(rep(0,12*btwnlons),kmlost2,rep(0,12*btwnlons))
      rlos<-length(kmlost) #the length for the derivplot
      zmean<-c(0,0,0)
      for (i in 3:(length(kmlost)-4)){
        p<-mean(kmlost[c(i-3,i-2,i-1,i,i+1,i+2,i+3)])
        zmean<-c(zmean,p)
      }
      finmena<-c(zmean,0,0,0)
      
      
      kmgain2<-gains[,2]*9.3^2
      kmgain<-c(rep(0,12*btwnlons),kmgain2,rep(0,12*btwnlons))
      zmean<-c(0,0,0)
      for (i in 3:(length(kmlost)-4)){
        p<-mean(kmgain[c(i-3,i-2,i-1,i,i+1,i+2,i+3)])
        zmean<-c(zmean,p)
      }
      fingain<-c(zmean,0,0,0)
      
      
      
      
      closeer3<-cbind(finmena,fingain)
      closeer2<-cbind(closeer3,rowSums(closeer3));colnames(closeer2)[3]<-'sum'
      closeer1<-cbind(closeer2,max(closeer2[,3])-closeer2[,2]);colnames(closeer1)[4]<-'g second'
      
      
      
      
      
      
      
      tmask_temp_reef_lon_ATL2<-mask_temp_reef_lon_ATL2
      tmask_temp_reef_lon_ATL<-mask_temp_reef_lon_ATL
      tlosthab2<-losthab2
      
      extent(tmask_temp_reef_lon_ATL2)<-c(0,360,-37,37)
      tmask_temp_reef_lon_ATL2<-raster::rotate(tmask_temp_reef_lon_ATL2)
      extent(tmask_temp_reef_lon_ATL2)<-c(0,360,-37,37)
      
      extent(tmask_temp_reef_lon_ATL)<-c(0,360,-37,37)
      tmask_temp_reef_lon_ATL<-raster::rotate(tmask_temp_reef_lon_ATL)
      extent(tmask_temp_reef_lon_ATL)<-c(0,360,-37,37)
      
      extent(tlosthab2)<-c(0,360,-37,37)
      tlosthab2<-raster::rotate(tlosthab2)
      extent(tlosthab2)<-c(0,360,-37,37)
      
    }
    
    ####################save load data comps as an workspace
    
    print('saving')
    #Save data
    setwd("E:/paleoreconstruction/paleoreconstruction")
    newdir<-getwd()
    subDir<-Species
    dir.create(file.path(newdir, subDir), showWarnings = FALSE)
    setwd(file.path(newdir, subDir))
    
    if(CONF==1){
      save.image(paste(paste(getwd(),Species,sep="/"),".RData",sep=""))  #saves a workspace file for later use
    }
    
    #save raster
    writeRaster((mask_temp_reef_lon_ATL > tr),paste(Species,CONF,"initial dist.nc"),overwrite=T)
    writeRaster((mask_temp_reef_lon_ATL2 > tr),paste(Species,CONF, "LGM","distribution.nc"),overwrite=T)
    
    #save stats
    write.csv(slick,paste(Species,"model stats",'.csv'))
    if(CONF==1){ write.table(t(slick),paste(getwd(),"/",Species," table.csv",sep=''),sep=",",row.names = F,col.names=row.names(slick),append=T)}
    if(CONF!=1){ write.table(t(slick),paste(getwd(),"/",Species," table.csv",sep=''),sep=",",row.names =F,col.names=F,append=T)}
    
    if(CONF==1){
    zz <- file(paste(Species,"stats auto.txt"),"w")
    sink(zz)

    print(paste("INPUT"))
    print(paste("#degrees C of minimum tolerance="  ,minTolerance))
    print(paste("#how far the species can disperse each year KM =",dispersalDist))
    print(paste("#the initial longitudinal buffer for dispersal uncertainty (degrees) =",LongitudinalBuffer) )
    print(paste("years until projected time =",ProjTime))
    print(paste("adaptation coef", AdaptationRange))
    print(paste("    ADAPTATION VALUES"))
    print(gam1e)
    print(summary(gam1))
    print(Lost)
    cat(paste(h1,h2,h3,h4,h5, sep='\n')   )
    sink();close(zz)
    }
    
    
    
    
    
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
    
    #load("C:/Users/Chris/Desktop/Comps/load data comps 10-1-2014.RData")
    
    
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
    
    
    
    tmask_temp_reef_lon_ATL2<-rotater(mask_temp_reef_lon_ATL2>tr,center=210)
    extent(tmask_temp_reef_lon_ATL2)<-c(-180,180,-37,37)
    tmask_temp_reef_lon_ATL<-rotater(mask_temp_reef_lon_ATL>tr,center=210)
    extent(tmask_temp_reef_lon_ATL)<-c(-180,180,-37,37)
    tlosthab2<-rotater(losthab2>0,center=210)
    extent(tlosthab2)<-c(-180,180,-37,37)
    
    LGMshelf<-raster('E:/paleoreconstruction/sea-level/LGMshelf.nc')
    extent(LGMshelf)<-c(c(-180,180,-70,70))
    tLGMshelf<-rotater(LGMshelf>0,center=210)
    extent(tLGMshelf)<-c(-180,180,-70,70)
    tLGMshelf<-crop(tLGMshelf,c(-180,180,-45,45))
    tLGMshelf[tLGMshelf==1]<-NA
    ################# FIGURES
    
    if(CONF==1){
      
      
      tiff(paste(Species,"synthesis plot2.png"),width=3900,height=1500, res = 300)
      plot.map("world", transf=F, center=210 , col="burlywood",bg="white",ylim=c(-45,45),fill=TRUE,mar=c(2,5,2,0),add=F) #center is still 0
      rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col ="lightsteelblue1")
      
      image(tLGMshelf,add=T,col=rev(c('lightsteelblue1','khaki1')),legend=F,maxpixels=3836160)
      
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
      image((tmask_temp_reef_lon_ATL2 > tr),add=T,maxpixels=3836160,col=rev(c("#00000000","red3")),legend=F)  ### new habitat 
      image((tmask_temp_reef_lon_ATL > tr),maxpixels=3836160,add=T,col=c("#00000000","blue"),legend=F)  ### currently inhabited
      #image((losthab2>0),add=T,col=c("red"),maxpixels=3836160)
      
      #color.legend(140,28,175,32,legend=round(seq(-845,336,length.out=3)), rect.col=c(color.palette4,color.palette3),cex=.975)
      legend("topright",c('LGM shelf',"LGM coral","modern coral"),col=c('khaki1',"red3","blue"),pch=c(16,16,16),box.col=NA,bg="#00000000",cex=1.1)
      
      #mtext("mm above sea-level", 1, line=-12.15, adj=0.99,cex=.975)
      
      
      text(-72.31,-25.336,'Indian Ocean',cex=1.2)
      text(87.3,3.4,'Pacific Ocean',cex=1.2)
      scalebar(d=36,xy=c(66,-39),label=c(0,'',4000),cex=.9,type='bar',divs=4,below="kilometers",adj=c(0.5,-1.1))
      dev.off()
      
      
      # png(paste(Species,"synthesis plot.png"),width=1323,height=418)#1323,418
      # 
      # 
      # layout( matrix( c(1,1,1,1,2,2),ncol=3),widths=c(2,2,1) )
      # par(mar = c(5, 5, 6, 0))
      # image(land,col=c("white","burlywood"),xlab="Longitude",ylab="Latitude",cex.main=4,cex.lab=1.8,cex.axis=1.2,ylim=c(-lonplot,lonplot),xlim=c(20,290))
      # image((tmask_temp_reef_lon_ATL2 > tr),add=T,col=rev(c("#00000000","limegreen")),legend=F)  ### new habitat 
      # image((tmask_temp_reef_lon_ATL > tr),add=T,col=c("#00000000","blue"),legend=F)  ### currently inhabited
      # image((losthab2>0),add=T,col=c("red"))
      # legend("topright",c("lost habitat","gained habitat","retained habitat"),col=c("red","green","blue"),pch=c(16,16,16),box.col=NA,bg="#00000000",cex=1.4)
      # box('plot',lty='solid')
      # compassRose(345,-30,cex=.5)
      # axis(4,tick=T,labels=F,tck=.02)
      # 
      # 
      # par(mar=c(5,0,6,5),font.lab=1,font.main=1)
      # 
      # plot(finmena,1:length(finmena),type="l",ylim=rev(range(c(1,rlos))),xlim=c(0,max(closeer1[,3])),ann = FALSE,axes=F)
      # polygon(finmena,1:length(finmena),col=rgb(.8,.1,.1,alpha=1))
      # polygon(closeer1[,4],1:length(finmena),col=rgb(.0,.8,.2,alpha=.5))
      # lines(closeer1[,4],1:length(finmena))
      # axis(1,cex.axis=1.2,at=(round(seq(0,max(closeer1[,3]),length.out=8),-3)))
      # axis(3,cex.axis=1.2,at=(seq(0,max(closeer1[,3]),length.out=8)),labels=c(rev(paste(round(seq(0,max(closeer1[,3]),length.out=8),-3)))))
      # axis(4,at=c(0,length(finmena)*.25,length(finmena)/2,length(finmena)*.75,length(finmena)),labels=c("40","20","0","-20","-40"))
      # title(xlab=expression(paste("Habitat lost (km"^"2",")")),cex.lab=1.8)
      # title(main=expression(paste(" Habitat gained (km"^"2",")")),cex.main=1.8,xpd=T)
      # legend(0,1070,"",pch=5,xpd=T,col=rgb(.8,.1,.1,alpha=.5),bg=rgb(.8,.1,.1,alpha=1),cex=.6)
      # legend(0,-140,"",pch=5,xpd=T,col=rgb(.0,.8,.2,alpha=.5),bg=rgb(.0,.8,.2,alpha=.5),cex=.6)
      # mtext("Latitude",4,cex=1.3,line=3)
      # dev.off()
    }
    
    
    ###### notify me of species being completed
    
    

    
  }
  
  try(send.mail(from = "sender@gmail.com",to = c("slipping24@gmail.com"),subject = paste(Species,CONF,"finished a spp"),body = "<html> Duck - <img src=\"http://beehivehairdresser.com/wp-content/uploads/2010/03/duck-butt.jpg\"></html>",html = TRUE, smtp = list(host.name = "smtp.gmail.com", port = 465, user.name = "slipping24", passwd = "?rocket24", ssl = TRUE),authenticate = TRUE,send = TRUE))
  
  print(Sys.time())
  
 # play(sin(1:20000/17)*2)
  
}


try(stopCluster(cl))
try(emailME<-send.mail(from = "sender@gmail.com",to = c("slipping24@gmail.com"),subject = "error in R home",body = "<html> Duck - <img src=\"http://beehivehairdresser.com/wp-content/uploads/2010/03/duck-butt.jpg\"></html>",html = TRUE, smtp = list(host.name = "smtp.gmail.com", port = 465, user.name = "slipping24", passwd = "?rocket24", ssl = TRUE),authenticate = TRUE,send = TRUE))

