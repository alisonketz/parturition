###
### 1/31/2018 
### Alison Ketz
### Movement data preprocessing using adehabitatLT package
###


###
### Preliminaries
###

rm(list=ls())

library(e1071)
library(geosphere)
library(lubridate)
library(Hmisc)
library(zoo)
library(stringr)
library(ggplot2)
library(dplyr)
library(ggmap)
library(adehabitatHR)
library(adehabitatLT)
library(maptools)
library(changepoint)
library(sp)
library(spatstat)#for "duplicate" function
library(readr)
library(mvtnorm)



# setwd("~/Documents/Parturition/180108_parturition")


###
### Load VIT data
###
###
### For running on windows machine
###

library(RODBC)
setwd("F:/Parturition/180108_parturition")

myconn=odbcConnectAccess('C:/Users/aketz/Documents/Data/SWDPPdeerDB.MDB')
d.vit <- sqlFetch(myconn, "VIT")
close(myconn)


###
### Load VIT data
###


# d.vit=mdb.get('~/Documents/Data/SWDPPdeerDB.MDB',tables= "VIT")
names(d.vit)=tolower(gsub('[[:punct:]]',"",names(d.vit)))
names(d.vit)=gsub(" ", "",names(d.vit), fixed = TRUE)
names(d.vit)[1]="lowtag"



d.vit$datedropped=as.character(d.vit$datedropped)
d.vit$juliandropped=yday(ymd(d.vit$datedropped))
d.vit$datefawnfound=as.character(d.vit$datefawnfound)
d.vit$julianfawnfound=yday(ymd(d.vit$datefawnfound))

n.vit=length(d.vit$lowtag)

#reorder vit data by lowtag/lowtag
d.vit=d.vit[order(d.vit$lowtag),]

#extract the individual ids
individs=d.vit$lowtag


#manipulate the julian day for the doe 6865 with best guess based on movement rates

d.vit$juliandropped[9] = 147

#par(mfrow=c(1,1))
#plot(d$julian[d$lowtag==6865],d$distance[d$lowtag==6865])

###
### Load data GPS location data
###

d = matrix(NA,nr=1,nc=13)
#for loop for reading in data, using vector of id's from the vit dataset
for(i in 1:n.vit){
    d.temp = read.csv(paste("F:/GPS_locations_ind/",d.vit$lowtag[i],".csv",sep=""),sep=",",header=TRUE,stringsAsFactors = FALSE, encoding = "UTF-8")
    d.temp$lowtag = d.vit$lowtag[i]
    names(d.temp)=tolower(names(d.temp))
    d=rbind(d,as.matrix(d.temp))
}
d=d[-1,]
d=data.frame(d,stringsAsFactors = FALSE)

for(j in 1:dim(d)[2]){
    d[,j]=str_trim(d[,j],side="both")
}

class.indx=c(5:7,9:12)
for(j in class.indx){
    d[,j]=as.numeric(d[,j])
}

d$lowtag=as.factor(d$lowtag)

head(d)

###
### double checking for duplicate locations
###

summary(duplicated(d))


#calculating julian day and omitting outside of parturition window
d$julian=yday(mdy_hms(d$date_time_local))

###
### increased fixes parturition window
###

start=yday(mdy("03/20/2017")) # May 6 beginning of parturition period
end=yday(mdy("07/07/2017")) #end of parturition period


#removing weird missing values of 5004 trajectory end....
#tail(d[d$lowtag==5004,],50)


###
### subset entire dataset to parturition window
###

d=d[d$julian>start & d$julian <= end,]


rm5004 = (dim(d[d$lowtag==5004,])[1]-15):dim(d[d$lowtag==5004,])[1]
d=d[-rm5004,]
tail(d[d$lowtag==5004,])

#or keep all data, but exclude first couple hundred observations
#for plotting all data, not just window of parturition

# d=data.frame(d %>% group_by(lowtag) %>% slice(160:(n()-300)))
# n.obs=table(d$lowtag)
# n.obs
# 
# indxing=c(1,rep(NA,n.vit-1),dim(d)[1])
# for(j in 2:n.vit){
#     indxing[j]=indxing[j-1]+n.obs[j-1]
# }
# d[indxing,]
# 
# d[d$lowtag==individs[1],]

### Adding Vit data to main GPS dataframe

d$dropped = 0

records=dim(d)[1]
records
for(i in 1:records){
    for(j in 1:n.vit){
        if(d$lowtag[i]==d.vit$lowtag[j]){
            if(d$julian[i]==d.vit$juliandropped[j])d$dropped[i]=1
        }
    }
}
sum(d$dropped)


###
### Converting date to POSIXct format
###

d$date_time_local=as.POSIXct(strptime(d$date_time_local,format="%m-%d-%Y %H:%M:%S"),tz="CST6CDT")
d$date_time_gmt=as.POSIXct(strptime(d$date_time_gmt,format="%m-%d-%Y %H:%M:%S"),tz="GMT")

dim(d)

#Create time lag between successive locations to censor data if needed.
timediff <- diff(d$date_time_local)
timediff
d=d[-1,]
d$timediff <-round(as.numeric(abs(timediff)))
rm=which(d$timediff>100)
d[rm,]
d=d[-rm,]

names(d)[1]="devname"
head(d)

###
### Dealing with missing data
###

#option 1, just remove
# d=d[!is.na(d$longitude),]

#option 2 impute
d$missing=0
for( i in 1:dim(d)[1]){
    if(is.na(d$longitude[i]))d$missing[i]=1
}

miss.per.ind=c()
for(j in 1:n.vit){
    miss.per.ind=c(miss.per.ind,sum(d$missing[d$lowtag==individs[j]]))
}
miss.per.ind

d=d[-c(1:3),]

for(i in 2:(dim(d)[1]-1)){
    if(is.na(d$longitude[i])){
        a=i-1
        while(is.na(d$longitude[a])){a=a-1}
        b=i+1
        while(is.na(d$longitude[b])){b=b+1}
        d[i,6:5] = midPoint(d[a,6:5],d[b,6:5])
    }
}


###
### Without projection of datum into R, can use geospatial package to calculate distance and bearings
###

bearing.out=bearing(cbind(d$longitude,d$latitude))
d$bearing=bearing.out
dist.out = distHaversine(d[,6:5])
d=d[-1,]
d$distance = dist.out

d=d[-c(dim(d)[1]-1,dim(d)[1]),]#remove last 2 entries which are NA and NaN

tail(d)


###
### Projections! 
###

# setup coordinates
coords = cbind(d$longitude, d$latitude)
sp = SpatialPoints(coords)

# make spatial data frame
# spdf = SpatialPointsDataFrame(coords, d)
spdf = SpatialPointsDataFrame(sp, d)

# EPSG strings
latlong = "+init=epsg:4326"
proj4string(spdf) = CRS(latlong)

d.sp.proj = spTransform(spdf, CRS("+proj=tmerc +lat_0=0 +lon_0=-90 +k=0.9996 +x_0=520000
                                   +y_0=-4480000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))

d=data.frame(d.sp.proj)

#20 and 21 are the coordinates in UTMs x y

###
### Animal paths
### Using the adelehabitatLT

d.traj <- as.ltraj(d[,20:21], date=d$date_time_local,id=d$lowtag)
plot(d.traj)

#Plot the trajectory for each individual - the numbers correspond to the ID in the d.traj object above
#blue is start
#red is end point

par(mfrow=c(2,2))
for(i in 1:n.vit){
    plot(d.traj[i])
}
###
### The last movements of individ1 is crazy
### must remove those points
plot(d.traj[4])




hist(d.traj[1], "dist", freq = TRUE)
hist(d.traj[2], "dist", freq = TRUE)

plotltr(d.traj[1], "dt/3600/24")

#converts traj object to data frame
dfdeer <- ld(d.traj)
dfdeer$id=as.character(dfdeer$id)
dfdeer=rbind(rep(NA,dim(dfdeer)[2]),dfdeer)

dfdeer=dfdeer[-dim(dfdeer)[1],]

d$rel.angle=dfdeer$rel.angle
d$dist.traj=dfdeer$dist
d$R2n=dfdeer$R2n
d$dx=dfdeer$dx
d$dy=dfdeer$dy

d$dt=dfdeer$dt


###
### Parameters returned by adehabitatLT object model. 
###

# • dx, dy, dt: these parameters measured at relocation i describe the increments
# of the x and y directions and time between the relocations i and
# i + 1. Such parameters are often used in the framework of stochastic differential
# equation modelling (e.g. Brillinger et al. 2004, Wiktorsson et al.
#                     2004);
# • dist: the distance between successive relocations is often used in animal
# movement analysis (e.g. Root and Kareiva 1984, Marsh and Jones 1988);

# • abs.angle: the absolute angle αi between the x direction and the step
# built by relocations i and i + 1 is sometimes used together with the parameter
# dist to fit movement models (e.g. Marsh and Jones 1988);

# • rel.angle: the relative angle βi measures the change of direction between
# the step built by relocations i − 1 and i and the step built by relocations
# i and i + 1 (often called “turning angle”). It is often used together with
# the parameter dist to fit movement models (e.g. Root and Kareiva 1984,
#                                            Marsh and Jones 1988);

# • R2n: the squared distance between the first relocation of the trajectory
# and the current relocation is often used to test some movements models
# (e.g. the correlated random walk, see the seminal paper of Kareiva and
# Shigesada, 1983).

###
### remove the initial index of each individual 
### because it makes no sense to calculate distances 
### between the locations of the individuals
###

remove.indx=which(is.na(d$dist.traj))
d[remove.indx,]
d=d[-remove.indx,]

###
###EDA plots
###

plot(d$dist.traj,d$distance)#identical
plot(d$bearing,d$rel.angle)#whoa mess
plot(d$altitude,d$dist.traj)
plot(d$distance,col = as.numeric(d$lowtag))
plot(d$altitude,col = as.numeric(d$lowtag))
plot(d$rel.angle,col = as.numeric(d$lowtag))
plot(d$bearing,col = as.numeric(d$lowtag))
plot(d$dx,col = as.numeric(d$lowtag))
plot(d$dy,col = as.numeric(d$lowtag))
plot(d$dt,col = as.numeric(d$lowtag))
plot(d$R2n[d$lowtag!=5004],col = as.numeric(d$lowtag[d$lowtag!=5004]))#exclude outlier 
plot(d$temp,col = as.numeric(d$lowtag))
plot(d$temp,d$distance,col = as.numeric(d$lowtag))
plot(d$distance,d$temp,col = as.numeric(d$lowtag))
plot(d$distance[d$lowtag!=5004],d$dropped[d$lowtag!=5004],col=as.numeric(d$lowtag[d$lowtag!=5004]))

###
### Calculating and plotting moving window metrics
###

n.obs=as.numeric(table(d$lowtag))
n.obs

#bookkeeping
vit.drop=matrix(NA,ncol=n.vit,nrow=max(n.obs))
julian.out=matrix(NA,ncol=n.vit,nrow=max(n.obs))

#min/max features
min.distance=matrix(NA,ncol=n.vit,nrow=max(n.obs))
max.angle=matrix(NA,ncol=n.vit,nrow=max(n.obs))

#mean features
mu.distance=matrix(NA,ncol=n.vit,nrow=max(n.obs))
mu.turn.angle=matrix(NA,ncol=n.vit,nrow=max(n.obs))
mu.bearing=matrix(NA,ncol=n.vit,nrow=max(n.obs))
mu.r2n=matrix(NA,ncol=n.vit,nrow=max(n.obs))
mu.altitude=matrix(NA,ncol=n.vit,nrow=max(n.obs))
mu.temp=matrix(NA,ncol=n.vit,nrow=max(n.obs))
mu.dx=matrix(NA,ncol=n.vit,nrow=max(n.obs))
mu.dy=matrix(NA,ncol=n.vit,nrow=max(n.obs))

#variation features
sig.distance=matrix(NA,ncol=n.vit,nrow=max(n.obs))
sig.turn.angle=matrix(NA,ncol=n.vit,nrow=max(n.obs))
sig.bearing=matrix(NA,ncol=n.vit,nrow=max(n.obs))
sig.r2n=matrix(NA,ncol=n.vit,nrow=max(n.obs))
sig.altitude=matrix(NA,ncol=n.vit,nrow=max(n.obs))
sig.temp=matrix(NA,ncol=n.vit,nrow=max(n.obs))
sig.dx=matrix(NA,ncol=n.vit,nrow=max(n.obs))
sig.dy=matrix(NA,ncol=n.vit,nrow=max(n.obs))

#correlation features
corr.dist.angle=matrix(NA,ncol=n.vit,nrow=max(n.obs))
corr.dist.bear=matrix(NA,ncol=n.vit,nrow=max(n.obs))
corr.dist.alt=matrix(NA,ncol=n.vit,nrow=max(n.obs))
corr.angle.bear=matrix(NA,ncol=n.vit,nrow=max(n.obs))
corr.angle.alt=matrix(NA,ncol=n.vit,nrow=max(n.obs))
corr.bear.alt=matrix(NA,ncol=n.vit,nrow=max(n.obs))


w=6
for(j in 1:n.vit){
    d.tmp=d[d$lowtag==individs[j],]
    for(i in (1+w):dim(d.tmp)[1]){
        #relative angle
        vit.drop[i,j]=d.tmp$dropped[i]
        julian.out[i,j]=d.tmp$julian[i]
        
        #min/max features
        min.distance[i,j]=min(d.tmp$distance[(i-w):i])
        max.angle[i,j]=max(d.tmp$rel.angle[(i-w):i])
        
        #mean features
        mu.distance[i,j]=mean(d.tmp$distance[(i-w):i])
        mu.turn.angle[i,j]=mean(d.tmp$rel.angle[(i-w):i])
        mu.bearing[i,j]=mean(d.tmp$bearing[(i-w):i])
        mu.r2n[i,j]=mean(d.tmp$R2n[(i-w):i])
        mu.altitude[i,j]=mean(d.tmp$altitude[(i-w):i])
        mu.temp[i,j]=mean(d.tmp$temp[(i-w):i])
        mu.dx[i,j]=mean(d.tmp$dx[(i-w):i])
        mu.dy[i,j]=mean(d.tmp$dy[(i-w):i])
        
        #variation features
        sig.distance[i,j]=sd(d.tmp$distance[(i-w):i])
        sig.turn.angle[i,j]=sd(d.tmp$rel.angle[(i-w):i])
        sig.bearing[i,j]=sd(d.tmp$bearing[(i-w):i])
        sig.r2n[i,j]=sd(d.tmp$R2n[(i-w):i])
        sig.altitude[i,j]=sd(d.tmp$altitude[(i-w):i])
        sig.temp[i,j]=sd(d.tmp$temp[(i-w):i])
        sig.dx[i,j]=sd(d.tmp$dx[(i-w):i])
        sig.dy[i,j]=sd(d.tmp$dy[(i-w):i])
        
        #correlation features
        corr.dist.angle[i,j]=corr(cbind(d.tmp$distance[(i-w):i],d.tmp$rel.angle[(i-w):i]))
        corr.dist.bear[i,j]=corr(cbind(d.tmp$distance[(i-w):i],d.tmp$bearing[(i-w):i]))
        corr.dist.alt[i,j]=corr(cbind(d.tmp$distance[(i-w):i],d.tmp$altitude[(i-w):i]))
        corr.angle.bear[i,j]=corr(cbind(d.tmp$rel.angle[(i-w):i],d.tmp$bearing[(i-w):i]))
        corr.angle.alt[i,j]=corr(cbind(d.tmp$rel.angle[(i-w):i],d.tmp$altitude[(i-w):i]))
        corr.bear.alt[i,j]=corr(cbind(d.tmp$bearing[(i-w):i],d.tmp$altitude[(i-w):i]))
        
                
        }
}
# angle.dist.corr=angle.dist.corr[,(w+2):dim(angle.dist.corr)[1]]
# vit.drop=vit.drop[,(w+2):dim(vit.drop)[2]]
# julian.out=julian.out[,(w+2):dim(julian.out)[2]]


vit.drop=vit.drop[(w+2):dim(vit.drop)[1],]
julian.out=julian.out[(w+2):dim(julian.out)[1],]

#min/max features
min.distance=min.distance[(w+2):dim(min.distance)[1],]
max.angle=max.angle[(w+2):dim(max.angle)[1],]

#mean features
mu.distance=mu.distance[(w+2):dim(mu.distance)[1],]
mu.turn.angle=mu.turn.angle[(w+2):dim(mu.turn.angle)[1],]
mu.bearing=mu.bearing[(w+2):dim(mu.bearing)[1],]
mu.r2n=mu.r2n[(w+2):dim(mu.r2n)[1],]
mu.altitude=mu.altitude[(w+2):dim(mu.altitude)[1],]
mu.temp=mu.temp[(w+2):dim(mu.temp)[1],]
mu.dx=mu.dx[(w+2):dim(mu.dx)[1],]
mu.dy=mu.dy[(w+2):dim(mu.dy)[1],]

mu.distance=apply(mu.distance,2,scale)
mu.turn.angle=apply(mu.turn.angle,2,scale)

head(mu.distance)
n.obs

#plotting density of mean distances
par(mfrow=c(3,2))
for(j in 1:n.vit){
    plot(density(mu.distance[,j],na.rm=TRUE))
    plot(mu.distance[,j],na.rm=TRUE)
}

par(mfrow=c(3,2))
for(j in 1:n.vit){
    plot(density(mu.turn.angle[,j],na.rm=TRUE))
    plot(mu.turn.angle[,j],na.rm=TRUE)
}


#variation features
sig.distance=sig.distance[(w+2):dim(sig.distance)[1],]
sig.turn.angle=sig.turn.angle[(w+2):dim(sig.turn.angle)[1],]
sig.bearing=sig.bearing[(w+2):dim(sig.bearing)[1],]
sig.r2n=sig.r2n[(w+2):dim(sig.r2n)[1],]
sig.altitude=sig.altitude[(w+2):dim(sig.altitude)[1],]
sig.temp=sig.temp[(w+2):dim(sig.temp)[1],]
sig.dx=sig.dx[(w+2):dim(sig.dx)[1],]
sig.dy=sig.dy[(w+2):dim(sig.dy)[1],]

#correlation features
corr.dist.angle=corr.dist.angle[(w+2):dim(corr.dist.angle)[1],]
corr.dist.bear=corr.dist.bear[(w+2):dim(corr.dist.bear)[1],]
corr.dist.alt=corr.dist.alt[(w+2):dim(corr.dist.alt)[1],]
corr.angle.bear=corr.angle.bear[(w+2):dim(corr.angle.bear)[1],]
corr.angle.alt=corr.angle.alt[(w+2):dim(corr.angle.alt)[1],]
corr.bear.alt=corr.bear.alt[(w+2):dim(corr.bear.alt)[1],]


vit.records = rep( list(list()), n.vit )
for(j in 1:n.vit){
    vit.records[[j]]=which(vit.drop[,j]==1)
}
vit.records


par(mfrow=c(3,2))
for(j in 1:n.vit){
    plot(julian.out[,j],corr.dist.angle[,j],main=individs[j])
}


###
### Plots min/max features over window
###

par(mfrow=c(2,1))
for(j in 1:n.vit){
    plot(julian.out[,j],min.distance[,j],main=d.vit$lowtag[j])
    abline(v=d.vit$juliandropped[j],col="darkred",lwd=2)
    plot(julian.out[,j],max.angle[,j])
    abline(v=d.vit$juliandropped[j],col="darkred",lwd=2)
}

###
### Plots mean features
###

par(mfrow=c(4,2))
for(j in 1:n.vit){
    plot(julian.out[,j],mu.distance[,j],main=d.vit$lowtag[j])
    abline(v=d.vit$juliandropped[j],col="darkred",lwd=2)
    plot(julian.out[,j],mu.turn.angle[,j])
    abline(v=d.vit$juliandropped[j],col="darkred",lwd=2)
    plot(julian.out[,j],mu.bearing[,j])
    abline(v=d.vit$juliandropped[j],col="darkred",lwd=2)
    plot(julian.out[,j],mu.r2n[,j])
    abline(v=d.vit$juliandropped[j],col="darkred",lwd=2)
    plot(julian.out[,j],mu.altitude[,j])
    abline(v=d.vit$juliandropped[j],col="darkred",lwd=2)
    plot(julian.out[,j],mu.temp[,j])
    abline(v=d.vit$juliandropped[j],col="darkred",lwd=2)
    plot(julian.out[,j],mu.dx[,j])
    abline(v=d.vit$juliandropped[j],col="darkred",lwd=2)
    plot(julian.out[,j],mu.dy[,j])
    abline(v=d.vit$juliandropped[j],col="darkred",lwd=2)
}


###
### Plots Variation Features
###
par(mfrow=c(4,2))
for(j in 1:n.vit){
    plot(julian.out[,j],sig.distance[,j],main=d.vit$lowtag[j])
    abline(v=d.vit$juliandropped[j],col="darkred",lwd=2)
    plot(julian.out[,j],sig.turn.angle[,j])
    abline(v=d.vit$juliandropped[j],col="darkred",lwd=2)
    plot(julian.out[,j],sig.bearing[,j])
    abline(v=d.vit$juliandropped[j],col="darkred",lwd=2)
    plot(julian.out[,j],sig.r2n[,j])
    abline(v=d.vit$juliandropped[j],col="darkred",lwd=2)
    plot(julian.out[,j],sig.altitude[,j])
    abline(v=d.vit$juliandropped[j],col="darkred",lwd=2)
    plot(julian.out[,j],sig.temp[,j])
    abline(v=d.vit$juliandropped[j],col="darkred",lwd=2)
    plot(julian.out[,j],sig.dx[,j])
    abline(v=d.vit$juliandropped[j],col="darkred",lwd=2)
    plot(julian.out[,j],sig.dy[,j])
    abline(v=d.vit$juliandropped[j],col="darkred",lwd=2)
}


par(mfrow=c(3,2))
for(j in 1:n.vit){
    plot(julian.out[,j],corr.dist.angle[,j],main=d.vit$lowtag[j])
    abline(v=d.vit$juliandropped[j],col="darkred",lwd=2)
    plot(julian.out[,j],corr.dist.bear[,j])
    abline(v=d.vit$juliandropped[j],col="darkred",lwd=2)
    plot(julian.out[,j],corr.dist.alt[,j])
    abline(v=d.vit$juliandropped[j],col="darkred",lwd=2)
    plot(julian.out[,j],corr.angle.bear[,j])
    abline(v=d.vit$juliandropped[j],col="darkred",lwd=2)
    plot(julian.out[,j],corr.angle.alt[,j])
    abline(v=d.vit$juliandropped[j],col="darkred",lwd=2)
    plot(julian.out[,j],corr.bear.alt[,j])
    abline(v=d.vit$juliandropped[j],col="darkred",lwd=2)
}


par(mfrow=c(3,2))
for(j in 1:n.vit){
    plot(julian.out[,j],corr.dist.angle[,j],main=d.vit$lowtag[j])
    abline(v=d.vit$juliandropped[j],col="darkred",lwd=2)
    plot(julian.out[,j],corr.dist.bear[,j])
    abline(v=d.vit$juliandropped[j],col="darkred",lwd=2)
    plot(julian.out[,j],corr.dist.alt[,j])
    abline(v=d.vit$juliandropped[j],col="darkred",lwd=2)
    plot(julian.out[,j],corr.angle.bear[,j])
    abline(v=d.vit$juliandropped[j],col="darkred",lwd=2)
    plot(julian.out[,j],corr.angle.alt[,j])
    abline(v=d.vit$juliandropped[j],col="darkred",lwd=2)
    plot(julian.out[,j],corr.bear.alt[,j])
    abline(v=d.vit$juliandropped[j],col="darkred",lwd=2)
}


#####################################################################################################
###
### Anomaly Detection analysis
###
#####################################################################################################

d$step = d$distance/d$timediff
d$logstep = log(d$step)

###
### center and scale variables for prediction
###

rescale=c(7,10,18,23,25,26,29)
names(d)[rescale]
# [1] "altitude"  "temp"      "bearing"   "rel.angle" "R2n"       "dx"        "step"     

scalealt=c()
scaletemp=c()
scalebear=c()
scaleangle=c()
scaleR2n=c()
scaledx=c()
scalestep=c()
logstepscale=c()

#center and scale variables for predicition
for(j in 1:n.vit){
    scalealt = c(scalealt,scale(d$altitude[d$lowtag==individs[j]]))
    scaletemp = c(scaletemp,scale(d$temp[d$lowtag==individs[j]]))
    scalebear = c(scalebear,scale(d$bearing[d$lowtag==individs[j]]))
    scaleangle = c(scaleangle,scale(d$rel.angle[d$lowtag==individs[j]]))
    scaleR2n = c(scaleR2n,scale(d$R2n[d$lowtag==individs[j]]))
    scaledx = c(scaledx,scale(d$dx[d$lowtag==individs[j]]))
    scalestep = c(scalestep,scale(d$step[d$lowtag==individs[j]]))
    logstepscale = c(logstepscale,scale(d$logstep[d$lowtag==individs[j]]))
    }

d$scalealt = scalealt
d$scaletemp = scaletemp
d$scaleangle = scaleangle
d$scalebear = scalebear
d$scaledx = scaledx
d$scalestep = scalestep
d$scaleR2n=scaleR2n
d$logscalestep=logstepscale



names(d)[30:37]
par(mfrow=c(2,2))
for(j in 30:37){
    hist(d[,j],breaks=50,main=names(d)[j])
}

#normal enough: altitude, temp, dx,r2n,logscalestep 
#bad: scaleangle, bearing



hist(d$scalestep)

d$scalestep = scale(d$step)
hist(d$scalestep)
head(d)

par(mfrow=c(3,1))
for(j in 1:n.vit){
    
    hist(log(d$step[d$lowtag==individs[j]]),breaks=50,main=individs[j])
    plot(density(log(d$step[d$lowtag==individs[j]])),main=individs[j])
    plot(d$julian[d$lowtag==individs[j]],log(d$step[d$lowtag==individs[j]]))
    abline(v=d.vit$juliandropped[j],col=2)
}


##############################################################################################
###
### Univariate Anomaly Detection with distance only
###
##############################################################################################
#parturition window = may 6 through jul 7 
part.window=yday("2017-05-06")


ind.mean=rep(NA,n.vit)
ind.sd=rep(NA,n.vit)
for(j in 1:n.vit){
    d.temp=d[d$lowtag==individs[j],]
    plot(density(d.temp$logstep[d.temp$julian<part.window]))
    ind.mean[j]=mean(d.temp$logstep[d.temp$julian<part.window])
    ind.sd[j]=sd(d.temp$logstep[d.temp$julian<part.window])
}

n.out=c()
for(j in 1:n.vit){
    d.temp=d[d$lowtag==individs[j] & d$julian>part.window,]
    n.temp=dim(d.temp)[1]
    n.out=c(n.out,n.temp)
}
n.out.max=max(n.out)

detect.density=matrix(NA,nrow=max(n.out),ncol=n.vit)
detect.quant=rep(NA,n.vit)
threshold.density=rep(NA,n.vit)
threshold.step=matrix(NA,nrow=10,ncol=n.vit)
alarm=rep(list(),n.vit)

par(mfrow=c(2,1))
for(j in 1:n.vit){
    d.temp=d[d$lowtag==individs[j] & d$julian>part.window,]
    n.temp=dim(d.temp)[1]
    # acf(d.temp$logstep)
    detect.quant[j]=quantile(d.temp$logstep,.07)
    for(i in 1:n.temp){
        threshold.density[j] = dnorm(detect.quant[j],ind.mean[j],ind.sd[j])
        detect.density[i,j] = dnorm(d.temp$logstep[i],ind.mean[j],ind.sd[j])
    }
    alarm[[j]]= d.temp[detect.density[,j]<threshold.density[j],]
    cat(sum(alarm[[j]]$dropped,na.rm=TRUE),'\t',j,'\n')
}

for(j in 1:n.vit){
    plot(alarm[[j]]$julian,alarm[[j]]$logstep)
    abline(v=d.vit$juliandropped[j])
}
par(mfrow=c(2,1))
for(j in 1:n.vit){
    plot(d$julian[d$lowtag==individs[j]],d$logstep[d$lowtag==individs[j]])
    points(alarm[[j]]$julian,alarm[[j]]$logstep,col=2)
    abline(v=d.vit$juliandropped[j])
    
}



##############################################################################################
###
### Evaluate Anomaly Detection within individual using ground truthing
### alarm = list of detected data points
### possible.hits = vector of the true number of vit drop ground truths
###
##############################################################################################

###
### number of true possible positives for vit drop
###

possible.hits=rep(NA,n.vit)
for(j in 1:n.vit){
    possible.hits[j] = sum(d$dropped[d$lowtag==individs[j]])
}
possible.hits

###
### Function for calculating metrics
### 

evaluate = function(alarm,possible.hits,n.vit=12){
    
    tp = rep(NA,n.vit)
    fp = rep(NA,n.vit)
    fn = rep(NA,n.vit)
    out.prec = rep(NA,n.vit)
    out.recall = rep(NA,n.vit)
    out.F1 = rep(NA,n.vit)
    for(j in 1:n.vit){
        # TRUE POSITIVE:  groundtruth data  says it's an anomaly and so algorithm does.
        tp[j] = sum(alarm[[j]]$dropped,na.rm=TRUE)
        
        # FALSE POSITIVE:  groundtruth data says it's not an anomaly, but algorithm says anomaly.
        fp[j] = sum(alarm[[j]]$dropped==0,na.rm=TRUE)
        
        # FALSE NEGATIVE: groundtruth data says it's an anomaly, but algorithm says not anomaly.
        fn[j] = possible.hits[j] - tp[j]
        
        # precision and recall
        out.prec[j] = tp[j]/(tp[j]+fp[j])
        out.recall[j] = tp[j]/(tp[j]+fn[j])
        
        # F1 value
        out.F1[j] = (2*out.prec[j]*out.recall[j])/(out.prec[j]+out.recall[j])
    }
    return(list(tp=tp,fp=fp,fn=fn,out.prec=out.prec,out.recall=out.recall,out.F1=out.F1))
    
}


out.eval=evaluate(alarm=alarm,possible.hits=possible.hits,n.vit=n.vit)
out.eval$tp

##############################################################################################
###
### Univariate Anomaly Detection Function
###
##############################################################################################

anomalyDetectUniv = function(n.out.max,n.vit,part.window,id,d,eps,im,is){
    detect.quant=rep(NA,n.vit)
    threshold.density=rep(NA,n.vit)
    detect.density=matrix(NA,nr=n.out.max,nc=n.vit)
    alarm=rep(list(),n.vit)
    epsilon=eps
    for(j in 1:n.vit){
        d.temp=d[d$lowtag==id[j] & d$julian>part.window,]
        n.temp=dim(d.temp)[1]
        detect.quant[j]=quantile(d.temp$logstep,epsilon)
        threshold.density[j] = dnorm(detect.quant[j],im[j],is[j])
        
        for(i in 1:n.temp){
            detect.density[i,j] = dnorm(d.temp$logstep[i],im[j],is[j])
        }
        alarm[[j]]= d.temp[detect.density[,j]<threshold.density[j],]
    }
    
    return(list(alarm=alarm,detect.density=detect.density,threshold.density=threshold.density,detect.quant=detect.quant))
}


##############################################################################################
###
### Multivariate Anomaly Detection within a function, 
###
##############################################################################################


anomalyDetect = function(n.vit,part.window=126,id,d,eps=.1,covs.indx){
    
    covs.indx =covs.indx
    numcovs = length(covs.indx)
    
    ind.mean=matrix(NA,nr=numcovs,nc=n.vit)
    ind.sd=matrix(NA,nr=numcovs,nc=n.vit)
    for(j in 1:n.vit){
        d.temp=d[d$lowtag==individs[j],]
        for(k in 1:numcovs){
            ind.mean[k,j]=mean(d.temp[d.temp$julian<part.window,covs.indx[k]],na.rm=TRUE)
            ind.sd[k,j]=sd(d.temp[d.temp$julian<part.window,covs.indx[k]],na.rm=TRUE)
        }
    }
    
    n.out=c()
    for(j in 1:n.vit){
        d.temp=d[d$lowtag==individs[j] & d$julian>part.window,]
        n.temp=dim(d.temp)[1]
        n.out=c(n.out,n.temp)
    }
    n.out.max=max(n.out)
    
    if(numcovs==1){
        ind.mean=c(ind.mean)
        ind.sd=c(ind.sd)
        epsilon=eps
        anomalyDetectUniv(n.out.max,n.vit=12,part.window=126,id=individs,d,eps=epsilon,im=ind.mean,is=ind.sd)
    }
    
    else{
        detect.quant=matrix(NA,nr=numcovs,nc=n.vit)
        epsilon=rep(eps,n.vit)
        threshold.density=rep(NA,n.vit)
        detect.density=matrix(NA,nr=n.out.max,nc=n.vit)
        alarm=rep(list(),n.vit)
        
        for(j in 1:n.vit){
            d.temp=d[d$lowtag==individs[j] & d$julian>part.window,]
            n.temp=dim(d.temp)[1]
            for(k in 1:numcovs){
                detect.quant[k,j]=quantile(d.temp[,covs.indx[k]],epsilon[j],na.rm=TRUE)
            }
            threshold.density[j] = dmvnorm(detect.quant[,j],ind.mean[,j],diag(ind.sd[,j]))
            for(i in 1:n.temp){
                detect.density[i,j] = dmvnorm(d.temp[i,covs.indx],ind.mean[,j],diag(ind.sd[,j]))
            }
            alarm[[j]]= d.temp[detect.density[,j]<threshold.density[j],]
        }        
    }#end else
    
    return(list(alarm=alarm,detect.density=detect.density,threshold.density=threshold.density,detect.quant=detect.quant))
}


###
### Function for calculating anomaly detection evaluation metrics
### Individual based
### 

evaluate.time.ind = function(alarm,possible.hits,n.vit=12,vitdropday){
    
    for(j in 1:n.vit){
        alarm[[j]]=alarm[[j]][alarm[[j]]$julian<=vitdropday[j],]
    }
    
    tp = rep(NA,n.vit)
    fp = rep(NA,n.vit)
    fn = rep(NA,n.vit)
    out.prec = rep(NA,n.vit)
    out.recall = rep(NA,n.vit)
    out.F1 = rep(NA,n.vit)
    for(j in 1:n.vit){
        # TRUE POSITIVE:  groundtruth data  says it's an anomaly and so algorithm does.
        tp[j] = sum(alarm[[j]]$dropped,na.rm=TRUE)
        
        # FALSE POSITIVE:  groundtruth data says it's not an anomaly, but algorithm says anomaly.
        fp[j] = sum(alarm[[j]]$dropped==0,na.rm=TRUE)
        
        # FALSE NEGATIVE: groundtruth data says it's an anomaly, but algorithm says not anomaly.
        fn[j] = possible.hits[j] - tp[j]
        
        # precision and recall
        out.prec[j] = tp[j]/(tp[j]+fp[j])
        out.recall[j] = tp[j]/(tp[j]+fn[j])
        
        # F1 value
        out.F1[j] = (2*out.prec[j]*out.recall[j])/(out.prec[j]+out.recall[j])
    }
    return(list(tp=tp,fp=fp,fn=fn,out.prec=out.prec,out.recall=out.recall,out.F1=out.F1))
    
}



starting=Sys.time()
testnum=15
eval=rep(list(),testnum)
for(k in 1:testnum){
    test = anomalyDetect(n.vit=12,part.window=part.window,id=individs,d,eps=k/100,covs.indx=c(7,10,23,26,30))
    eval[[k]]=evaluate.time.ind(alarm=test$alarm,possible.hits=possible.hits,n.vit=12,vitdropday = d.vit$juliandropped)
}

ending=Sys.time()
ending-starting

###
### evaluation metrics and plots
###

#epsilon values it tried
epsilon = seq(.01,.15,by=.01)


precision=matrix(NA,nc=n.vit)
for(k in 1:testnum){
    precision=rbind(precision,eval[[k]]$out.prec)
}
precision=precision[-1,]
par(mfrow=c(3,2))
for(j in 1:n.vit){
    plot(epsilon,precision[,j],main=d.vit$lowtag[j],ylim=c(0,1))
}


RECALL=matrix(NA,nc=n.vit)
for(k in 1:testnum){
    RECALL=rbind(RECALL,eval[[k]]$out.recall)
}
RECALL=RECALL[-1,]
par(mfrow=c(3,2))
for(j in 1:n.vit){
    plot(epsilon,RECALL[,j],main=d.vit$lowtag[j],ylim=c(0,1))
}



F1score=matrix(NA,nc=n.vit)
for(k in 1:testnum){
    F1score=rbind(F1score,eval[[k]]$out.F1)
}
F1score=F1score[-1,]
par(mfrow=c(3,2))
for(j in 1:n.vit){
    plot(epsilon,F1score[,j],main=d.vit$lowtag[j],ylim=c(0,1))
}




TP=matrix(NA,nc=n.vit)
for(k in 1:testnum){
    TP=rbind(TP,eval[[k]]$tp)
}
TP=TP[-1,]
par(mfrow=c(3,2))
for(j in 1:n.vit){
    plot(epsilon,TP[,j],main=d.vit$lowtag[j],ylim=c(0,10))
}

FP=matrix(NA,nc=n.vit)
for(k in 1:testnum){
    FP=rbind(FP,eval[[k]]$fp)
}
FP=FP[-1,]
par(mfrow=c(3,2))
for(j in 1:n.vit){
    plot(epsilon,FP[,j],main=d.vit$lowtag[j])
}

FN=matrix(NA,nc=n.vit)
for(k in 1:testnum){
    FN=rbind(FN,eval[[k]]$fn)
}
FN=FN[-1,]
par(mfrow=c(3,2))
for(j in 1:n.vit){
    plot(epsilon,FN[,j],main=d.vit$lowtag[j],ylim=c(0,12))
}

par(mfrow=c(3,2))
FP=cbind(FP,epsilon)
for(j in 1:n.vit){
    plot(FP[,j],TP[,j],main=d.vit$lowtag[j],type="n")
    text(FP[,j],TP[,j],labels = FP[,n.vit+1])
}

par(mfrow=c(3,2))
for(j in 1:n.vit){
    plot(FP[,j],FN[,j],main=d.vit$lowtag[j],type="n")
    text(FP[,j],FN[,j],labels = FP[,n.vit+1])
}

#########################################################################################
###
### Function for calculating anomaly detection evaluation metrics
### Across entire population
### 
#########################################################################################


evaluate.time.pop = function(alarm,possible.hits,n.vit=12,vitdropday){
    
    for(j in 1:n.vit){
        alarm[[j]]=alarm[[j]][alarm[[j]]$julian<=vitdropday[j],]
    }
    
    tp = rep(NA,n.vit)
    fp = rep(NA,n.vit)
    fn = rep(NA,n.vit)

    for(j in 1:n.vit){
        
        # TRUE POSITIVE:  groundtruth data  says it's an anomaly and so algorithm does.
        tp[j] = sum(alarm[[j]]$dropped,na.rm=TRUE)
        
        # FALSE POSITIVE:  groundtruth data says it's not an anomaly, but algorithm says anomaly.
        fp[j] = sum(alarm[[j]]$dropped==0,na.rm=TRUE)
        
        # FALSE NEGATIVE: groundtruth data says it's an anomaly, but algorithm says not anomaly.
        fn[j] = possible.hits[j] - tp[j]
        
    }
    
    tp.pop = sum(tp)
    fp.pop = sum(fp)
    fn.pop = sum(fn)

    # precision and recall
    out.prec = tp.pop/(tp.pop+fp.pop)
    out.recall = tp.pop/(tp.pop+fn.pop)
    
    # F1 value
    out.F1 = (2*out.prec*out.recall)/(out.prec+out.recall)
    
        
    return(list(tp=tp,fp=fp,fn=fn,out.prec=out.prec,out.recall=out.recall,out.F1=out.F1,tp.pop=tp.pop,fp.pop=fp.pop,fn.pop=fn.pop))
    
}


###
### Run the anomaly detection across different values of epsilon
###




names(d)
starting=Sys.time()
testnum=15
eval=rep(list(),testnum)
for(k in 1:testnum){
    test = anomalyDetect(n.vit=12,part.window=part.window,id=individs,d,eps=k/100,covs.indx=31:37)
    eval[[k]]=evaluate.time.pop(alarm=test$alarm,possible.hits=possible.hits,n.vit=12,vitdropday = d.vit$juliandropped)
}

ending=Sys.time()
ending-starting




###Output evaluation metrics for optimizing epsilon


#epsilon values it tried
epsilon = seq(.01,.15,by=.01)

par(mfrow=c(3,2))
precision=rep(NA,n.vit)
for(k in 1:testnum){
    precision[k]=eval[[k]]$out.prec
}
plot(epsilon,precision,main="Precision")


RECALL.pop=rep(NA,n.vit)
for(k in 1:testnum){
    RECALL.pop[k]=eval[[k]]$out.recall
}
plot(epsilon,RECALL.pop,main="RECALL.pop/Sensitivity")

F1.pop=rep(NA,n.vit)
for(k in 1:testnum){
    F1.pop[k]=eval[[k]]$out.F1
}
plot(epsilon,F1.pop,main="F1 Score")


TP.pop=rep(NA,n.vit)
for(k in 1:testnum){
    TP.pop[k]=eval[[k]]$tp.pop
}
plot(epsilon,TP.pop,main="TP.pop")

FN.pop=rep(NA,n.vit)
for(k in 1:testnum){
    FN.pop[k]=eval[[k]]$fn.pop
}
plot(epsilon,FN.pop,main="FN.pop (120 total possible)")



FP.pop=rep(NA,n.vit)
for(k in 1:testnum){
    FP.pop[k]=eval[[k]]$fp.pop
}
plot(epsilon,FP.pop,main="FP.pop")


par(mfrow=c(1,1))
plot(TP.pop,FP.pop,main="TP.pop vs FP.pop")


##############################################################################################
###
### Multivariate Anomaly Detection on 1 individual Using Episolon = .12
### Changing predictors
###
##############################################################################################


anomalyDetectUniv = function(n.out.max,n.vit,part.window,id,d,eps,im,is){
    if(n.vit==1){
        epsilon=eps
        d.temp=d[d$julian>part.window,]
        n.temp=dim(d.temp)[1]
        detect.quant=quantile(d.temp$logstep,epsilon)
        threshold.density = dnorm(detect.quant,im,is)
        for(i in 1:n.temp){
            detect.density[i] = dnorm(d.temp$logstep[i],im,is)
        }
        alarm= list(d.temp[detect.density<threshold.density,])
    }else{
        detect.quant=rep(NA,n.vit)
        threshold.density=rep(NA,n.vit)
        detect.density=matrix(NA,nr=n.out.max,nc=n.vit)
        alarm=rep(list(),n.vit)
        epsilon=eps
        for(j in 1:n.vit){
            d.temp=d[d$lowtag==id[j] & d$julian>part.window,]
            n.temp=dim(d.temp)[1]
            detect.quant[j]=quantile(d.temp$logstep,epsilon)
            threshold.density[j] = dnorm(detect.quant[j],im[j],is[j])
            
            for(i in 1:n.temp){
                detect.density[i,j] = dnorm(d.temp$logstep[i],im[j],is[j])
            }
            alarm[[j]]= d.temp[detect.density[,j]<threshold.density[j],]
        }
    }
    return(list(alarm=alarm,detect.density=detect.density,threshold.density=threshold.density,detect.quant=detect.quant))
}


anomalyDetect = function(n.vit,part.window=126,id,d,eps,covs.indx){
    
    if(n.vit==1){
        covs.indx =covs.indx
        numcovs = length(covs.indx)
        
        if(length(eps)==1)eps=rep(eps,numcovs)
        
        ind.mean=rep(NA,numcovs)
        ind.sd=rep(NA,numcovs)
        d.temp=d
        for(k in 1:numcovs){
            ind.mean[k]=mean(d.temp[d.temp$julian<part.window,covs.indx[k]],na.rm=TRUE)
            ind.sd[k]=sd(d.temp[d.temp$julian<part.window,covs.indx[k]],na.rm=TRUE)
        }
        
        d.temp=d[d$julian>part.window,]
        n.out=dim(d.temp)[1]
        n.out.max=n.out
        
        if(numcovs==1){
            ind.mean=c(ind.mean)
            ind.sd=c(ind.sd)
            epsilon=eps
            anomalyDetectUniv(n.out,n.vit=1,part.window=part.window,id=id,d,eps=epsilon,im=ind.mean,is=ind.sd)
        }
        else{
            detect.quant=rep(NA,numcovs)
            if(length(eps)==1)eps=rep(eps,numcovs)
            epsilon=eps
            detect.density=rep(NA,n.out)
            alarm=rep(list(),n.vit)
            d.temp=d[d$julian>part.window,]
            n.temp=dim(d.temp)[1]
            for(k in 1:numcovs){
                detect.quant[k]=quantile(d.temp[,covs.indx[k]],epsilon[k],na.rm=TRUE)
            }
            threshold.density = dmvnorm(detect.quant,ind.mean,diag(ind.sd))
            for(i in 1:n.temp){
                detect.density[i] = dmvnorm(d.temp[i,covs.indx],ind.mean,diag(ind.sd))
            }
            alarm= list(d.temp[detect.density<threshold.density,])
        }#end else numcovs
    }#end if n.vit=1
    else {
        covs.indx =covs.indx
        numcovs = length(covs.indx)
        if(length(eps)==1)epsilon=rep(eps,numcovs)
        ind.mean=matrix(NA,nr=numcovs,nc=n.vit)
        ind.sd=matrix(NA,nr=numcovs,nc=n.vit)
        for(j in 1:n.vit){
            d.temp=d[d$lowtag==id[j],]
            for(k in 1:numcovs){
                ind.mean[k,j]=mean(d.temp[d.temp$julian<part.window,covs.indx[k]],na.rm=TRUE)
                ind.sd[k,j]=sd(d.temp[d.temp$julian<part.window,covs.indx[k]],na.rm=TRUE)
            }
        }
        
        n.out=c()
        for(j in 1:n.vit){
            d.temp=d[d$lowtag==individs[j] & d$julian>part.window,]
            n.temp=dim(d.temp)[1]
            n.out=c(n.out,n.temp)
        }
        n.out.max=max(n.out)
        
        if(numcovs==1){
            ind.mean=c(ind.mean)
            ind.sd=c(ind.sd)
            epsilon=eps
            anomalyDetectUniv(n.out.max,n.vit=12,part.window=126,id=id,d,eps=epsilon,im=ind.mean,is=ind.sd)
        }
        
        else{
            detect.quant=matrix(NA,nr=numcovs,nc=n.vit)
            epsilon=rep(eps,numcovs)
            threshold.density=rep(NA,n.vit)
            detect.density=matrix(NA,nr=n.out.max,nc=n.vit)
            alarm=rep(list(),n.vit)
            
            for(j in 1:n.vit){
                d.temp=d[d$lowtag==individs[j] & d$julian>part.window,]
                n.temp=dim(d.temp)[1]
                for(k in 1:numcovs){
                    detect.quant[k,j]=quantile(d.temp[,covs.indx[k]],epsilon[k],na.rm=TRUE)
                }
                threshold.density[j] = dmvnorm(detect.quant[,j],ind.mean[,j],diag(ind.sd[,j]))
                for(i in 1:n.temp){
                    detect.density[i,j] = dmvnorm(d.temp[i,covs.indx],ind.mean[,j],diag(ind.sd[,j]))
                }
                alarm[[j]]= d.temp[detect.density[,j]<threshold.density[j],]
            }#endfor        
        }#end else
        
    }#end group n.vit>1 else
    
    return(list(alarm=alarm,detect.density=detect.density,threshold.density=threshold.density,detect.quant=detect.quant))
}



evaluate.time.pop = function(alarm,possible.hits,n.vit=12,vitdropday){

    if(n.vit==1){
        alarm=alarm[[1]][alarm[[1]]$julian<=vitdropday,]

                # TRUE POSITIVE:  groundtruth data  says it's an anomaly and so algorithm does.
        tp = sum(alarm$dropped,na.rm=TRUE)
            
            # FALSE POSITIVE:  groundtruth data says it's not an anomaly, but algorithm says anomaly.
        fp = sum(alarm$dropped==0,na.rm=TRUE)
            
        # FALSE NEGATIVE: groundtruth data says it's an anomaly, but algorithm says not anomaly.
        fn = possible.hits - tp
            
        tp.pop = tp
        fp.pop = fp
        fn.pop = fn
        
        # precision and recall
        out.prec = tp.pop/(tp.pop+fp.pop)
        out.recall = tp.pop/(tp.pop+fn.pop)
        
        # F1 value
        out.F1 = (2*out.prec*out.recall)/(out.prec+out.recall)
    }
    else{
        for(j in 1:n.vit){
            alarm[[j]]=alarm[[j]][alarm[[j]]$julian<=vitdropday[j],]
        }
        
        tp = rep(NA,n.vit)
        fp = rep(NA,n.vit)
        fn = rep(NA,n.vit)
        
        for(j in 1:n.vit){
            
            # TRUE POSITIVE:  groundtruth data  says it's an anomaly and so algorithm does.
            tp[j] = sum(alarm[[j]]$dropped,na.rm=TRUE)
            
            # FALSE POSITIVE:  groundtruth data says it's not an anomaly, but algorithm says anomaly.
            fp[j] = sum(alarm[[j]]$dropped==0,na.rm=TRUE)
            
            # FALSE NEGATIVE: groundtruth data says it's an anomaly, but algorithm says not anomaly.
            fn[j] = possible.hits[j] - tp[j]
            
        }
        
        tp.pop = sum(tp)
        fp.pop = sum(fp)
        fn.pop = sum(fn)
        
        # precision and recall
        out.prec = tp.pop/(tp.pop+fp.pop)
        out.recall = tp.pop/(tp.pop+fn.pop)
        
        # F1 value
        out.F1 = (2*out.prec*out.recall)/(out.prec+out.recall)
    }#endelse
    
    return(list(tp=tp,fp=fp,fn=fn,out.prec=out.prec,out.recall=out.recall,out.F1=out.F1,tp.pop=tp.pop,fp.pop=fp.pop,fn.pop=fn.pop))
    
}


out.ind.5732=anomalyDetect(1,part.window=126,id=c(5732),d[d$lowtag==5732,],eps=.12,31:37)
out.ind.5004=anomalyDetect(1,part.window=126,id=c(5004),d[d$lowtag==5004,],eps=.12,31:37)
evaluate.time.pop(alarm=out.ind.5732$alarm,possible.hits=possible.hits[5],n.vit=1,vitdropday = d.vit$juliandropped[5])
evaluate.time.pop(alarm=out.ind.5004$alarm,possible.hits=possible.hits[1],n.vit=1,vitdropday = d.vit$juliandropped[1])

covs.indx=31:37
epsnum=15
covsnum=length(covs.indx)


df=d[d$lowtag!=individs[j],]
train.cov=array(NA,c(epsnum,3))
for(k in 1:epsnum){# loop over epsilon values
    fit.train = anomalyDetect(n.vit=11,part.window=part.window,id=individs[-j],d=df,eps=rep(k/100,covsnum),covs.indx=31:37)# fit model
    eval.temp=evaluate.time.pop(alarm=fit.train$alarm,possible.hits=possible.hits[-j],n.vit=11,vitdropday = d.vit$juliandropped[-j])
    train.cov[k,1]=eval.temp$out.prec
    train.cov[k,2]=eval.temp$out.recall
    train.cov[k,3]=eval.temp$out.F1# calculate eval = recall
}
which(train.cov,1,max)# find best episilon for each covariate





testnum=15
test=rep(list(),n.vit)
for(j in 1:n.vit){
    
    df=d[d$lowtag!=individs[j],]
    eval=rep(list(),testnum)
    for(k in 1:testnum){
        anomalyDetect(n.vit=11,part.window=part.window,id=id,d=df,eps=k/100,covs.indx=31:37)
        fit.train = anomalyDetect(n.vit=11,part.window=part.window,id=individs[-j],d=data,eps=k/100,covs.indx=31:37)
        eval=evaluate.time.pop(alarm=train$alarm,possible.hits=possible.hits[-j],n.vit=11,vitdropday = d.vit$juliandropped[-j])
    }
    
    
}
for(k in 1:testnum){
    test = anomalyDetect(n.vit=11,part.window=part.window,id=individs,d,eps=k/100,covs.indx=31:37)
    eval[[k]]=evaluate.time.pop(alarm=test$alarm,possible.hits=possible.hits,n.vit=12,vitdropday = d.vit$juliandropped)
}
loo.eval[[j]]=eval[[k]]


ending=Sys.time()
ending-starting



##############################################################################################
###
### Multivariate Anomaly Detection Using Episolon = .12
### Changing predictors
###
##############################################################################################

#possible variables are 31:37
variables_to_fit_Mat = expand.grid(c(TRUE,FALSE), c(TRUE,FALSE),
                                    c(TRUE,FALSE), c(TRUE,FALSE),
                                    c(TRUE,FALSE), c(TRUE,FALSE),
                                    c(TRUE,FALSE))
variables_to_fit_Mat = variables_to_fit_Mat[-(dim(variables_to_fit_Mat)[1]),]# removes pointless line with no predictors
names(variables_to_fit_Mat) = names(d)[31:37]

num.models = dim(variables_to_fit_Mat)[1]
fit.all = matrix(NA,nr = num.models,nc = 6)
names(fit.all) = c("Precision","Recall","F1","tp","fp","fn")

starting=Sys.time()
for(i in 1:num.models){
    test = anomalyDetect(n.vit=12,part.window=part.window,id=individs,d,eps=.12,covs.indx=(30+which(variables_to_fit_Mat[i,]==TRUE)))
    evaluated=evaluate.time.pop(alarm=test$alarm,possible.hits=possible.hits,n.vit=12,vitdropday = d.vit$juliandropped)
    fit.all[i,1] = evaluated$out.prec
    fit.all[i,2] = evaluated$out.recall
    fit.all[i,3] = evaluated$out.F1
    fit.all[i,4] = evaluated$tp.pop
    fit.all[i,5] = evaluated$fp.pop
    fit.all[i,6] = evaluated$fn.pop
}
ending=Sys.time()
ending-starting

par(mfrow=c(3,2))
for(i in 1:6){
    plot(fit.all[,i])
}
plot(fit.all[,1])

checkup=apply(fit.all,2,function(x){which(x==max(x))})
checkup[[6]] = which(fit.all[,6]==min(fit.all[,6]))


variables_to_fit_Mat[checkup[[1]],]
variables_to_fit_Mat[checkup[[2]],]
variables_to_fit_Mat[checkup[[3]],]
variables_to_fit_Mat[checkup[[4]],]
variables_to_fit_Mat[checkup[[5]],]
variables_to_fit_Mat[checkup[[6]],]

#finds the 5 biggest values, want to minimize FP/FN rates
checkup=apply(fit.all,2,function(x){tail(order(x),5)})

checkup[,5] = head(order(fit.all[,5]),5)
checkup[,6] =  head(order(fit.all[,6]),5)


ranking=rbind(apply(variables_to_fit_Mat[checkup[,1],],2,sum),apply(variables_to_fit_Mat[checkup[,2],],2,sum),
              apply(variables_to_fit_Mat[checkup[,3],],2,sum),apply(variables_to_fit_Mat[checkup[,4],],2,sum),
              apply(variables_to_fit_Mat[checkup[,5],],2,sum),apply(variables_to_fit_Mat[checkup[,6],],2,sum)
        )

apply(ranking,2,sum)

ranking2=rbind(apply(variables_to_fit_Mat[checkup[,1],],2,sum),apply(variables_to_fit_Mat[checkup[,2],],2,sum),
              apply(variables_to_fit_Mat[checkup[,3],],2,sum)
)

apply(ranking2,2,sum)

names(variables_to_fit_Mat)[c(TRUE,TRUE,TRUE,FALSE,TRUE,FALSE,FALSE)] #final model?
final.out = anomalyDetect(n.vit=12,part.window=part.window,id=individs,d,eps=.10,covs.indx=c(31,32,33,35))
evaluated=evaluate.time.pop(alarm=final.out$alarm,possible.hits=possible.hits,n.vit=12,vitdropday = d.vit$juliandropped)
final.out
evaluated


par(mfrow=c(4,1))
for(j in 1:n.vit){
    # plot(d$julian[d$lowtag==individs[j]],d$logstep[d$lowtag==individs[j]],main=d.vit$lowtag[j])
    # points(final.out$alarm[[j]]$julian,final.out$alarm[[j]]$logstep,col=2)
    # abline(v=d.vit$juliandropped[j])
    plot(d$julian[d$lowtag==individs[j]],d$scalealt[d$lowtag==individs[j]],main=d.vit$lowtag[j])
    points(final.out$alarm[[j]]$julian,final.out$alarm[[j]]$scalealt,col=2)
    abline(v=d.vit$juliandropped[j])

    plot(d$julian[d$lowtag==individs[j]],d$scaletemp[d$lowtag==individs[j]],main=d.vit$lowtag[j])
    points(final.out$alarm[[j]]$julian,final.out$alarm[[j]]$scaletemp,col=2)
    abline(v=d.vit$juliandropped[j])
    
    plot(d$julian[d$lowtag==individs[j]],d$scaleangle[d$lowtag==individs[j]],main=d.vit$lowtag[j])
    points(final.out$alarm[[j]]$julian,final.out$alarm[[j]]$scaleangle,col=2)
    abline(v=d.vit$juliandropped[j])
    
    plot(d$julian[d$lowtag==individs[j]],d$scaledx[d$lowtag==individs[j]],main=d.vit$lowtag[j])
    points(final.out$alarm[[j]]$julian,final.out$alarm[[j]]$scaledx,col=2)
    abline(v=d.vit$juliandropped[j])
    }


names(variables_to_fit_Mat) #final model?
final.out2 = anomalyDetect(n.vit=12,part.window=part.window,id=individs,d,eps=.10,covs.indx=c(31,32,35,37))
evaluated2=evaluate.time.pop(alarm=final.out$alarm,possible.hits=possible.hits,n.vit=12,vitdropday = d.vit$juliandropped)
final.out2
evaluated
evaluated2


par(mfrow=c(4,1))
for(j in 1:n.vit){
    # plot(d$julian[d$lowtag==individs[j]],d$logstep[d$lowtag==individs[j]],main=d.vit$lowtag[j])
    # points(final.out2$alarm[[j]]$julian,final.out2$alarm[[j]]$logstep,col=2)
    # abline(v=d.vit$juliandropped[j])
    plot(d$julian[d$lowtag==individs[j]],d$scalealt[d$lowtag==individs[j]],main=d.vit$lowtag[j])
    points(final.out2$alarm[[j]]$julian,final.out2$alarm[[j]]$scalealt,col=2)
    abline(v=d.vit$juliandropped[j])
    
    plot(d$julian[d$lowtag==individs[j]],d$scaletemp[d$lowtag==individs[j]],main=d.vit$lowtag[j])
    points(final.out2$alarm[[j]]$julian,final.out2$alarm[[j]]$scaletemp,col=2)
    abline(v=d.vit$juliandropped[j])
    
    plot(d$julian[d$lowtag==individs[j]],d$scaledx[d$lowtag==individs[j]],main=d.vit$lowtag[j])
    points(final.out2$alarm[[j]]$julian,final.out2$alarm[[j]]$scaledx,col=2)
    abline(v=d.vit$juliandropped[j])
    
    plot(d$julian[d$lowtag==individs[j]],d$logscalestep[d$lowtag==individs[j]],main=d.vit$lowtag[j])
    points(final.out2$alarm[[j]]$julian,final.out2$alarm[[j]]$logscalestep,col=2)
    abline(v=d.vit$juliandropped[j])
    
}

###
### If only use Recall/sensitivity as criteria
###

apply(variables_to_fit_Mat[checkup[,2],],2,sum)
par(mfrow=c(1,1))
plot(fit.all[,2])
points(checkup[,2],fit.all[checkup[,2],2],col=2)
variables_to_fit_Mat[checkup[,2],]#scaletemp,scaleanlge,scaledx

#model3
names(variables_to_fit_Mat)
final.out3 = anomalyDetect(n.vit=12,part.window=part.window,id=individs,d,eps=.1,covs.indx=c(32,33,35))
evaluated3=evaluate.time.pop(alarm=final.out$alarm,possible.hits=possible.hits,n.vit=12,vitdropday = d.vit$juliandropped)
evaluated
evaluated3


par(mfrow=c(4,1))
for(j in 1:n.vit){
    # plot(d$julian[d$lowtag==individs[j]],d$logstep[d$lowtag==individs[j]],main=d.vit$lowtag[j])
    # points(final.out2$alarm[[j]]$julian,final.out2$alarm[[j]]$logstep,col=2)
    # abline(v=d.vit$juliandropped[j])
    plot(d$julian[d$lowtag==individs[j]],d$scalealt[d$lowtag==individs[j]],main=d.vit$lowtag[j])
    points(final.out2$alarm[[j]]$julian,final.out2$alarm[[j]]$scalealt,col=2)
    abline(v=d.vit$juliandropped[j])
    
    plot(d$julian[d$lowtag==individs[j]],d$scaletemp[d$lowtag==individs[j]],main=d.vit$lowtag[j])
    points(final.out2$alarm[[j]]$julian,final.out2$alarm[[j]]$scaletemp,col=2)
    abline(v=d.vit$juliandropped[j])
    
    plot(d$julian[d$lowtag==individs[j]],d$scaledx[d$lowtag==individs[j]],main=d.vit$lowtag[j])
    points(final.out2$alarm[[j]]$julian,final.out2$alarm[[j]]$scaledx,col=2)
    abline(v=d.vit$juliandropped[j])
    
    plot(d$julian[d$lowtag==individs[j]],d$logscalestep[d$lowtag==individs[j]],main=d.vit$lowtag[j])
    points(final.out2$alarm[[j]]$julian,final.out2$alarm[[j]]$logscalestep,col=2)
    abline(v=d.vit$juliandropped[j])
    
}


# ##############################################################################################
# ###
# ### Multivariate Anomaly Detection within a function, and using evalaute funcition to optimize epsilon
# ###
# ##############################################################################################
# 
# par(mfrow=c(2,1))
# 
# for(j in 1:n.vit){
#     plot(alarm[[j]]$julian,alarm[[j]]$logstep)
#     abline(v=d.vit$juliandropped[j])
# }
# par(mfrow=c(2,1))
# for(j in 1:n.vit){
#     plot(d$julian[d$lowtag==individs[j]],d$logstep[d$lowtag==individs[j]])
#     points(alarm[[j]]$julian,alarm[[j]]$logstep,col=2)
#     abline(v=d.vit$juliandropped[j])
#     
# }
# 
# 
# 
# 
# 
# ###################################################################################################
# ###
# ### Examining Time of Day, to see if movements change more depending on time of day
# ###
# ###################################################################################################
# 
# 
# d$hour=hour(d$date_time_local)
# par(mfrow=c(1,1))
# for(j in 1:n.vit){
#     plot(d$hour[d$lowtag==individs[j]],d$logstep[d$lowtag==individs[j]])
# }
# 
# ###
# ### selection threshold
# ###
# 
# # yval = groundtruth
# #pval = validation set
# 
# selectThreshold=function(yval, pval){
#     
#     bestEpsilon = 0
#     bestF1 = 0
#     F1 = 0
#     
#     check=seq(min(pval),max(pval),by=1000)
#     
#     for (epsilon in check){
#         predictions = (pval < epsilon)
#         tp = sum(predictions == 1 & yval == 1)
#         fp = sum(predictions == 1 & yval == 0)
#         fn = sum(predictions == 0 & yval == 1)
#         prec = tp/(tp+fp)
#         rec = tp/(tp+fn)
#         F1 = (2*prec*rec)/(prec+rec)
#         if (F1 > bestF1){
#             bestF1 = F1
#             bestEpsilon = epsilon   
#         }
#     }
#     
#     return(list(bestEpsilon=bestEpsilon, bestF1=bestF1)
# }
# 
# 
# # par(mfrow=c(3,1))
# # for(j in 1:n.vit){
# #     plot(julian.out[,j],sig.distance[,j])
# #     abline(v=d.vit$juliandropped[j])
# #     hist(sig.distance[,j])
# #     plot(density(sig.distance[,j]),na.rm=TRUE)
# #     abline(v=min(sig.distance[,j]),col=2)
# # }
# 
# 
# min(d$julian)
# #which index starts each individual's run?
# 
# time.min=c(1,rep(NA,n.vit-1),dim(d)[1])
# for(j in 2:n.vit){
#     time.min[j]=time.min[j-1]+n.obs[j-1]
# }
# time.min
# par(mfrow=c(2,2))
# for(j in 1:n.vit){
#     d.temp=d[d$lowtag==individs[j],]
#     plot(density(d.temp$step[d.temp$julian<126]))
#     
# }
# 
# for(i in time.min[j]:time.min[j+1])
#     mean(d$distance[d$lowtag==individs[j]])
# 
# +
#     d.temp$step[d.temp$julian[i,j]<126]
# d[which(!is.finite(d$step)),]
# 
# 
# 
# for(j in 1:n.vit){
#     d$step[d$lowtag==individs[j]]
# }
# 

#####################################################################################################
###
### Changepoint analysis
###
#####################################################################################################


# apply(angle.dist.corr,2,sd,na.rm=TRUE)
# 
# sd(angle.dist.corr,na.rm=TRUE)
# par(mfrow=c(1,1))
# change.out=list(rep(list(),n.vit))
# birth.pred.indx=c()
# birth.pred=c()
# for(i in 1:n.vit){
#     change.out[[i]]=cpt.var(angle.dist.corr[i,!is.na(angle.dist.corr[i,])])
#     plot(change.out[[i]])
#     abline(v=head(vit.records[[i]],1),col="darkblue")
#     abline(v=tail(vit.records[[i]],1),col="darkblue")
#     birth.pred.indx[i]=cpts(change.out[[i]])
#     birth.pred[i]=d$julian[d$lowtag==individs[i]][birth.pred.indx[i]]
# }


#####################################################################################################
###
### Supervised Learning
###
###
#####################################################################################################


mu.dist.acf=apply(mu.distance[1,1:n.obs[1]],acf,na.rm=TRUE)


which(is.na(mu.distance),arr.ind = T)
par(mfrow=c(3,2))
for(j in 1:n.vit){
    acf(d$dist[d$lowtag==individs[j]],main=individs[j])
}



###
### naiveBayes classifier
###

#possible predictors
# vit.drop
# julian.out
# 
# #min/max features
# min.distance
# max.angle
# 
# #mean features
# mu.distance
# mu.turn.angle
# mu.bearing
# mu.altitude
# mu.temp
# mu.dx
# mu.dy
# 
# #variation features
# sig.distance
# sig.turn.angle
# sig.bearing
# sig.r2n
# sig.altitude
# sig.temp
# sig.dx
# sig.dy
# 
# #correlation features
# corr.dist.angle
# corr.dist.bear
# corr.dist.alt
# corr.angle.bear
# corr.angle.alt
# corr.bear.alt

#by individual
vit.drop=apply(vit.drop,2,as.factor)
pred.nb=rep(list(),n.vit)

for(j in 2:n.vit){
    train=data.frame(cbind(mu.distance[,j],mu.turn.angle[,j]))
    test=data.frame(cbind(mu.distance[,1],mu.turn.angle[,1]))
    names(train)=names(test)=c('mu_distance','mu_turn_angle')
    head(test)
    out.nb[[j]] = naiveBayes(train,vit.drop[,j])
    pred.nb[[j]]=predict(out.nb[[j]],test,type="raw")
}

###
### Population level
###

head(corr.dist.angle)
length(c(t(corr.dist.angle)))

# 
# df.pop=matrix(NA,dim(mu.distance)[1]*dim(mu.distance)[2],2+2+7+8+6)
# 
# df.pop[,1]=c(vit.drop)
# df.pop[,2]=c(julian.out)
# # 
# # #min/max features
# df.pop[,3]=c(t(min.distance))
# df.pop[,4]=c(t(max.angle))
#  
# # #mean features
# df.pop[,5]=c(t(mu.distance))
# df.pop[,6]=c(t(mu.turn.angle))
# df.pop[,7]=c(t(mu.bearing))
# df.pop[,8]=c(t(mu.altitude))
# df.pop[,9]=c(t(mu.temp))
# df.pop[,10]=c(t(mu.dx))
# df.pop[,11]=c(t(mu.dy))
# # 
# # #variation features
# df.pop[,12]=c(t(sig.distance))
# df.pop[,13]=c(t(sig.turn.angle))
# df.pop[,14]=c(t(sig.bearing))
# df.pop[,15]=c(t(sig.r2n))
# df.pop[,16]=c(t(sig.altitude))
# df.pop[,17]=c(t(sig.temp))
# df.pop[,18]=c(t(sig.dx))
# df.pop[,19]=c(t(sig.dy))
# 
# #correlation features
# df.pop[,20]=c(t(corr.dist.angle))
# df.pop[,21]=c(t(corr.dist.bear))
# df.pop[,22]=c(t(corr.dist.alt))
# df.pop[,23]=c(t(corr.angle.bear))
# df.pop[,24]=c(t(corr.angle.alt))
# df.pop[,25]=c(t(corr.bear.alt))

                  
###
### naiveBayes across windows
###

# vit.drop
# julian.out
# 
# #min/max features
 
# #mean features
# mu.distance
# mu.turn.angle
# mu.bearing
# mu.altitude
# mu.temp
# mu.dx
# mu.dy
# 
# #variation features
# sig.distance
# sig.turn.angle
# sig.bearing
# sig.r2n
# sig.altitude
# sig.temp
# sig.dx
# sig.dy
# 
# #correlation features
# corr.dist.angle
# corr.dist.bear
# corr.dist.alt
# corr.angle.bear
# corr.angle.alt
# corr.bear.alt

vit.drop=apply(vit.drop,2,as.factor)
# replace.01=function(x){
#     ifelse(x==0,'n','y')
# }
# 
# for(i in 1:dim(vit.drop)[1]){
#     for(j in 1:dim(vit.drop)[2]){
#         vit.drop[i,j]=replace.01(vit.drop[i,j])
#     }
# }

###
### logistic regression
###
window=dim(mu.distance)[1]
out.glm=rep(list(),window)
pred.glm=rep(list(),window)

for(w in 1:window){
    train=data.frame(cbind(vit.drop[w,-1],mu.distance[w,-1],mu.turn.angle[w,-1]))
    test=data.frame(cbind(vit.drop[w,1],mu.distance[w,1],mu.turn.angle[w,1]))
    names(train)=names(test)=c('vit.drop','mu_distance','mu_turn_angle')
    out.glm[[w]]=glm(vit.drop~.,data=train,family="binomial")
    pred.glm[[w]]=predict(out.glm[[w]], test, type="response")
}
w

out.glm[[w]]


plot(as.numeric(pred.glm))







window=dim(mu.distance)[2]
out.nb=rep(list(),window)
pred.nb=rep(list(),window)
fitit=c()
for(w in 1:window){
        train=data.frame(cbind(mu.distance[w,],mu.turn.angle[w,]))
        test=data.frame(cbind(vit.drop[w,],mu.distance[w,],mu.turn.angle[w,]))
        names(train)=names(test)=c('mu_distance','mu_turn_angle')
        svm(x=train,y=vit.drop[w,],data=train)
        
        out.nb[[w]] = naiveBayes(train,as.factor(vit.drop[w,-1]))
        pred.nb[[w]]=predict(out.nb[[w]],list(test))
        fitit=c(fitit,w)
        
        svm()
    }

}
fitit
out.nb    
naiveBayes(train,as.factor(vit.drop[w,-1]))

head(pred.nb)






###
### Random Forest
###

library(randomForest)
names(d)

d.rf=d[,c(15,31:37)]
scalealt_imput=rfImpute(dropped ~ ., scalealt,data=d.rf)

which(is.na(d[31:37]),arr.ind = TRUE)

# Create the forest.

output.forest <- randomForest(dropped ~ scalealt+scaletemp+
                              scaleangle+scalebear+scaledx+logscalestep, 
                              data = d)

