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
library(caret)
library(mvtnorm)
library(raster)
library(dplyr)
library(randomForest)
library(randomForestExplainer)


setwd("~/Documents/Parturition/180108_parturition")


###
### Load VIT data
###

d.vit=mdb.get('~/Documents/Data/SWDPPdeerDB.MDB',tables= "VIT")
names(d.vit)=tolower(gsub('[[:punct:]]',"",names(d.vit)))
names(d.vit)[10]="lowtag"
d.vit$datedropped=as.character(d.vit$datedropped)
d.vit$juliandropped=yday(mdy_hms(d.vit$datedropped))
head(d.vit)
d.vit$datefawnfound=as.character(d.vit$datefawnfound)
d.vit$julianfawnfound=yday(mdy_hms(d.vit$datefawnfound))

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
#for loop for reading in data, using vector of lowtag's from the vit dataset
for(i in 1:n.vit){
    d.temp = read.table(paste("/home/aketz/Documents/Data/GPS_locations_ind/",d.vit$lowtag[i],".csv",sep=""),sep=",",header=TRUE,stringsAsFactors = FALSE)
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

###
### subset entire dataset to parturition window
###

d=d[d$julian>start & d$julian <= end,]

#removing last day of parturition of 5004, outlier movements
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


###
### double checking for duplicate locations
###

summary(duplicated(d))

###
### Adding Vit data to main GPS dataframe
###

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



#Create time lag between successive locations to censor data if needed.
time.diff <- diff(d$date_time_local)
d=d[-1,]
d$timediff <-round(as.numeric(abs(time.diff)))
rm=which(d$timediff>10)
d=d[-rm,]
names(d)[1]="devname"

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
### removed those points up above
# plot(d.traj[1])

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
### Description of parameters returned by adehabitatLT object model. 
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


########################################################################################################
###
### remove the initial index of each individual 
### because it makes no sense to calculate distances 
### between the locations of the individuals
###
#######################################################################################################

remove.indx=which(is.na(d$dist.traj))
d[remove.indx,]
d=d[-remove.indx,]

#######################################################################################################
###
### Adding Landcover data
###
#######################################################################################################3

### Convert primary GPS location data back to spatial points data frame

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

# extent of points data
extent.d.sp=extent(d.sp.proj)
d.sp.crs = crs(d.sp.proj)

###
### Import NLCD
###

#load cover data
nlcd = raster("~/Documents/Data/nlcd_2011_landcover_2011_edition_2014_10_10/nlcd_2011_landcover_2011_edition_2014_10_10.img")

# convert lat/lon to appropriate projection
crs_args = nlcd@crs@projargs
d.sp.proj.nlcd=spTransform(d.sp.proj,crs_args)

#clip the NLCD to the d.sp.proj size + small 1% buffer
b = bbox(d.sp.proj.nlcd)
b[,1]=b[,1]*.99
b[,2]=b[,2]*1.01
nlcd.crop=crop(nlcd,b)

#plot of nlcd with movements of deer on top
par(mfrow=c(1,1))
plot(nlcd.crop)
plot(d.sp.proj.nlcd,add=T)    

#extract land cover data for each point, given buffer size
landcover = extract(nlcd.crop, d.sp.proj.nlcd)

# generate land cover number to name conversions
num.codes = unique(unlist(landcover))
cover.names = nlcd@data@attributes[[1]]$NLCD.2011.Land.Cover.Class[num.codes + 1]
levels(cover.names)[1] = NA 
conversions = data.frame(num.codes, cover.names)
conversions = na.omit(conversions)
conversions = conversions[order(conversions$num.codes),]
conversions

# # summarize each site's data by proportion of each cover type
summ = lapply(landcover, function(x){prop.table(table(x))})

# convert to data frame
nlcd.points = data.frame(d.sp.proj.nlcd,cover2 = names(unlist(summ)))

# create cover name column
nlcd.points$cover = nlcd.points$cover2
levels(nlcd.points$cover) = conversions$cover.names
table(nlcd.points$cover)
coords=coordinates(nlcd.points$longitude,nlcd.points$latitude)

#convert and project nlcd.points back to spatial points data frame, and then back to just data frame
coords = cbind(nlcd.points$longitude, nlcd.points$latitude)
sp = SpatialPoints(coords)

# make spatial data frame
# spdf = SpatialPointsDataFrame(coords, d)
spdf = SpatialPointsDataFrame(sp, nlcd.points)

# EPSG strings
latlong = "+init=epsg:4326"
proj4string(spdf) = CRS(latlong)
d.nlcd = spTransform(spdf, CRS("+proj=tmerc +lat_0=0 +lon_0=-90 +k=0.9996 +x_0=520000
                                  +y_0=-4480000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
d=data.frame(d.nlcd)
rm.columns=c(29:32,34:36)
d=d[,-rm.columns]

head(d)
levels(d$cover)

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

par(mfrow=c(3,2))
for(j in 1:n.vit){
    plot(julian.out[1:(length(julian.out[,j])-7),j],sig.distance[,j],main=individs[j])
    abline(v=d.vit$juliandropped[j],col=2)
    abline(v=julian.out[which(sig.distance[,j]==min(sig.distance[1:200,j],na.rm=TRUE)),j])
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

names(d)[32:39]
par(mfrow=c(2,2))
for(j in 30:37){
    hist(d[,j],breaks=50,main=names(d)[j])
}

#normal enough: altitude, temp, dx,r2n,logscalestep 
#bad: scaleangle, bearing

par(mfrow=c(3,1))
for(j in 1:n.vit){
    hist(log(d$step[d$lowtag==individs[j]]),breaks=50,main=individs[j])
    plot(density(log(d$step[d$lowtag==individs[j]])),main=individs[j])
    plot(d$julian[d$lowtag==individs[j]],log(d$step[d$lowtag==individs[j]]))
    abline(v=d.vit$juliandropped[j],col=2)
}


###
### Paturition window
###

part.window=yday("2017-05-06")

###
### Total number of possible anomalies to detect
###

possible.hits=rep(NA,n.vit)
for(j in 1:n.vit){
    possible.hits[j] = sum(d$dropped[d$lowtag==individs[j]])
}
possible.hits


##############################################################################################
###
### Expanding Anomaly Detection function to deal 
### with multiple individuals as well as single
###
##############################################################################################
pw=part.window
ph=possible.hits

source("anomalyDetectUniv.R")
source("anomalyDetect.R")
source("evaluate_ad.R")

out.test=anomalyDetect(n.vit=12,part.window=126,id=individs,d=d,eps=.12,covs.indx=32:39)
    
out.ind.5732=anomalyDetect(1,part.window=pw,id=c(5732),d[d$lowtag==5732,],eps=.12,32:39)
out.ind.5004=anomalyDetect(1,part.window=126,id=c(5004),d[d$lowtag==5004,],eps=.12,32:39)
evaluate.time.pop(alarm=out.ind.5732$alarm,possible.hits=possible.hits[5],n.vit=1,vitdropday = d.vit$juliandropped[5])
evaluate.time.pop(alarm=out.ind.5004$alarm,possible.hits=possible.hits[1],n.vit=1,vitdropday = d.vit$juliandropped[1])

covs.indx=32:39
epsnum=15
covsnum=length(covs.indx)
loo.eval=rep(list(),3*n.vit)
ci=32:39
starting=Sys.time()
for(m in 1:3){
    for(j in 1:n.vit){
        df=d[d$lowtag!=individs[j],]
        train.cov=array(NA,c(covsnum,epsnum,3))
        for(h in 1:covsnum){
            for(k in 1:epsnum){# loop over epsilon values
                fit.train = anomalyDetect(n.vit=11,part.window=pw,id=individs[-j],d=df,eps=k/100,covs.indx=ci[h])# fit model
                eval.temp=evaluate.time.pop(alarm=fit.train$alarm,possible.hits=ph[-j],n.vit=11,vitdropday = d.vit$juliandropped[-j])
                train.cov[h,k,1]=eval.temp$out.prec
                train.cov[h,k,2]=eval.temp$out.recall
                train.cov[h,k,3]=eval.temp$out.F1# calculate eval = recall
            }
        }
        compare=apply(train.cov,c(1,3),max)
        decide=matrix(NA,nr=covsnum,ncol=3)
        for(h in 1:covsnum){
            for(i in 1:3){
                decide[h,i]=min(which(train.cov[h,,i]==compare[h,i]))
            }
        }
        decide=decide/100
        d.ind=d[d$lowtag==individs[j],]
        fit.test= anomalyDetect(n.vit=1,part.window=pw,id=individs[j],d=d.ind,eps=decide[,m],covs.indx=ci)# fit model
        eval.test=evaluate.time.pop(alarm=fit.test$alarm,possible.hits=ph[j],n.vit=1,vitdropday = d.vit$juliandropped[j])
        loo.indx=j+(m-1)*12
        loo.eval[[loo.indx]]=eval.test
    }    
}
ending=Sys.time()
total.time=ending-starting
total.time


check.prec=rep(NA,3*n.vit)
check.recall=rep(NA,3*n.vit)
check.F1=rep(NA,3*n.vit)
check.TP=rep(NA,3*n.vit)
check.FP=rep(NA,3*n.vit)
check.FN=rep(NA,3*n.vit)

loo.eval[[25:36]]

loo.eval
length(loo.eval[[25]])

for(j in 1:(3*n.vit)){
    check.prec[j]=loo.eval[[j]]$out.prec
    check.recall[j]=loo.eval[[j]]$out.recall
    check.F1[j]=loo.eval[[j]]$out.F1
    check.TP[j]=loo.eval[[j]]$tp
    check.FP[j]=loo.eval[[j]]$fp
    check.FN[j]=loo.eval[[j]]$fn
    }
    

par(mfrow=c(3,1))
test.crit=c(rep(1,n.vit),rep(2,n.vit),rep(3,n.vit))
test.ind=rep(1:n.vit,3)
pdf("Anomaly_eval_1.pdf")
plot(test.ind,check.prec,col=test.crit,main="Precision",ylab="Precision",xlab="Individual")
legend(x=5,y=.23,fill=c(1:3),legend=c("Precision","Recall","F1"))
plot(test.ind,check.recall,col=test.crit,main="Recall",ylab="Recall",xlab="Individual")
plot(test.ind,check.F1,col=test.crit,main="F1",ylab="F1",xlab="Individual")
dev.off()
pdf("Anomaly_eval_1.pdf")
plot(test.ind,check.TP,col=test.crit,main="TP",ylab="TP",xlab="Individual")
plot(test.ind,check.FP,col=test.crit,main="FP",ylab="FP",xlab="Individual")
plot(test.ind,check.FN,col=test.crit,main="FN",ylab="FN",xlab="Individual")
legend("topright",fill=c(1:3),legend=c("Precision","Recall","F1"))
dev.off()


save.image("anomaly.Rdata")

par(mfrow=c(3,1))

pdf("checks/Precision_check.pdf")
plot(check.prec,main="Precision when different criteria")
abline(v=c(12.5,24.5))
text(5,.25,"Precision")
text(18,.25,"Recall")
text(28,.25,"F1")
dev.off()

pdf("checks/Recall_check.pdf")
plot(check.recall,main="Recall when different criteria")
abline(v=c(12.5,24.5))
text(5,.8,"Precision")
text(18,.8,"Recall")
text(28,.8,"F1")
dev.off()

pdf("checks/F1_check.pdf")
plot(check.F1,main="F1 when different criteria")
abline(v=c(12.5,24.5))
text(5,.25,"Precision")
text(18,.25,"Recall")
text(28,.25,"F1")
dev.off()

pdf("checks/TP_check.pdf")
plot(check.TP,main="TP when diff criteria")
abline(v=c(12.5,24.5))
text(5,10,"Precision")
text(18,10,"Recall")
text(28,10,"F1")
dev.off()

pdf("checks/FP_check.pdf")
plot(check.FP,main="FP")
abline(v=c(12.5,24.5))
text(5,120,"Precision")
text(18,120,"Recall")
text(28,120,"F1")
dev.off()

pdf("checks/FN_check.pdf")
plot(check.FN,main="FN")
abline(v=c(12.5,24.5))
text(5,11,"Precision")
text(18,11,"Recall")
text(28,11,"F1")
dev.off()


# 
# ##############################################################################################
# ###
# ### Multivariate Anomaly Detection Using Episolon = .12
# ### Changing predictors
# ###
# ##############################################################################################
# 
# #possible variables are 31:37
# variables_to_fit_Mat = expand.grid(c(TRUE,FALSE), c(TRUE,FALSE),
#                                    c(TRUE,FALSE), c(TRUE,FALSE),
#                                    c(TRUE,FALSE), c(TRUE,FALSE),
#                                    c(TRUE,FALSE))
# variables_to_fit_Mat = variables_to_fit_Mat[-(dim(variables_to_fit_Mat)[1]),]# removes pointless line with no predictors
# names(variables_to_fit_Mat) = names(d)[31:37]
# 
# num.models = dim(variables_to_fit_Mat)[1]
# fit.all = matrix(NA,nr = num.models,nc = 6)
# names(fit.all) = c("Precision","Recall","F1","tp","fp","fn")
# 
# starting=Sys.time()
# for(i in 1:num.models){
#     test = anomalyDetect(n.vit=12,part.window=part.window,id=individs,d,eps=.12,covs.indx=(30+which(variables_to_fit_Mat[i,]==TRUE)))
#     evaluated=evaluate.time.pop(alarm=test$alarm,possible.hits=possible.hits,n.vit=12,vitdropday = d.vit$juliandropped)
#     fit.all[i,1] = evaluated$out.prec
#     fit.all[i,2] = evaluated$out.recall
#     fit.all[i,3] = evaluated$out.F1
#     fit.all[i,4] = evaluated$tp.pop
#     fit.all[i,5] = evaluated$fp.pop
#     fit.all[i,6] = evaluated$fn.pop
# }
# ending=Sys.time()
# ending-starting
# 
# par(mfrow=c(3,2))
# for(i in 1:6){
#     plot(fit.all[,i])
# }
# plot(fit.all[,1])
# 
# checkup=apply(fit.all,2,function(x){which(x==max(x))})
# checkup[[6]] = which(fit.all[,6]==min(fit.all[,6]))
# 
# 
# variables_to_fit_Mat[checkup[[1]],]
# variables_to_fit_Mat[checkup[[2]],]
# variables_to_fit_Mat[checkup[[3]],]
# variables_to_fit_Mat[checkup[[4]],]
# variables_to_fit_Mat[checkup[[5]],]
# variables_to_fit_Mat[checkup[[6]],]
# 
# #finds the 5 biggest values, want to minimize FP/FN rates
# checkup=apply(fit.all,2,function(x){tail(order(x),5)})
# 
# checkup[,5] = head(order(fit.all[,5]),5)
# checkup[,6] =  head(order(fit.all[,6]),5)
# 
# 
# ranking=rbind(apply(variables_to_fit_Mat[checkup[,1],],2,sum),apply(variables_to_fit_Mat[checkup[,2],],2,sum),
#               apply(variables_to_fit_Mat[checkup[,3],],2,sum),apply(variables_to_fit_Mat[checkup[,4],],2,sum),
#               apply(variables_to_fit_Mat[checkup[,5],],2,sum),apply(variables_to_fit_Mat[checkup[,6],],2,sum)
# )
# 
# apply(ranking,2,sum)
# 
# ranking2=rbind(apply(variables_to_fit_Mat[checkup[,1],],2,sum),apply(variables_to_fit_Mat[checkup[,2],],2,sum),
#                apply(variables_to_fit_Mat[checkup[,3],],2,sum)
# )
# 
# apply(ranking2,2,sum)
# 
# names(variables_to_fit_Mat)[c(TRUE,TRUE,TRUE,FALSE,TRUE,FALSE,FALSE)] #final model?
# final.out = anomalyDetect(n.vit=12,part.window=part.window,id=individs,d,eps=.10,covs.indx=c(31,32,33,35))
# evaluated=evaluate.time.pop(alarm=final.out$alarm,possible.hits=possible.hits,n.vit=12,vitdropday = d.vit$juliandropped)
# final.out
# evaluated
# 
# 
# par(mfrow=c(4,1))
# for(j in 1:n.vit){
#     # plot(d$julian[d$lowtag==individs[j]],d$logstep[d$lowtag==individs[j]],main=d.vit$lowtag[j])
#     # points(final.out$alarm[[j]]$julian,final.out$alarm[[j]]$logstep,col=2)
#     # abline(v=d.vit$juliandropped[j])
#     plot(d$julian[d$lowtag==individs[j]],d$scalealt[d$lowtag==individs[j]],main=d.vit$lowtag[j])
#     points(final.out$alarm[[j]]$julian,final.out$alarm[[j]]$scalealt,col=2)
#     abline(v=d.vit$juliandropped[j])
#     
#     plot(d$julian[d$lowtag==individs[j]],d$scaletemp[d$lowtag==individs[j]],main=d.vit$lowtag[j])
#     points(final.out$alarm[[j]]$julian,final.out$alarm[[j]]$scaletemp,col=2)
#     abline(v=d.vit$juliandropped[j])
#     
#     plot(d$julian[d$lowtag==individs[j]],d$scaleangle[d$lowtag==individs[j]],main=d.vit$lowtag[j])
#     points(final.out$alarm[[j]]$julian,final.out$alarm[[j]]$scaleangle,col=2)
#     abline(v=d.vit$juliandropped[j])
#     
#     plot(d$julian[d$lowtag==individs[j]],d$scaledx[d$lowtag==individs[j]],main=d.vit$lowtag[j])
#     points(final.out$alarm[[j]]$julian,final.out$alarm[[j]]$scaledx,col=2)
#     abline(v=d.vit$juliandropped[j])
# }
# 
# 
# names(variables_to_fit_Mat) #final model?
# final.out2 = anomalyDetect(n.vit=12,part.window=part.window,id=individs,d,eps=.10,covs.indx=c(31,32,35,37))
# evaluated2=evaluate.time.pop(alarm=final.out$alarm,possible.hits=possible.hits,n.vit=12,vitdropday = d.vit$juliandropped)
# final.out2
# evaluated
# evaluated2
# 
# 
# par(mfrow=c(4,1))
# for(j in 1:n.vit){
#     # plot(d$julian[d$lowtag==individs[j]],d$logstep[d$lowtag==individs[j]],main=d.vit$lowtag[j])
#     # points(final.out2$alarm[[j]]$julian,final.out2$alarm[[j]]$logstep,col=2)
#     # abline(v=d.vit$juliandropped[j])
#     plot(d$julian[d$lowtag==individs[j]],d$scalealt[d$lowtag==individs[j]],main=d.vit$lowtag[j])
#     points(final.out2$alarm[[j]]$julian,final.out2$alarm[[j]]$scalealt,col=2)
#     abline(v=d.vit$juliandropped[j])
#     
#     plot(d$julian[d$lowtag==individs[j]],d$scaletemp[d$lowtag==individs[j]],main=d.vit$lowtag[j])
#     points(final.out2$alarm[[j]]$julian,final.out2$alarm[[j]]$scaletemp,col=2)
#     abline(v=d.vit$juliandropped[j])
#     
#     plot(d$julian[d$lowtag==individs[j]],d$scaledx[d$lowtag==individs[j]],main=d.vit$lowtag[j])
#     points(final.out2$alarm[[j]]$julian,final.out2$alarm[[j]]$scaledx,col=2)
#     abline(v=d.vit$juliandropped[j])
#     
#     plot(d$julian[d$lowtag==individs[j]],d$logscalestep[d$lowtag==individs[j]],main=d.vit$lowtag[j])
#     points(final.out2$alarm[[j]]$julian,final.out2$alarm[[j]]$logscalestep,col=2)
#     abline(v=d.vit$juliandropped[j])
#     
# }
# 
# ###
# ### If only use Recall/sensitivity as criteria
# ###
# 
# apply(variables_to_fit_Mat[checkup[,2],],2,sum)
# par(mfrow=c(1,1))
# plot(fit.all[,2])
# points(checkup[,2],fit.all[checkup[,2],2],col=2)
# variables_to_fit_Mat[checkup[,2],]#scaletemp,scaleanlge,scaledx
# 
# #model3
# names(variables_to_fit_Mat)
# final.out3 = anomalyDetect(n.vit=12,part.window=part.window,id=individs,d,eps=.1,covs.indx=c(32,33,35))
# evaluated3=evaluate.time.pop(alarm=final.out$alarm,possible.hits=possible.hits,n.vit=12,vitdropday = d.vit$juliandropped)
# evaluated
# evaluated3
# 
# 
# par(mfrow=c(4,1))
# for(j in 1:n.vit){
#     # plot(d$julian[d$lowtag==individs[j]],d$logstep[d$lowtag==individs[j]],main=d.vit$lowtag[j])
#     # points(final.out2$alarm[[j]]$julian,final.out2$alarm[[j]]$logstep,col=2)
#     # abline(v=d.vit$juliandropped[j])
#     plot(d$julian[d$lowtag==individs[j]],d$scalealt[d$lowtag==individs[j]],main=d.vit$lowtag[j])
#     points(final.out2$alarm[[j]]$julian,final.out2$alarm[[j]]$scalealt,col=2)
#     abline(v=d.vit$juliandropped[j])
#     
#     plot(d$julian[d$lowtag==individs[j]],d$scaletemp[d$lowtag==individs[j]],main=d.vit$lowtag[j])
#     points(final.out2$alarm[[j]]$julian,final.out2$alarm[[j]]$scaletemp,col=2)
#     abline(v=d.vit$juliandropped[j])
#     
#     plot(d$julian[d$lowtag==individs[j]],d$scaledx[d$lowtag==individs[j]],main=d.vit$lowtag[j])
#     points(final.out2$alarm[[j]]$julian,final.out2$alarm[[j]]$scaledx,col=2)
#     abline(v=d.vit$juliandropped[j])
#     
#     plot(d$julian[d$lowtag==individs[j]],d$logscalestep[d$lowtag==individs[j]],main=d.vit$lowtag[j])
#     points(final.out2$alarm[[j]]$julian,final.out2$alarm[[j]]$logscalestep,col=2)
#     abline(v=d.vit$juliandropped[j])
#     
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
# #width of the window from above
# #relative angle
# vit.drop[i,j]=d.tmp$dropped[i]
# julian.out[i,j]=d.tmp$julian[i]
# 
# #min/max features
# min.distance[i,j]=min(d.tmp$distance[(i-w):i])
# max.angle[i,j]=max(d.tmp$rel.angle[(i-w):i])
# 
# #mean features
# mu.distance[i,j]=mean(d.tmp$distance[(i-w):i])
# mu.turn.angle[i,j]=mean(d.tmp$rel.angle[(i-w):i])
# mu.bearing[i,j]=mean(d.tmp$bearing[(i-w):i])
# mu.r2n[i,j]=mean(d.tmp$R2n[(i-w):i])
# mu.altitude[i,j]=mean(d.tmp$altitude[(i-w):i])
# mu.temp[i,j]=mean(d.tmp$temp[(i-w):i])
# mu.dx[i,j]=mean(d.tmp$dx[(i-w):i])
# mu.dy[i,j]=mean(d.tmp$dy[(i-w):i])
# 
# #variation features
# sig.distance[i,j]=sd(d.tmp$distance[(i-w):i])
# sig.turn.angle[i,j]=sd(d.tmp$rel.angle[(i-w):i])
# sig.bearing[i,j]=sd(d.tmp$bearing[(i-w):i])
# sig.r2n[i,j]=sd(d.tmp$R2n[(i-w):i])
# sig.altitude[i,j]=sd(d.tmp$altitude[(i-w):i])
# sig.temp[i,j]=sd(d.tmp$temp[(i-w):i])
# sig.dx[i,j]=sd(d.tmp$dx[(i-w):i])
# sig.dy[i,j]=sd(d.tmp$dy[(i-w):i])
# 
# #correlation features
# corr.dist.angle[i,j]=corr(cbind(d.tmp$distance[(i-w):i],d.tmp$rel.angle[(i-w):i]))
# corr.dist.bear[i,j]=corr(cbind(d.tmp$distance[(i-w):i],d.tmp$bearing[(i-w):i]))
# corr.dist.alt[i,j]=corr(cbind(d.tmp$distance[(i-w):i],d.tmp$altitude[(i-w):i]))
# corr.angle.bear[i,j]=corr(cbind(d.tmp$rel.angle[(i-w):i],d.tmp$bearing[(i-w):i]))
# corr.angle.alt[i,j]=corr(cbind(d.tmp$rel.angle[(i-w):i],d.tmp$altitude[(i-w):i]))
# corr.bear.alt[i,j]=corr(cbind(d.tmp$bearing[(i-w):i],d.tmp$altitude[(i-w):i]))
n.avg = n.obs-w-1
dim(julian.out)
dim(mu.distance)

par(mfrow=c(2,1))
acf(mu.distance[1:n.avg[1],1])
plot(mu.distance[1:n.avg[1],1])
plot(julian.out[1:n.avg[1],1],mu.distance[1:n.avg[1],1])
abline(v=d.vit$juliandropped[1],col=2)


#make a population level data frame to fit random forest
df.ma = data.frame(cbind(c(vit.drop),c(julian.out),
                         c(min.distance),
                         c(max.angle),
                         c(mu.distance),
                         c(mu.turn.angle),
                         c(mu.altitude),
                         c(mu.temp),
                         c(mu.dx),
                         c(mu.dy),
                         c(sig.distance),
                         c(sig.turn.angle),
                         c(sig.altitude),
                         c(sig.temp),
                         c(sig.dx),
                         c(sig.dy)),stringsAsFactors = FALSE)
remove =which(is.na(df.ma),arr.ind = T)[,1]
df.ma = df.ma[-remove,]

names(df.ma) = c("vitdrop","julian.out","min.distance","max.angle","mu.distance","mu.turn.angle","mu.altitude","mu.temp",
                 "mu.dx","mu.dy","sig.distance","sig.turn.angle","sig.altitude","sig.temp","sig.dx","sig.dy")
for(h in 2:dim(df.ma)[2]){
    class(df.ma[,h])="numeric"
}
df.ma$vitdrop = as.factor(df.ma$vitdrop)

fit.rf <- randomForest(vitdrop ~ max.angle + mu.distance + mu.turn.angle + mu.altitude + mu.temp + mu.dx +
                        mu.dy + sig.distance + sig.turn.angle + sig.altitude + sig.temp + sig.dx +sig.dy,
                    data=df.ma, 
                    localImp=TRUE, 
                    ntree=501)
fit.rf
summary(fit.rf)
plot(fit.rf)
varImpPlot(fit.rf)
length(mu.distance[1,])
class(mu.distance)
n.obs

min_depth_frame = min_depth_distribution(fit.rf)
plot_min_depth_distribution(min_depth_frame)
importance_frame = measure_importance(fit.rf)
importance_frame
plot_multi_way_importance(importance_frame, size_measure = "no_of_nodes")

plot_multi_way_importance(importance_frame, size_measure = "p_value")

plot_importance_ggpairs(importance_frame)
plot_importance_rankings(importance_frame)
vars = important_variables(importance_frame, k = 5, measures = c("mean_min_depth", "no_of_trees"))
interactions_frame <- min_depth_interactions(fit.rf, vars)
plot_min_depth_interactions(interactions_frame)
interactions_frame <- min_depth_interactions(fit.rf, vars, mean_sample = "relevant_trees", uncond_mean_sample = "relevant_trees")
plot_min_depth_interactions(interactions_frame)
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE, echo = FALSE)

explain_forest(fit.rf, interactions = TRUE, data = df.ma)

mu.dist.acf=apply(mu.distance[1,1:n.avg[1]],acf,na.rm=TRUE)


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


fitit
out.nb    
naiveBayes(train,as.factor(vit.drop[w,-1]))

head(pred.nb)


#####################################################################################################
###
### Moving window analysis svm
###
###
#####################################################################################################
fit.svm=svm(vitdrop~max.angle + mu.distance + mu.turn.angle + mu.altitude + mu.temp + mu.dx +
        mu.dy + sig.distance + sig.turn.angle + sig.altitude + sig.temp + sig.dx +sig.dy,data=df.ma)
pred.svm <- predict(fit.svm, df.ma)
points(df.ma$mu.distance, pred.svm, col = "blue", pch=4)


corr.bear.alt=matrix(NA,ncol=n.vit,nrow=max(n.obs))
dim(mu.distance)[1]
svm.window=rep(list(),dim(mu.distance)[1])
predict.window=rep(list(),dim(mu.distance)[1])
glm.window=rep(list(),dim(mu.distance)[1])

w=6
for(j in 1:n.vit){
        data.temp =data.frame(cbind(vit.drop[1:n.avg[j],j],mu.distance[1:n.avg[j],j],sig.distance[1:n.avg[j],j],mu.dx[1:n.avg[j],j]),stringsAsFactors = FALSE)
        names(data.temp)=c("vitdrop","mu.distance","sig.distance","mu.dx")
        class(data.temp$mu.distance)="numeric"
        class(data.temp$sig.distance)="numeric"
        class(data.temp$mu.dx)="numeric"
        data.temp$vitdrop=as.factor(data.temp$vitdrop)
        data.temp$mu.distance=scale(data.temp$mu.distance)
        data.temp$mu.distance.sq=(data.temp$mu.distance)^2
        svm.window[[j]] = svm(vitdrop~ mu.distance+sig.distance+ mu.distance.sq,data=data.temp,probability=TRUE)
        predict.window[[j]] = predict(svm.window[[j]],data.temp)
        glm.window[[j]]= glm(vitdrop~ mu.distance+sig.distance+mu.distance.sq,data=data.temp,family="binomial")
}
lapply(glm.window,summary)
glm.predict=lapply(glm.window,predict)
plot(glm.predict[[1]])

par(mfrow=c(6,1))
for(j in 1:n.vit){
    plot(julian.out[1:n.avg[j],j],predict.window[[j]])
}



dev.off()
lapply(predict.window,sum)

class(predict.window[[1]])
names(d)[32:39]
try=glm(d$dropped~scalealt   +  scaletemp +  scalebear   + scaledx  +  scaleR2n+logscalestep+cover,data=d,family=binomial(link="logit"),na.action(na.pass))
summary(try)
step(try) 

glm.loo=rep(list(),n.vit)
glm.loo.predict=rep(list(),n.vit)

for(j in 1:n.vit){
    d.loo=d[d$lowtag!=individs[j],]
    glm.loo[[j]]=glm(d.loo$dropped~scalealt+scaleR2n+logscalestep,data=d.loo,family=binomial(link="logit"),na.action(na.pass))
    d.new=d[d$lowtag==individs[j],]
    glm.loo.predict[[j]]=predict(glm.loo[[j]],d.new[,c(32,38,39)],type="response")
    }

par(mfrow=c(2,1))
for(j in 1:n.vit){
    plot(d$julian[d$lowtag==individs[j]],glm.loo.predict[[j]],main=d.vit$lowtag[j])
    abline(v=d.vit$juliandropped[j],col=2)
    plot(d$julian[d$lowtag==individs[j]],sigmoid(glm.loo.predict[[j]]),main=d.vit$lowtag[j])
    abline(v=d.vit$juliandropped[j],col=2)
    }

plot(d$julian[d$lowtag==individs[j]],d$dropped[d$lowtag==individs[j]])

plot()

for(j in 1:n.vit){
    length(glm.loo[[j]])
    length(glm.loo.predict[[j]])
    data.temp=d[d$lowtag==individs[j],]
    data.temp$julian
    }
    

for()
