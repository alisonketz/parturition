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
### Load previous runs
###

load("anomaly_v8_day.Rdata")

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
### Add sigma_day as a feature, for the within day variation of distances covered
###

sig.temp=rep(NA,dim(d)[1])
for(i in 1:dim(d)[1]){
    sig.temp[i]=sd(d$distance[d$lowtag==d$lowtag[i] & d$julian==d$julian[i]],na.rm=TRUE)
}
d$sigma.dist.day=sig.temp

###
### Add sigma_dx_day for within day variation in dx
###

sig.dx=rep(NA,dim(d)[1])
for(i in 1:dim(d)[1]){
    sig.dx[i]=sd(d$dx[d$lowtag==d$lowtag[i] & d$julian==d$julian[i]],na.rm=TRUE)
}
d$sigma.dx.day=sig.dx

###
### Add sigma_dy_day for within day variation in dy
###

sig.dy=rep(NA,dim(d)[1])
for(i in 1:dim(d)[1]){
    sig.dy[i]=sd(d$dy[d$lowtag==d$lowtag[i] & d$julian==d$julian[i]],na.rm=TRUE)
}
d$sigma.dy.day=sig.dy


###
### Add sigma_turn_day for within day variation in turn angle
###

sig.angle=rep(NA,dim(d)[1])
for(i in 1:dim(d)[1]){
    sig.angle[i]=sd(d$rel.angle[d$lowtag==d$lowtag[i] & d$julian==d$julian[i]],na.rm=TRUE)
}
d$sigma.angle.day=sig.angle

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
plot(d$sigma.dist.day,col=d$lowtag)
plot(d$sigma.dx.day,col=d$lowtag)
plot(d$sigma.dy.day,col=d$lowtag)



d$step = d$distance/d$timediff
# d$logstep = log(d$step)



##############################################################################################
###
### creating a day based time step data frame, and running anomaly detection 
###
##############################################################################################

# features to calculate within day means and standard deviations, and assign to daily
#lowtag, julian, dropped,step,altittude, temp, bearing, rel.angle, R2n, dx,dy
#lowtag, julian, dropped,
# mean and sd of : step,altitude, temp, bearing, rel.angle, R2n, dx,dy
#matrix dim should  = 8*2+3 = 19


d$lowtag = as.character(d$lowtag)
d.day.temp=matrix(NA,nr = 1,nc=19)#this will be too big, will trim after
for(j in 1:n.vit){
    d.temp=d[d$lowtag == individs[j],]
    julian.temp = unique(d.temp$julian)
    temp.mat = matrix(NA,nr = length(julian.temp),nc = 19)
    for(i in 1:length(julian.temp)){
        temp.mat[i,1] = d.temp$lowtag[1]
        temp.mat[i,2] = julian.temp[i]
        temp.mat[i,3] = ifelse(sum(d.temp$dropped[d.temp$julian == julian.temp[i]])>1,1,0)
        temp.mat[i,4] = mean(d.temp$step[d.temp$julian == julian.temp[i]],na.rm=TRUE)#mu.step
        temp.mat[i,5] = sd(d.temp$step[d.temp$julian == julian.temp[i]],na.rm=TRUE)#sig.step
        temp.mat[i,6] = mean(d.temp$altitude[d.temp$julian == julian.temp[i]],na.rm=TRUE)#mu.altitude
        temp.mat[i,7] = sd(d.temp$altitude[d.temp$julian == julian.temp[i]],na.rm=TRUE)#sig.altitude
        temp.mat[i,8] = mean(d.temp$temp[d.temp$julian == julian.temp[i]],na.rm=TRUE)#mu.temp
        temp.mat[i,9] = sd(d.temp$temp[d.temp$julian == julian.temp[i]],na.rm=TRUE)#sig.temp
        temp.mat[i,10] = mean(d.temp$bearing[d.temp$julian == julian.temp[i]],na.rm=TRUE)#mu.bearing
        temp.mat[i,11] = sd(d.temp$bearing[d.temp$julian == julian.temp[i]],na.rm=TRUE)#sig.bearing
        temp.mat[i,12] = mean(d.temp$rel.angle[d.temp$julian == julian.temp[i]],na.rm=TRUE)#mu.rel.angle
        temp.mat[i,13] = sd(d.temp$rel.angle[d.temp$julian == julian.temp[i]],na.rm=TRUE)#sig.rel.angle
        temp.mat[i,14] = mean(d.temp$R2n[d.temp$julian == julian.temp[i]],na.rm=TRUE)#mu.R2n
        temp.mat[i,15] = sd(d.temp$R2n[d.temp$julian == julian.temp[i]],na.rm=TRUE)#sig.R2n
        temp.mat[i,16] = mean(d.temp$dx[d.temp$julian == julian.temp[i]],na.rm=TRUE)#mu.dx
        temp.mat[i,17] = sd(d.temp$dx[d.temp$julian == julian.temp[i]],na.rm=TRUE)#sig.dx
        temp.mat[i,18] = mean(d.temp$dy[d.temp$julian == julian.temp[i]],na.rm=TRUE)#mu.dy
        temp.mat[i,19] = sd(d.temp$dy[d.temp$julian == julian.temp[i]],na.rm=TRUE)#sig.dy
    }
    d.day.temp = rbind(d.day.temp,temp.mat)
}
d.day.temp=d.day.temp[-1,]


d.tres.temp=matrix(NA,nr = 1,nc=19)#this will be too big, will trim after
for(j in 1:n.vit){
    d.temp=d[d$lowtag == individs[j],]
    julian.temp = unique(d.temp$julian)
    temp.mat = matrix(NA,nr = length(julian.temp),nc = 19)
    for(i in 1:length(julian.temp)){
        temp.mat[i,1] = d.temp$lowtag[1]
        temp.mat[i,2] = julian.temp[i]
        temp.mat[i,3] = ifelse(sum(d.temp$dropped[d.temp$julian == julian.temp[i]])>1,1,0)
        temp.mat[i,4] = mean(d.temp$step[d.temp$julian == julian.temp[i] | d.temp$julian == (julian.temp[i]-1) |d.temp$julian == (julian.temp[i]-2)],na.rm=TRUE)#mu.step
        temp.mat[i,5] = sd(d.temp$step[d.temp$julian == julian.temp[i] | d.temp$julian == (julian.temp[i]-1) |d.temp$julian == (julian.temp[i]-2)],na.rm=TRUE)#sig.step
        temp.mat[i,6] = mean(d.temp$altitude[d.temp$julian == julian.temp[i] | d.temp$julian == (julian.temp[i]-1) |d.temp$julian == (julian.temp[i]-2)],na.rm=TRUE)#mu.altitude
        temp.mat[i,7] = sd(d.temp$altitude[d.temp$julian == julian.temp[i] | d.temp$julian == (julian.temp[i]-1) |d.temp$julian == (julian.temp[i]-2)],na.rm=TRUE)#sig.altitude
        temp.mat[i,8] = mean(d.temp$temp[d.temp$julian == julian.temp[i] | d.temp$julian == (julian.temp[i]-1) |d.temp$julian == (julian.temp[i]-2)],na.rm=TRUE)#mu.temp
        temp.mat[i,9] = sd(d.temp$temp[d.temp$julian == julian.temp[i] | d.temp$julian == (julian.temp[i]-1) |d.temp$julian == (julian.temp[i]-2)],na.rm=TRUE)#sig.temp
        temp.mat[i,10] = mean(d.temp$bearing[d.temp$julian == julian.temp[i] | d.temp$julian == (julian.temp[i]-1) |d.temp$julian == (julian.temp[i]-2)],na.rm=TRUE)#mu.bearing
        temp.mat[i,11] = sd(d.temp$bearing[d.temp$julian == julian.temp[i] | d.temp$julian == (julian.temp[i]-1) |d.temp$julian == (julian.temp[i]-2)],na.rm=TRUE)#sig.bearing
        temp.mat[i,12] = mean(d.temp$rel.angle[d.temp$julian == julian.temp[i] | d.temp$julian == (julian.temp[i]-1) |d.temp$julian == (julian.temp[i]-2)],na.rm=TRUE)#mu.rel.angle
        temp.mat[i,13] = sd(d.temp$rel.angle[d.temp$julian == julian.temp[i] | d.temp$julian == (julian.temp[i]-1) |d.temp$julian == (julian.temp[i]-2)],na.rm=TRUE)#sig.rel.angle
        temp.mat[i,14] = mean(d.temp$R2n[d.temp$julian == julian.temp[i] | d.temp$julian == (julian.temp[i]-1) |d.temp$julian == (julian.temp[i]-2)],na.rm=TRUE)#mu.R2n
        temp.mat[i,15] = sd(d.temp$R2n[d.temp$julian == julian.temp[i] | d.temp$julian == (julian.temp[i]-1) |d.temp$julian == (julian.temp[i]-2)],na.rm=TRUE)#sig.R2n
        temp.mat[i,16] = mean(d.temp$dx[d.temp$julian == julian.temp[i] | d.temp$julian == (julian.temp[i]-1) |d.temp$julian == (julian.temp[i]-2)],na.rm=TRUE)#mu.dx
        temp.mat[i,17] = sd(d.temp$dx[d.temp$julian == julian.temp[i] | d.temp$julian == (julian.temp[i]-1) |d.temp$julian == (julian.temp[i]-2)],na.rm=TRUE)#sig.dx
        temp.mat[i,18] = mean(d.temp$dy[d.temp$julian == julian.temp[i] | d.temp$julian == (julian.temp[i]-1) |d.temp$julian == (julian.temp[i]-2)],na.rm=TRUE)#mu.dy
        temp.mat[i,19] = sd(d.temp$dy[d.temp$julian == julian.temp[i] | d.temp$julian == (julian.temp[i]-1) |d.temp$julian == (julian.temp[i]-2)],na.rm=TRUE)#sig.dy
    }
    d.tres.temp = rbind(d.tres.temp,temp.mat)
}
d.tres.temp=d.tres.temp[-1,]
d.tres.temp=d.tres.temp[,4:19]

head(d.tres.temp)

d.day.tres=data.frame(cbind(d.day.temp,d.tres.temp),stringsAsFactors = FALSE)

names(d.day.tres)=c("lowtag","julian","dropped","mu.step.day","sig.step.day","mu.altitude.day","sig.altitude.day","mu.temp.day","sig.temp.day",
                    "mu.bearing.day","sig.bearing.day","mu.rel.angle.day","sig.rel.angle.day","mu.R2n.day","sig.R2n.day","mu.dx.day","sig.dx.day","mu.dy.day","sig.dy.day","mu.step.tres","sig.step.tres","mu.altitude.tres","sig.altitude.tres","mu.temp.tres","sig.temp.tres",
                    "mu.bearing.tres","sig.bearing.tres","mu.rel.angle.tres","sig.rel.angle.tres","mu.R2n.tres","sig.R2n.tres","mu.dx.tres","sig.dx.tres","mu.dy.tres","sig.dy.tres")




dim(d.day.tres)


for(i in 3:35){
    d.day.tres[,i]=as.numeric(d.day.tres[,i])
}
head(d.day.tres)

###
### center and scale variables for prediction
###

for(j in 1:n.vit){
    d.day.tres[d.day.tres$lowtag == individs[j],4:35]=apply(d.day.tres[d.day.tres$lowtag == individs[j],4:35],2,scale)
}

d.day.tres$julian = as.integer(d.day.tres$julian)


###
### Paturition window
###

part.window=yday("2017-05-06")

###
### Total number of possible anomalies to detect
###

possible.hits=rep(NA,n.vit)
for(j in 1:n.vit){
    possible.hits[j] = sum(d.day.tres$dropped[d.day.tres$lowtag==individs[j]])
}
possible.hits


head()


##############################################################################################
###
### Expanding Anomaly Detection function to deal 
### with multiple individuals as well as single
### optimizing epsilon of features
###
##############################################################################################

d.day.tres$lowtag=as.factor(d.day.tres$lowtag)
pw=part.window
ph=possible.hits
nv=n.vit
ci=4:35#covariate indexes

source("anomalyDetectUniv.R")
source("anomalyDetect.R")
source("evaluate_ad.R")

out.test=anomalyDetect(n.vit=nv,part.window=pw,id=individs,d=d.day.tres,eps=rep(.01,length(ci)),covs.indx=ci)
out.test$detect.density
out.test$threshold.density
out.test$alarm

eval.temp=evaluate.time.pop(alarm=out.test$alarm,possible.hits=ph,n.vit=nv,vitdropday = d.vit$juliandropped)
eval.temp
# out.ind.5732=anomalyDetect(1,part.window=pw,id=c(5732),d[d$lowtag==5732,],eps=.12,33:41)
# out.ind.5004=anomalyDetect(1,part.window=126,id=c(5004),d[d$lowtag==5004,],eps=.12,33:41)
# evaluate.time.pop(alarm=out.ind.5732$alarm,possible.hits=possible.hits[5],n.vit=1,vitdropday = d.vit$juliandropped[5])
# evaluate.time.pop(alarm=out.ind.5004$alarm,possible.hits=possible.hits[1],n.vit=1,vitdropday = d.vit$juliandropped[1])


epsnum=15 #number of epsilon quantiles tested
covsnum=length(ci)
loo.eval=rep(list(),3*n.vit)# there are 3 criteria we can use to tune epsilon for anomaly detection
n.loo=n.vit-1 #number in leave one out training
decide=matrix(NA,nr=covsnum,ncol=3)
starting=Sys.time()
for(m in 1:3){
    for(j in 1:n.vit){
        df=d.day.tres[d.day.tres$lowtag!=individs[j],]
        train.cov=array(NA,c(covsnum,epsnum,3))
        for(h in 1:covsnum){
            for(k in 1:epsnum){# loop over epsilon values
                fit.train = anomalyDetect(n.vit=n.loo,part.window=pw,id=individs[-j],d=df,eps=k/100,covs.indx=ci[h])# fit model
                eval.temp=evaluate.time.pop(alarm=fit.train$alarm,possible.hits=ph[-j],n.vit=n.loo,vitdropday = d.vit$juliandropped[-j])
                train.cov[h,k,1]=eval.temp$out.prec
                train.cov[h,k,2]=eval.temp$out.recall
                train.cov[h,k,3]=eval.temp$out.F1# calculate eval = recall
            }
        }
        compare=apply(train.cov,c(1,3),max,na.rm=TRUE)
        for(h in 1:covsnum){
            for(i in 1:3){
                decide[h,i]=min(which(train.cov[h,,i]==compare[h,i]))
            }
        }
        decide=decide/100
        d.ind=d.day.tres[d.day.tres$lowtag==individs[j],]
        fit.test= anomalyDetect(n.vit=1,part.window=pw,id=individs[j],d=d.ind,eps=decide[,m],covs.indx=ci)# fit model
        eval.test=evaluate.time.pop(alarm=fit.test$alarm,possible.hits=ph[j],n.vit=1,vitdropday = d.vit$juliandropped[j])
        loo.indx=j+(m-1)*12
        loo.eval[[loo.indx]]=eval.test
    }    
}
ending=Sys.time()
total.time=ending-starting
total.time

for(j in 1:5){
    beep(sound=3)
}

 
check.prec=rep(NA,3*n.vit)
check.recall=rep(NA,3*n.vit)
check.F1=rep(NA,3*n.vit)
check.TP=rep(NA,3*n.vit)
check.FP=rep(NA,3*n.vit)
check.FN=rep(NA,3*n.vit)


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
# pdf("Anomaly_eval_1.pdf")
plot(test.ind,check.prec,col=test.crit,main="Precision",ylab="Precision",xlab="Individual",ylim=c(0,1))
legend(x=1.5,y=.9,fill=c(1:3),legend=c("Precision","Recall","F1"))
plot(test.ind,check.recall,col=test.crit,main="Recall",ylab="Recall",xlab="Individual")
plot(test.ind,check.F1,col=test.crit,main="F1",ylab="F1",xlab="Individual")
# dev.off()
# pdf("Anomaly_eval_2.pdf")
plot(test.ind,check.TP,col=test.crit,main="TP",ylab="TP",xlab="Individual")
plot(test.ind,check.FP,col=test.crit,main="FP",ylab="FP",xlab="Individual")
plot(test.ind,check.FN,col=test.crit,main="FN",ylab="FN",xlab="Individual")
legend("topright",fill=c(1:3),legend=c("Precision","Recall","F1"))
# dev.off()


par(mfrow=c(3,1))
# pdf("checks/Precision_check.pdf")
plot(check.prec,main="Precision when different criteria")
abline(v=c(12.5,24.5))
text(5,.25,"Precision")
text(18,.25,"Recall")
text(28,.25,"F1")
# dev.off()

# pdf("checks/Recall_check.pdf")
plot(check.recall,main="Recall when different criteria")
abline(v=c(12.5,24.5))
text(5,.8,"Precision")
text(18,.8,"Recall")
text(28,.8,"F1")
# dev.off()

# pdf("checks/F1_check.pdf")
plot(check.F1,main="F1 when different criteria")
abline(v=c(12.5,24.5))
text(5,.25,"Precision")
text(18,.25,"Recall")
text(28,.25,"F1")
# dev.off()

# pdf("checks/TP_check.pdf")
plot(check.TP,main="TP when diff criteria")
abline(v=c(12.5,24.5))
text(5,10,"Precision")
text(18,10,"Recall")
text(28,10,"F1")
# dev.off()

# pdf("checks/FP_check.pdf")
plot(check.FP,main="FP")
abline(v=c(12.5,24.5))
text(5,120,"Precision")
text(18,120,"Recall")
text(28,120,"F1")
# dev.off()

# pdf("checks/FN_check.pdf")
plot(check.FN,main="FN")
abline(v=c(12.5,24.5))
text(5,11,"Precision")
text(18,11,"Recall")
text(28,11,"F1")
# dev.off()


# ##############################################################################################
# ###
# ### Multivariate Anomaly Detection Using F1 criteria with above epsilon
# ### changing for all featurs predictors
# ###
# ##############################################################################################
names(d.day.tres)
#only plotting 4/16 covariates

###precision final model
fit.prec = anomalyDetect(n.vit=nv,part.window=pw,id=individs,d=d.day.tres,eps=decide[,1],covs.indx=ci)# fit model
eval.prec=evaluate.time.pop(alarm=fit.prec$alarm,possible.hits=ph,n.vit=nv,vitdropday = d.vit$juliandropped)
par(mfrow=c(4,1))
for(j in 1:n.vit){
    pdf(paste("outputplots_anomaly10/",individs[j],"_prec_part1.pdf",sep=""),width=6,height=9)
    par(mfrow=c(4,1))
    plot(d.day.tres$julian[d.day.tres$lowtag==individs[j]],d.day.tres$mu.step[d.day.tres$lowtag==individs[j]],main=paste(individs[j],"mu.step"))
    points(fit.prec$alarm[[j]]$julian,fit.prec$alarm[[j]]$mu.step,col=2)
    abline(v=d.vit$juliandropped[j])
    plot(d.day.tres$julian[d.day.tres$lowtag==individs[j]],d.day.tres$sig.step[d.day.tres$lowtag==individs[j]],main="sig.step")
    points(fit.prec$alarm[[j]]$julian,fit.prec$alarm[[j]]$sig.step,col=2)
    abline(v=d.vit$juliandropped[j])
    plot(d.day.tres$julian[d.day.tres$lowtag==individs[j]],d.day.tres$mu.altitude[d.day.tres$lowtag==individs[j]],main="mu.altitude")
    points(fit.prec$alarm[[j]]$julian,fit.prec$alarm[[j]]$mu.altitude,col=2)
    abline(v=d.vit$juliandropped[j])
    plot(d.day.tres$julian[d.day.tres$lowtag==individs[j]],d.day.tres$sig.altitude[d.day.tres$lowtag==individs[j]],main="sig.altitude")
    points(fit.prec$alarm[[j]]$julian,fit.prec$alarm[[j]]$sig.altitude,col=2)
    abline(v=d.vit$juliandropped[j])
    dev.off()
}
for(j in 1:n.vit){
    pdf(paste("outputplots_anomaly10/",individs[j],"_prec_part2.pdf",sep=""),width=6,height=9)
    par(mfrow=c(4,1))
    plot(d.day.tres$julian[d.day.tres$lowtag==individs[j]],d.day.tres$mu.temp[d.day.tres$lowtag==individs[j]],main=paste(individs[j],"mu.temp"))
    points(fit.prec$alarm[[j]]$julian,fit.prec$alarm[[j]]$mu.temp,col=2)
    abline(v=d.vit$juliandropped[j])
    plot(d.day.tres$julian[d.day.tres$lowtag==individs[j]],d.day.tres$sig.temp[d.day.tres$lowtag==individs[j]],main="sig.temp")
    points(fit.prec$alarm[[j]]$julian,fit.prec$alarm[[j]]$sig.temp,col=2)
    abline(v=d.vit$juliandropped[j])
    plot(d.day.tres$julian[d.day.tres$lowtag==individs[j]],d.day.tres$mu.bearing[d.day.tres$lowtag==individs[j]],main="mu.bearing")
    points(fit.prec$alarm[[j]]$julian,fit.prec$alarm[[j]]$mu.bearing,col=2)
    abline(v=d.vit$juliandropped[j])
    plot(d.day.tres$julian[d.day.tres$lowtag==individs[j]],d.day.tres$sig.bearing[d.day.tres$lowtag==individs[j]],main="sig.bearing")
    points(fit.prec$alarm[[j]]$julian,fit.prec$alarm[[j]]$sig.bearing,col=2)
    abline(v=d.vit$juliandropped[j])
    dev.off()
}

for(j in 1:n.vit){
    pdf(paste("outputplots_anomaly10/",individs[j],"_prec_part3.pdf",sep=""),width=6,height=9)
    par(mfrow=c(4,1))
    plot(d.day.tres$julian[d.day.tres$lowtag==individs[j]],d.day.tres$mu.rel.angle[d.day.tres$lowtag==individs[j]],main=paste(individs[j],"mu.rel.angle"))
    points(fit.prec$alarm[[j]]$julian,fit.prec$alarm[[j]]$mu.rel.angle,col=2)
    abline(v=d.vit$juliandropped[j])
    plot(d.day.tres$julian[d.day.tres$lowtag==individs[j]],d.day.tres$sig.rel.angle[d.day.tres$lowtag==individs[j]],main="sig.rel.angle")
    points(fit.prec$alarm[[j]]$julian,fit.prec$alarm[[j]]$sig.rel.angle,col=2)
    abline(v=d.vit$juliandropped[j])
    plot(d.day.tres$julian[d.day.tres$lowtag==individs[j]],d.day.tres$mu.R2n[d.day.tres$lowtag==individs[j]],main="mu.R2n")
    points(fit.prec$alarm[[j]]$julian,fit.prec$alarm[[j]]$mu.R2n,col=2)
    abline(v=d.vit$juliandropped[j])
    plot(d.day.tres$julian[d.day.tres$lowtag==individs[j]],d.day.tres$sig.R2n[d.day.tres$lowtag==individs[j]],main="sig.R2n")
    points(fit.prec$alarm[[j]]$julian,fit.prec$alarm[[j]]$sig.R2n,col=2)
    abline(v=d.vit$juliandropped[j])
    dev.off()
}
for(j in 1:n.vit){
    pdf(paste("outputplots_anomaly10/",individs[j],"_prec_part4.pdf",sep=""),width=6,height=9)
    par(mfrow=c(4,1))
    plot(d.day.tres$julian[d.day.tres$lowtag==individs[j]],d.day.tres$mu.dx[d.day.tres$lowtag==individs[j]],main=paste(individs[j],"mu.dx"))
    points(fit.prec$alarm[[j]]$julian,fit.prec$alarm[[j]]$mu.dx,col=2)
    abline(v=d.vit$juliandropped[j])
    plot(d.day.tres$julian[d.day.tres$lowtag==individs[j]],d.day.tres$sig.dx[d.day.tres$lowtag==individs[j]],main="sig.dx")
    points(fit.prec$alarm[[j]]$julian,fit.prec$alarm[[j]]$sig.dx,col=2)
    abline(v=d.vit$juliandropped[j])
    plot(d.day.tres$julian[d.day.tres$lowtag==individs[j]],d.day.tres$mu.dy[d.day.tres$lowtag==individs[j]],main="mu.dy")
    points(fit.prec$alarm[[j]]$julian,fit.prec$alarm[[j]]$mu.dy,col=2)
    abline(v=d.vit$juliandropped[j])
    plot(d.day.tres$julian[d.day.tres$lowtag==individs[j]],d.day.tres$sig.dy[d.day.tres$lowtag==individs[j]],main="sig.dy")
    points(fit.prec$alarm[[j]]$julian,fit.prec$alarm[[j]]$sig.dy,col=2)
    abline(v=d.vit$juliandropped[j])
    dev.off()
}



###recall final model
fit.recall = anomalyDetect(n.vit=nv,part.window=pw,id=individs,d=d.day.tres,eps=decide[,2],covs.indx=ci)# fit model
eval.recall=evaluate.time.pop(alarm=fit.recall$alarm,possible.hits=ph,n.vit=nv,vitdropday = d.vit$juliandropped)
par(mfrow=c(4,1))
for(j in 1:n.vit){
    pdf(paste("outputplots_anomaly10/",individs[j],"_recall_part1.pdf",sep=""),width=6,height=9)
    par(mfrow=c(4,1))
    plot(d.day.tres$julian[d.day.tres$lowtag==individs[j]],d.day.tres$mu.step[d.day.tres$lowtag==individs[j]],main=paste(individs[j],"mu.step"))
    points(fit.recall$alarm[[j]]$julian,fit.recall$alarm[[j]]$mu.step,col=2)
    abline(v=d.vit$juliandropped[j])
    plot(d.day.tres$julian[d.day.tres$lowtag==individs[j]],d.day.tres$sig.step[d.day.tres$lowtag==individs[j]],main="sig.step")
    points(fit.recall$alarm[[j]]$julian,fit.recall$alarm[[j]]$sig.step,col=2)
    abline(v=d.vit$juliandropped[j])
    plot(d.day.tres$julian[d.day.tres$lowtag==individs[j]],d.day.tres$mu.altitude[d.day.tres$lowtag==individs[j]],main="mu.altitude")
    points(fit.recall$alarm[[j]]$julian,fit.recall$alarm[[j]]$mu.altitude,col=2)
    abline(v=d.vit$juliandropped[j])
    plot(d.day.tres$julian[d.day.tres$lowtag==individs[j]],d.day.tres$sig.altitude[d.day.tres$lowtag==individs[j]],main="sig.altitude")
    points(fit.recall$alarm[[j]]$julian,fit.recall$alarm[[j]]$sig.altitude,col=2)
    abline(v=d.vit$juliandropped[j])
    dev.off()
}
for(j in 1:n.vit){
    pdf(paste("outputplots_anomaly10/",individs[j],"_recall_part2.pdf",sep=""),width=6,height=9)
    par(mfrow=c(4,1))
    plot(d.day.tres$julian[d.day.tres$lowtag==individs[j]],d.day.tres$mu.temp[d.day.tres$lowtag==individs[j]],main=paste(individs[j],"mu.temp"))
    points(fit.recall$alarm[[j]]$julian,fit.recall$alarm[[j]]$mu.temp,col=2)
    abline(v=d.vit$juliandropped[j])
    plot(d.day.tres$julian[d.day.tres$lowtag==individs[j]],d.day.tres$sig.temp[d.day.tres$lowtag==individs[j]],main="sig.temp")
    points(fit.recall$alarm[[j]]$julian,fit.recall$alarm[[j]]$sig.temp,col=2)
    abline(v=d.vit$juliandropped[j])
    plot(d.day.tres$julian[d.day.tres$lowtag==individs[j]],d.day.tres$mu.bearing[d.day.tres$lowtag==individs[j]],main="mu.bearing")
    points(fit.recall$alarm[[j]]$julian,fit.recall$alarm[[j]]$mu.bearing,col=2)
    abline(v=d.vit$juliandropped[j])
    plot(d.day.tres$julian[d.day.tres$lowtag==individs[j]],d.day.tres$sig.bearing[d.day.tres$lowtag==individs[j]],main="sig.bearing")
    points(fit.recall$alarm[[j]]$julian,fit.recall$alarm[[j]]$sig.bearing,col=2)
    abline(v=d.vit$juliandropped[j])
    dev.off()
}

for(j in 1:n.vit){
    pdf(paste("outputplots_anomaly10/",individs[j],"_recall_part3.pdf",sep=""),width=6,height=9)
    par(mfrow=c(4,1))
    plot(d.day.tres$julian[d.day.tres$lowtag==individs[j]],d.day.tres$mu.rel.angle[d.day.tres$lowtag==individs[j]],main=paste(individs[j],"mu.rel.angle"))
    points(fit.recall$alarm[[j]]$julian,fit.recall$alarm[[j]]$mu.rel.angle,col=2)
    abline(v=d.vit$juliandropped[j])
    plot(d.day.tres$julian[d.day.tres$lowtag==individs[j]],d.day.tres$sig.rel.angle[d.day.tres$lowtag==individs[j]],main="sig.rel.angle")
    points(fit.recall$alarm[[j]]$julian,fit.recall$alarm[[j]]$sig.rel.angle,col=2)
    abline(v=d.vit$juliandropped[j])
    plot(d.day.tres$julian[d.day.tres$lowtag==individs[j]],d.day.tres$mu.R2n[d.day.tres$lowtag==individs[j]],main="mu.R2n")
    points(fit.recall$alarm[[j]]$julian,fit.recall$alarm[[j]]$mu.R2n,col=2)
    abline(v=d.vit$juliandropped[j])
    plot(d.day.tres$julian[d.day.tres$lowtag==individs[j]],d.day.tres$sig.R2n[d.day.tres$lowtag==individs[j]],main="sig.R2n")
    points(fit.recall$alarm[[j]]$julian,fit.recall$alarm[[j]]$sig.R2n,col=2)
    abline(v=d.vit$juliandropped[j])
    dev.off()
}
for(j in 1:n.vit){
    pdf(paste("outputplots_anomaly10/",individs[j],"_recall_part4.pdf",sep=""),width=6,height=9)
    par(mfrow=c(4,1))
    plot(d.day.tres$julian[d.day.tres$lowtag==individs[j]],d.day.tres$mu.dx[d.day.tres$lowtag==individs[j]],main=paste(individs[j],"mu.dx"))
    points(fit.recall$alarm[[j]]$julian,fit.recall$alarm[[j]]$mu.dx,col=2)
    abline(v=d.vit$juliandropped[j])
    plot(d.day.tres$julian[d.day.tres$lowtag==individs[j]],d.day.tres$sig.dx[d.day.tres$lowtag==individs[j]],main="sig.dx")
    points(fit.recall$alarm[[j]]$julian,fit.recall$alarm[[j]]$sig.dx,col=2)
    abline(v=d.vit$juliandropped[j])
    plot(d.day.tres$julian[d.day.tres$lowtag==individs[j]],d.day.tres$mu.dy[d.day.tres$lowtag==individs[j]],main="mu.dy")
    points(fit.recall$alarm[[j]]$julian,fit.recall$alarm[[j]]$mu.dy,col=2)
    abline(v=d.vit$juliandropped[j])
    plot(d.day.tres$julian[d.day.tres$lowtag==individs[j]],d.day.tres$sig.dy[d.day.tres$lowtag==individs[j]],main="sig.dy")
    points(fit.recall$alarm[[j]]$julian,fit.recall$alarm[[j]]$sig.dy,col=2)
    abline(v=d.vit$juliandropped[j])
    dev.off()
}



### F1 final model
fit.F1 = anomalyDetect(n.vit=nv,part.window=pw,id=individs,d=d.day.tres,eps=decide[,3],covs.indx=ci)# fit model
eval.F1=evaluate.time.pop(alarm=fit.F1$alarm,possible.hits=ph,n.vit=nv,vitdropday = d.vit$juliandropped)
for(j in 1:n.vit){
    pdf(paste("outputplots_anomaly10/",individs[j],"_recall_part1.pdf",sep=""),width=6,height=9)
    par(mfrow=c(4,1))
    plot(d.day.tres$julian[d.day.tres$lowtag==individs[j]],d.day.tres$mu.step[d.day.tres$lowtag==individs[j]],main=paste(individs[j],"mu.step"))
    points(fit.recall$alarm[[j]]$julian,fit.recall$alarm[[j]]$mu.step,col=2)
    abline(v=d.vit$juliandropped[j])
    plot(d.day.tres$julian[d.day.tres$lowtag==individs[j]],d.day.tres$sig.step[d.day.tres$lowtag==individs[j]],main="sig.step")
    points(fit.recall$alarm[[j]]$julian,fit.recall$alarm[[j]]$sig.step,col=2)
    abline(v=d.vit$juliandropped[j])
    plot(d.day.tres$julian[d.day.tres$lowtag==individs[j]],d.day.tres$mu.altitude[d.day.tres$lowtag==individs[j]],main="mu.altitude")
    points(fit.recall$alarm[[j]]$julian,fit.recall$alarm[[j]]$mu.altitude,col=2)
    abline(v=d.vit$juliandropped[j])
    plot(d.day.tres$julian[d.day.tres$lowtag==individs[j]],d.day.tres$sig.altitude[d.day.tres$lowtag==individs[j]],main="sig.altitude")
    points(fit.recall$alarm[[j]]$julian,fit.recall$alarm[[j]]$sig.altitude,col=2)
    abline(v=d.vit$juliandropped[j])
    dev.off()
}
for(j in 1:n.vit){
    pdf(paste("outputplots_anomaly10/",individs[j],"_recall_part2.pdf",sep=""),width=6,height=9)
    par(mfrow=c(4,1))
    plot(d.day.tres$julian[d.day.tres$lowtag==individs[j]],d.day.tres$mu.temp[d.day.tres$lowtag==individs[j]],main=paste(individs[j],"mu.temp"))
    points(fit.recall$alarm[[j]]$julian,fit.recall$alarm[[j]]$mu.temp,col=2)
    abline(v=d.vit$juliandropped[j])
    plot(d.day.tres$julian[d.day.tres$lowtag==individs[j]],d.day.tres$sig.temp[d.day.tres$lowtag==individs[j]],main="sig.temp")
    points(fit.recall$alarm[[j]]$julian,fit.recall$alarm[[j]]$sig.temp,col=2)
    abline(v=d.vit$juliandropped[j])
    plot(d.day.tres$julian[d.day.tres$lowtag==individs[j]],d.day.tres$mu.bearing[d.day.tres$lowtag==individs[j]],main="mu.bearing")
    points(fit.recall$alarm[[j]]$julian,fit.recall$alarm[[j]]$mu.bearing,col=2)
    abline(v=d.vit$juliandropped[j])
    plot(d.day.tres$julian[d.day.tres$lowtag==individs[j]],d.day.tres$sig.bearing[d.day.tres$lowtag==individs[j]],main="sig.bearing")
    points(fit.recall$alarm[[j]]$julian,fit.recall$alarm[[j]]$sig.bearing,col=2)
    abline(v=d.vit$juliandropped[j])
    dev.off()
}

for(j in 1:n.vit){
    pdf(paste("outputplots_anomaly10/",individs[j],"_recall_part3.pdf",sep=""),width=6,height=9)
    par(mfrow=c(4,1))
    plot(d.day.tres$julian[d.day.tres$lowtag==individs[j]],d.day.tres$mu.rel.angle[d.day.tres$lowtag==individs[j]],main=paste(individs[j],"mu.rel.angle"))
    points(fit.recall$alarm[[j]]$julian,fit.recall$alarm[[j]]$mu.rel.angle,col=2)
    abline(v=d.vit$juliandropped[j])
    plot(d.day.tres$julian[d.day.tres$lowtag==individs[j]],d.day.tres$sig.rel.angle[d.day.tres$lowtag==individs[j]],main="sig.rel.angle")
    points(fit.recall$alarm[[j]]$julian,fit.recall$alarm[[j]]$sig.rel.angle,col=2)
    abline(v=d.vit$juliandropped[j])
    plot(d.day.tres$julian[d.day.tres$lowtag==individs[j]],d.day.tres$mu.R2n[d.day.tres$lowtag==individs[j]],main="mu.R2n")
    points(fit.recall$alarm[[j]]$julian,fit.recall$alarm[[j]]$mu.R2n,col=2)
    abline(v=d.vit$juliandropped[j])
    plot(d.day.tres$julian[d.day.tres$lowtag==individs[j]],d.day.tres$sig.R2n[d.day.tres$lowtag==individs[j]],main="sig.R2n")
    points(fit.recall$alarm[[j]]$julian,fit.recall$alarm[[j]]$sig.R2n,col=2)
    abline(v=d.vit$juliandropped[j])
    dev.off()
}
for(j in 1:n.vit){
    pdf(paste("outputplots_anomaly10/",individs[j],"_recall_part4.pdf",sep=""),width=6,height=9)
    par(mfrow=c(4,1))
    plot(d.day.tres$julian[d.day.tres$lowtag==individs[j]],d.day.tres$mu.dx[d.day.tres$lowtag==individs[j]],main=paste(individs[j],"mu.dx"))
    points(fit.recall$alarm[[j]]$julian,fit.recall$alarm[[j]]$mu.dx,col=2)
    abline(v=d.vit$juliandropped[j])
    plot(d.day.tres$julian[d.day.tres$lowtag==individs[j]],d.day.tres$sig.dx[d.day.tres$lowtag==individs[j]],main="sig.dx")
    points(fit.recall$alarm[[j]]$julian,fit.recall$alarm[[j]]$sig.dx,col=2)
    abline(v=d.vit$juliandropped[j])
    plot(d.day.tres$julian[d.day.tres$lowtag==individs[j]],d.day.tres$mu.dy[d.day.tres$lowtag==individs[j]],main="mu.dy")
    points(fit.recall$alarm[[j]]$julian,fit.recall$alarm[[j]]$mu.dy,col=2)
    abline(v=d.vit$juliandropped[j])
    plot(d.day.tres$julian[d.day.tres$lowtag==individs[j]],d.day.tres$sig.dy[d.day.tres$lowtag==individs[j]],main="sig.dy")
    points(fit.recall$alarm[[j]]$julian,fit.recall$alarm[[j]]$sig.dy,col=2)
    abline(v=d.vit$juliandropped[j])
    dev.off()
}




cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

pdf(paste("outputplots_anomaly10/","AD_eval_F1_num1.pdf",sep=""),width=6,height=9)
par(mfrow=c(3,1))   
for(j in 1:3){
    tab = table(fit.F1$alarm[[j]]$julian)
    plot(tab,main=paste("Anomaly Detection",individs[j]),xlab = "Julian day",ylab = "Number of hits",xlim=c(126,188))
    abline(v=d.vit$juliandropped[j],col=cbPalette[4],lty=2)
    }
dev.off()

pdf(paste("outputplots_anomaly10/","AD_eval_F1_num2.pdf",sep=""),width=6,height=9)
par(mfrow=c(3,1))   
for(j in c(4,6)){
    tab = table(fit.F1$alarm[[j]]$julian)
    plot(tab,main=paste("Anomaly Detection",individs[j]),xlab = "Julian day",ylab = "Number of hits",xlim=c(126,188))
    abline(v=d.vit$juliandropped[j],col=cbPalette[4],lty=2)
}
dev.off()

pdf(paste("outputplots_anomaly10/","AD_eval_F1_num3.pdf",sep=""),width=6,height=9)
par(mfrow=c(3,1))   
for(j in 7:9){
    tab = table(fit.F1$alarm[[j]]$julian)
    plot(tab,main=paste("Anomaly Detection",individs[j]),xlab = "Julian day",ylab = "Number of hits",xlim=c(126,188))
    abline(v=d.vit$juliandropped[j],col=cbPalette[4],lty=2)
}
dev.off()

pdf(paste("outputplots_anomaly10/","AD_eval_F1_num4.pdf",sep=""),width=6,height=9)
par(mfrow=c(3,1))   
for(j in 10:12){
    tab = table(fit.F1$alarm[[j]]$julian)
    plot(tab,main=paste("Anomaly Detection",individs[j]),xlab = "Julian day",ylab = "Number of hits",xlim=c(126,188))
    abline(v=d.vit$juliandropped[j],col=cbPalette[4],lty=2)
}
dev.off()
 


pdf(paste("outputplots_anomaly10/","AD_eval_recall_num1.pdf",sep=""),width=6,height=9)
par(mfrow=c(3,1))   
for(j in 1:3){
    tab = table(fit.recall$alarm[[j]]$julian)
    plot(tab,main=paste("Anomaly Detection",individs[j]),xlab = "Julian day",ylab = "Number of hits",xlim=c(126,188))
    abline(v=d.vit$juliandropped[j],col=cbPalette[4],lty=2)
}
dev.off()

pdf(paste("outputplots_anomaly10/","AD_eval_recall_num2.pdf",sep=""),width=6,height=9)
par(mfrow=c(3,1))   
for(j in 4:6){
    tab = table(fit.recall$alarm[[j]]$julian)
    plot(tab,main=paste("Anomaly Detection",individs[j]),xlab = "Julian day",ylab = "Number of hits",xlim=c(126,188))
    abline(v=d.vit$juliandropped[j],col=cbPalette[4],lty=2)
}
dev.off()

pdf(paste("outputplots_anomaly10/","AD_eval_recall_num3.pdf",sep=""),width=6,height=9)
par(mfrow=c(3,1))   
for(j in 7:9){
    tab = table(fit.recall$alarm[[j]]$julian)
    plot(tab,main=paste("Anomaly Detection",individs[j]),xlab = "Julian day",ylab = "Number of hits",xlim=c(126,188))
    abline(v=d.vit$juliandropped[j],col=cbPalette[4],lty=2)
}
dev.off()

pdf(paste("outputplots_anomaly10/","AD_eval_recall_num4.pdf",sep=""),width=6,height=9)
par(mfrow=c(3,1))   
for(j in 10:12){
    tab = table(fit.recall$alarm[[j]]$julian)
    plot(tab,main=paste("Anomaly Detection",individs[j]),xlab = "Julian day",ylab = "Number of hits",xlim=c(126,188))
    abline(v=d.vit$juliandropped[j],col=cbPalette[4],lty=2)
}
dev.off()


pdf(paste("outputplots_anomaly10/","AD_eval_prec_num1.pdf",sep=""),width=6,height=9)
par(mfrow=c(3,1))   
for(j in 1:3){
    tab = table(fit.prec$alarm[[j]]$julian)
    plot(tab,main=paste("Anomaly Detection",individs[j]),xlab = "Julian day",ylab = "Number of hits",xlim=c(126,188))
    abline(v=d.vit$juliandropped[j],col=cbPalette[4],lty=2)
}
dev.off()

pdf(paste("outputplots_anomaly10/","AD_eval_prec_num2.pdf",sep=""),width=6,height=9)
par(mfrow=c(3,1))   
for(j in c(4,6)){
    tab = table(fit.prec$alarm[[j]]$julian)
    plot(tab,main=paste("Anomaly Detection",individs[j]),xlab = "Julian day",ylab = "Number of hits",xlim=c(126,188))
    abline(v=d.vit$juliandropped[j],col=cbPalette[4],lty=2)
}
dev.off()

pdf(paste("outputplots_anomaly10/","AD_eval_prec_num3.pdf",sep=""),width=6,height=9)
par(mfrow=c(3,1))   
for(j in 7:9){
    tab = table(fit.prec$alarm[[j]]$julian)
    plot(tab,main=paste("Anomaly Detection",individs[j]),xlab = "Julian day",ylab = "Number of hits",xlim=c(126,188))
    abline(v=d.vit$juliandropped[j],col=cbPalette[4],lty=2)
}
dev.off()

pdf(paste("outputplots_anomaly10/","AD_eval_prec_num4.pdf",sep=""),width=6,height=9)
par(mfrow=c(3,1))   
for(j in 10:12){
    tab = table(fit.prec$alarm[[j]]$julian)
    plot(tab,main=paste("Anomaly Detection",individs[j]),xlab = "Julian day",ylab = "Number of hits",xlim=c(126,188))
    abline(v=d.vit$juliandropped[j],col=cbPalette[4],lty=2)
}
dev.off()


###
### Evaluating predictions
###

#how many days are anomalies prior to vit drop day
#col1 is total number of days with hits prior to drop+1
#col2 is the total number of days prior to vit drop
#col3 is the ratio of hits to total possible days
#col4 is the logical, whether the the algorithm hit +/- 1 day of vit drop
#col5 is the logical, whether the the algorithm hit +/- 2 days of vit dropth
#col6 is the total number of hits within 1 days
#col7 is the total number of hits within 2 days
#col8 is the total number of hits on vit drop day
#col9 is the total number of hits prior to vit drop day
hit.day.prec=matrix(NA,nr=n.vit,nc = 9)
for(j in 1:n.vit){
    hit.day.prec[j,1]=length(unique(fit.prec$alarm[[j]]$julian[fit.prec$alarm[[j]]$julian <(d.vit$juliandropped[j]-1)]))
    hit.day.prec[j,2]=(d.vit$juliandropped[j]-1)-126+1
    hit.day.prec[j,3]=hit.day.prec[j,1]/hit.day.prec[j,2]
    hit.day.prec[j,4]=sum(fit.prec$alarm[[j]]$julian <=d.vit$juliandropped[j]+1 & fit.prec$alarm[[j]]$julian >=d.vit$juliandropped[j]-1)>0
    hit.day.prec[j,5]=sum(fit.prec$alarm[[j]]$julian <=d.vit$juliandropped[j]+2 & fit.prec$alarm[[j]]$julian >=d.vit$juliandropped[j]-2)>0
    hit.day.prec[j,6]=sum(fit.prec$alarm[[j]]$julian <=d.vit$juliandropped[j]+1 & fit.prec$alarm[[j]]$julian >=d.vit$juliandropped[j]-1)
    hit.day.prec[j,7]=sum(fit.prec$alarm[[j]]$julian <=d.vit$juliandropped[j]+2 & fit.prec$alarm[[j]]$julian >=d.vit$juliandropped[j]-2)
    hit.day.prec[j,8]=sum(fit.prec$alarm[[j]]$julian ==d.vit$juliandropped[j])
    hit.day.prec[j,9]=sum(fit.prec$alarm[[j]]$julian <(d.vit$juliandropped[j]-1))
    }
hit.day.prec

write.csv(hit.day.prec,file="daily_summary_hitday_prec.csv")

hit.day.recall=matrix(NA,nr=n.vit,nc = 9)
#col1 is total number of days with hits prior to day before vit drop
#col2 is the total number of days prior to vit drop
#col3 is the ratio of hits to total possible days
#col4 is the logical, whether the the algorithm hit +/- 1 day of vit drop
#col5 is the logical, whether the the algorithm hit +/- 3 days of vit dropth
#col6 is the total number of hits within 1 days
#col7 is the total number of hits within 3 days
#col8 is the total number of hits on vit drop day
#col9 is the total number of hits prior to vit drop day
for(j in 1:n.vit){
    hit.day.recall[j,1]=length(unique(fit.recall$alarm[[j]]$julian[fit.recall$alarm[[j]]$julian <(d.vit$juliandropped[j]-1)]))
    hit.day.recall[j,2]=(d.vit$juliandropped[j]-1)-126+1
    hit.day.recall[j,3]=hit.day.recall[j,1]/hit.day.recall[j,2]
    hit.day.recall[j,4]=sum(fit.recall$alarm[[j]]$julian <=d.vit$juliandropped[j]+1 & fit.recall$alarm[[j]]$julian >=d.vit$juliandropped[j]-1)>0
    hit.day.recall[j,5]=sum(fit.recall$alarm[[j]]$julian <=d.vit$juliandropped[j]+2 & fit.recall$alarm[[j]]$julian >=d.vit$juliandropped[j]-2)>0
    hit.day.recall[j,6]=sum(fit.recall$alarm[[j]]$julian <=d.vit$juliandropped[j]+1 & fit.recall$alarm[[j]]$julian >=d.vit$juliandropped[j]-1)
    hit.day.recall[j,7]=sum(fit.recall$alarm[[j]]$julian <=d.vit$juliandropped[j]+2 & fit.recall$alarm[[j]]$julian >=d.vit$juliandropped[j]-2)
    hit.day.recall[j,8]=sum(fit.recall$alarm[[j]]$julian ==d.vit$juliandropped[j])
    hit.day.recall[j,9]=sum(fit.recall$alarm[[j]]$julian <(d.vit$juliandropped[j]-1))
}
hit.day.recall

write.csv(hit.day.recall,file="daily_summary_hitday_recall.csv")

fit.recall$alarm[[5]]


Metric=c("total number of days with hits prior to day before vit drop",
         "total number of days prior to vit drop",
         "ratio of hits to total possible days",
         "logical- whether the the algorithm hit +/- 1 day of vit drop",
         "logical- whether the the algorithm hit +/- 3 days of vit dropth",
         "total number of hits within 1 days",
         "total number of hits within 3 days",
         "total number of hits on vit drop day",
         "total number of hits prior to vit drop day")

hit.day.recall[,3]=round(hit.day.recall[,3],2)

out.recall=as.data.frame(cbind(Metric,t(hit.day.recall)),stringsAsFactors=FALSE)re
names(out.recall)=NULL
library(xtable)
xtable(out.recall,include.rownames=FALSE)

hit.day.recall

hit.day.F1=matrix(NA,nr=n.vit,nc = 9)
#col1 is total number of days with hits prior to day before vit drop
#col2 is the total number of days prior to vit drop
#col3 is the ratio of hits to total possible days
#col4 is the logical, whether the the algorithm hit +/- 1 day of vit drop
#col5 is the logical, whether the the algorithm hit +/- 3 days of vit dropth
#col6 is the total number of hits within 1 days
#col7 is the total number of hits within 3 days
#col8 is the total number of hits on vit drop day
#col9 is the total number of hits prior to vit drop day
for(j in 1:n.vit){
    hit.day.F1[j,1]=length(unique(fit.F1$alarm[[j]]$julian[fit.F1$alarm[[j]]$julian <(d.vit$juliandropped[j]-1)]))
    hit.day.F1[j,2]=(d.vit$juliandropped[j]-1)-126+1
    hit.day.F1[j,3]=hit.day.F1[j,1]/hit.day.F1[j,2]
    hit.day.F1[j,4]=sum(fit.F1$alarm[[j]]$julian <=d.vit$juliandropped[j]+1 & fit.F1$alarm[[j]]$julian >=d.vit$juliandropped[j]-1)>0
    hit.day.F1[j,5]=sum(fit.F1$alarm[[j]]$julian <=d.vit$juliandropped[j]+2 & fit.F1$alarm[[j]]$julian >=d.vit$juliandropped[j]-2)>0
    hit.day.F1[j,6]=sum(fit.F1$alarm[[j]]$julian <=d.vit$juliandropped[j]+1 & fit.F1$alarm[[j]]$julian >=d.vit$juliandropped[j]-1)
    hit.day.F1[j,7]=sum(fit.F1$alarm[[j]]$julian <=d.vit$juliandropped[j]+2 & fit.F1$alarm[[j]]$julian >=d.vit$juliandropped[j]-2)
    hit.day.F1[j,8]=sum(fit.F1$alarm[[j]]$julian ==d.vit$juliandropped[j])
    hit.day.F1[j,9]=sum(fit.F1$alarm[[j]]$julian <(d.vit$juliandropped[j]-1))
}
hit.day.F1

write.csv(hit.day.F1,file="daily_summary_hitday_F1.csv")

##############################################################################################

save.image("anomaly_v10_day.Rdata")

