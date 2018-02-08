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
### increased fixes parturition window
###

start=yday(mdy("05/05/2017")) # May 1 beginning of parturitionn period
end=yday(mdy("07/07/2017")) #end of parturition period


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
### subset entire dataset to parturition window
###

d=d[d$julian>=start & d$julian <= end,]

#or keep all data, but exclude first couple hundred observations
#for plotting all data, not just window of parturition

# d=data.frame(d %>% group_by(lowtag) %>% slice(200:(n()-200)))
table(d$lowtag)

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



#Create time lag between successive locations to censor data if needed.
timediff <- diff(d$date_time_local)
d=d[-1,]
d$timediff <-round(as.numeric(abs(timediff)))
rm=which(d$timediff>10)
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
plot(d.traj[1])




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
