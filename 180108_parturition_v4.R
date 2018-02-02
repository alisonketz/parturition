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
###
### weird unix windows conversion
###


d=data.frame(d,stringsAsFactors = FALSE)

for(j in 1:dim(d)[2]){
    d[,j]=str_trim(d[,j],side="both")
}

class.indx=c(5:7,9:12)
for(j in class.indx){
    d[,j]=as.numeric(d[,j])
}

d$lowtag=as.numeric(d$lowtag)
d=d[order(d$lowtag),]
d$lowtag=as.factor(d$lowtag)

###
### double checking for duplicate locations
###

summary(duplicated(d))

#calculating julian day and omitting outside of parturition window
d$julian=yday(mdy_hms(d$date_time_local))

#subset entire dataset to parturition window
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
head(d$latitude)
head(d$longitude)
head(d$missing)
sum(d$missing)

dim(d)

miss.per.ind=c()
for(j in 1:n.vit){
    miss.per.ind=c(miss.per.ind,sum(d$missing[d$lowtag==individs[j]]))
}


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
d=d[-1,]
dist.out = distHaversine(d[,6:5])
d$distance = dist.out

#remove the initial index of each individual 
#because it makes no sense to calculate distances 
#between the locations of the individuals

remove.indx = rep(NA,n.vit)
for(j in 1:n.vit){
    remove.indx[j] = min(which(d$lowtag==individs[j]))
}
remove.indx = remove.indx[-1]

d=d[-remove.indx]

d=d[-c(dim(d)[1]-1,dim(d)[1]),]

tail(d)


###
### Projections! 
###

# points from scratch
coords = cbind(x, y)
sp = SpatialPoints(coords)
# make spatial data frame
spdf = SpatialPointsDataFrame(coords, data)
spdf = SpatialPointsDataFrame(sp, data)
# promote data frame to spatial
coordinates(data) = cbind(x, y)
coordinates(data) = ~lon + lat


d.sp = 

newData <- spTransform(d[6,5], CRS("+proj=tmerc +lat_0=0 +lon_0=-90 +k=0.9996 +x_0=520000
                                   +y_0=-4480000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))





# 
# 
# 
# ###
# ### Looking at moving window correlation of turn angle and distance
# ###
# 
# n.obs=as.numeric(table(d$id))
# n.obs
# 
# 
# #relative angle      
# angle.dist.corr=matrix(NA,nr=n.vit,nc=max(n.obs))
# vit.drop=matrix(NA,nr=n.vit,nc=max(n.obs))
# julian.out=matrix(NA,nr=n.vit,nc=max(n.obs))
# window=18
# for(j in 1:n.vit){
#     d.tmp=d[d$lowtag==individs[j],]
#     for(i in (1+window):dim(d.tmp)[1]){
#         angle.dist.corr[j,i]=corr(cbind(d.tmp$rel.angle[(i-window):i],d.tmp$dist[i-window:i]))
#         vit.drop[j,i]=d.tmp[i,15]
#         julian.out[j,i]=d.tmp$julian[i] 
#         }    
# }
# angle.dist.corr=angle.dist.corr[,(window+2):dim(angle.dist.corr)[2]]
# vit.drop=vit.drop[,(window+2):dim(vit.drop)[2]]
# julian.out=julian.out[,(window+2):dim(julian.out)[2]]
# 
# 
# 
# 
# vit.records = rep( list(list()), n.vit ) 
# for(i in 1:n.vit){
#     vit.records[[i]]=which(vit.drop[i,]==1)
# }
# vit.records
# 
# 
# par(mfrow=c(4,3))
# for(j in 1:n.vit){
#     plot(angle.dist.corr[j,])
#     abline(v=head(vit.records[[j]],1),col="darkred",lwd=2)
#     abline(v=tail(vit.records[[j]],1),col="darkred",lwd=2)
#     }
# 
# # angle.dist.corr=angle.dist.corr[,-1]
# 
# 
# #absolute angle      
# absangle.dist.corr=matrix(NA,nr=n.vit,nc=max(n.obs))
# 
# window=18
# for(j in 1:n.vit){
#     d.tmp=d[d$lowtag==individs[j],]
#     for(i in (1+window):dim(d.tmp)[1]){
#         absangle.dist.corr[j,i]=corr(cbind(d.tmp$abs.angle[(i-window):i],d.tmp$dist[i-window:i]))
#     }    
# }
# absangle.dist.corr=absangle.dist.corr[,(window+2):dim(absangle.dist.corr)[2]]
# 
# 
# par(mfrow=c(4,3))
# for(i in 1:n.vit){
#     plot(absangle.dist.corr[i,])
#     abline(v=head(vit.records[[i]],1),col="darkred",lwd=2)
#     abline(v=tail(vit.records[[i]],1),col="darkred",lwd=2)
# }
# 
# head(t(angle.dist.corr))
# 
# which(angle.dist.corr==apply(angle.dist.corr,1,min,na.rm=TRUE),arr.ind = TRUE)
# which(absangle.dist.corr==apply(absangle.dist.corr,1,min,na.rm=TRUE),arr.ind = TRUE)
# 
# 
# 
# 
# ###
# ###
# ###
# 
# plot(loess(c(angle.dist.corr)~c(julian.out)))
# 
# 
# apply(angle.dist.corr,2,sd,na.rm=TRUE)
# 
# sd(angle.dist.corr,na.rm=TRUE)
# par(mfrow=c(4,3))
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
# 
# 
# 
# 
# 
# birth.pred.indx
# birth.pred
# par(mfrow=c(1,1))
# plot(d$julian)
# max(d$julian)
# 
# min(d$julian)
# min(d$date_time_gmt)
# dates_min=list()
# dates_julian_min=c()
# for(i in 1:n.vit){
#     dates_min[i]=as.character(min(d$date_time_local[d$lowtag==individs[i]]))
#     dates_julian_min[i]=min(d$julian[d$id==individs[i]])
#     }
# 
# dates_min
# dates_julian_min
# 
# predict.out=as.data.frame(cbind(individs,birth.pred,dates_julian_min),stringsAsFactors = FALSE)
# 
# save(predict.out,file="predict.out.Rdata")
# 

      
###
### Calculating haversine distance, when missing values, use distance of points on either end
###
      
      
  
###
### Calculating step rates
###
#for all data since beginning of time

table(d$lowtag)
individs=levels(d$lowtag) #individua
n.obs=as.numeric(table(d$lowtag))

step.dist=matrix(NA,nr=n.vit,nc=max(n.obs)-1)
vit.drop=matrix(NA,nr=n.vit,nc=max(n.obs)-1)
julian.out=matrix(NA,nr=n.vit,nc=max(n.obs)-1)
time.step=matrix(NA,nr=n.vit,nc=max(n.obs)-1)

for(i in 1:n.vit){
  data.tmp=d[d$lowtag==individs[i],]
  for(j in 2:n.obs[i]){
      if(!is.na(data.tmp[j,5])){
          step.dist[i,j-1]=distVincentyEllipsoid(data.tmp[j-1,6:5],data.tmp[j,6:5]) # Haversine distance in meters
          vit.drop[i,j-1]=data.tmp[j-1,15]
          julian.out[i,j-1]=data.tmp$julian[j-1]
          time.step[i,j-1]=round(data.tmp[j,4]-data.tmp[j-1,4])
      }
      else{
          step.dist[i,j-1]=distVincentyEllipsoid(data.tmp[j-1,6:5],data.tmp[j+1,6:5]) # Haversine distance in meters
          vit.drop[i,j-1]=data.tmp[j-1,15]
          julian.out[i,j-1]=data.tmp$julian[j-1]
          time.step[i,j-1]=round(data.tmp[j+1,4]-data.tmp[j-1,4])
      }
   
  }
}

vit.julian=d.vit$juliandropped[order(d.vit$lowtag)]
fawn.julian=d.vit$julianfawnfound[order(d.vit$lowtag)]

 
par(mfrow=c(4,3))
for(i in 1:12){
    hist(time.step[i,])
}

par(mfrow=c(4,3))
for(i in 1:12){
    hist(step.dist[i,])
}

      
par(mfrow=c(4,3))
for(i in 1:n.vit){
    plot(julian.out[i,],step.dist[i,],xlab="Julian Day")
    abline(v=head(vit.julian[i],1),col="darkred",lwd=2)
    abline(v=head(fawn.julian[i],1),col="darkblue",lwd=2)
    }
hist(step.dist[i,],breaks = 50)

#changing vit.julian because we can tell that the VIT fell out early. Consequently, guessing for the purpose of this training set

vit.julian[9]=146
par(mfrow=c(4,3))
for(i in 1:n.vit){
    plot(julian.out[i,],step.dist[i,],xlab="Julian Day")
    abline(v=head(vit.julian[i],1),col="darkred",lwd=2)
    abline(v=head(fawn.julian[i],1),col="darkblue",lwd=2)
}


min(step.dist,na.rm=T)

event.window=rep(list(),n.vit)
event.julian=rep(list(),n.vit)
for(i in 1:n.vit){
    temp.indx=which(julian.out[i,]==vit.julian[i])
    event.window[[i]]=step.dist[i,(min(temp.indx)-30):(max(temp.indx)+30)]
    event.julian[[i]]=julian.out[i,(min(temp.indx)-30):(max(temp.indx)+30)]    
    
}


event.window
par(mfrow=c(4,3))
for(i in 1:n.vit){
    plot(event.julian[[i]],event.window[[i]],xlab="Julian Day",main=individs[i])
    abline(v=head(vit.julian[i],1),col="darkred",lwd=2)
    abline(v=head(fawn.julian[i],1),col="darkblue",lwd=2)
}

step.dist[1,]
