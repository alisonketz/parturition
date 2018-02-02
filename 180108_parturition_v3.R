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


d.vit$juliandropped



###
### Load data GPS location data
###


start=yday(mdy("05/05/2017")) # May 1 beginning of parturitionn period
end=yday(mdy("07/07/2017")) #end of parturition period

d = matrix(NA,nr=1,nc=13)
#for loop for reading in data, using vector of id's from the vit dataset
for(i in 1:n.vit){
    d.temp = read.table(paste("/home/aketz/Documents/Data/GPS_locations_ind/",d.vit$lowtag[i],".csv",sep=""),sep=",",header=TRUE,stringsAsFactors = FALSE)
    d.temp$id = d.vit$lowtag[i]
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

table(d$id)
d$id=as.numeric(d$id)
d=d[order(d$id),]

d$id=as.factor(d$id)

###
### Converting lat/long data into UTMs
###

same.onestep=c()
same.twostep=c()
same.threestep=c()
for(j in 4:dim(d)[1]){
    if(!is.na(d$latitude[j]) & !is.na(d$latitude[j-1])){
        if(d$latitude[j]==d$latitude[j-1] & d$longitude[j]==d$longitude[j-1]){
            same.onestep=c(same.onestep,j)
            # if(d$latitude[j]==d$latitude[j-2] & d$longitude[j]==d$longitude[j-2]){
            #     same.twostep=c(same.twostep,j)
            #     if(d$latitude[j]==d$latitude[j-3] & d$longitude[j]==d$longitude[j-3]){
            #         same.threestep=c(same.threestep,j)
            #     }
            # }
        }
    }
}
same.onestep


#Function
LongLatToUTM<-function(x,y,zone){
    xy <- data.frame(ID = 1:length(x), X = x, Y = y)
    browser()
    coordinates(xy) <- c("X", "Y")
    proj4string(xy) <- CRS("+proj=longlat +datum=WGS84")  ## for example
    res <- spTransform(xy, CRS(paste("+proj=utm +zone=",zone," ellps=WGS84",sep='')))
    return(as.data.frame(res))
}

# Example
x<-c( -94.99729,-94.99726,-94.99457,-94.99458,-94.99729)
y<-c( 29.17112, 29.17107, 29.17273, 29.17278, 29.17112)
xy=LongLatToUTM(d$longitude,d$latitude,15)



#calculating julian day and omitting outside of parturition window
d$julian=yday(mdy_hms(d$date_time_gmt))

#subset entire dataset to parturition window
d=d[d$julian>start & d$julian < end,]

#or keep all data, but exclude first couple hundred observations
#for plotting all data, not just window of parturition

# d=data.frame(d %>% group_by(id) %>% slice(200:(n()-200)))
table(d$id)

### Adding Vit data to main GPS dataframe

d$dropped = 0

records=dim(d)[1]
records
for(i in 1:records){
    for(j in 1:n.vit){
        if(d$id[i]==d.vit$lowtag[j]){
            if(d$julian[i]==d.vit$juliandropped[j])d$dropped[i]=1
        }
    }
}
sum(d$dropped)



###
### Formatting for adehabitateLT ltraj object

d$date_time_local=as.POSIXct(strptime(d$date_time_local,format="%m-%d-%Y %H:%M:%S"),tz="CST")

#d=d[!is.na(d$longitude),]


#Using a spatial data frame object

#all missing data points are missing for both lats and longs
# 
sum(is.na(d$longitude))
sum(which(is.na(d$longitude))==which(is.na(d$latitude)))

head(d[is.na(d$longitude),])

# 
# d.spdf=d
# d.spdf= SpatialPointsDataFrame(coords=as.data.frame(cbind(d$longitude,d$latitude)),data=d.spdf,proj4string =
#                                    CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
# 
# plot(d.spdf)
# d.traj <- as.ltraj(coordinates(d.spdf), date=d.spdf$date_time_local,id=d.spdf$id)
# head(d.spdf)


###
### Minimum convex polygon
###
# d.mcp <- mcp(d.spdf[, "id"], percent = 100,unin = "m", unout = "km2")
# 
# # You can plot the home ranges - use the 'col' argument to specify colors
# # for different individuals - there are 7 individuals tracked, so we can
# # specify we want to use colors 1:7, based on standard R colors. You could
# # customize this by specifying indivual colors (one per individual ID) add
# # the argument 'axes=TRUE' to add x/y axes based on the data projection
# library(RColorBrewer)
# 
# 
# plot(d.mcp, col = brewer.pal(12,"Paired",color))
# # Custom-specified colors with axes

# 
# ###
# ### Kernal Density Estimate
# ###
# 
# d.Khref <- kernelUD(d.spdf[, "id"], h = "href")
# # You can plot the KDE for each individual separately:
# image(d.Khref)
# #can output to raster
# 


###
### Dealing with NA's
###
# 
# for(i in 2:records){
#     if(is.na(d[i,5]))d[i,5]=d[i-1,5]
#     if(is.na(d[i,6]))d[i,6]=d[i-1,6]
#     if(is.na(d[i,7]))d[i,7]=d[i-1,7]
#     }

###
### Animal paths
###

d.traj <- as.ltraj(as.data.frame(cbind(d$longitude,d$latitude)), date=d$date_time_local,id=d$id)

# d.traj <- as.ltraj(coordinates(d.spdf), date=d.spdf$date_time_local,id=d.spdf$id)
# View what this looks like
# d.traj
# refda <- strptime("00:00", "%H:%M", tz="US/Central")
# 
# d.traj.set=setNA(d.traj,refda)

#Plot the trajectory for each individual - the numbers correspond to the ID in the d.traj object above
#blue is start
#red is end point

par(mfrow=c(2,2))
for(i in 1:n.vit){
    plot(d.traj[i])
}
d.traj[1]
summary(d$id)
head(d.traj)

hist(d.traj[1], "dist", freq = TRUE)


plot(d.traj[2], "dist", freq = TRUE)

plotltr(d.traj[1], "dt/3600/24")

#converts traj object to data frame
dfdeer <- ld(d.traj)
dfdeer$id=as.character(dfdeer$id)
levels(as.factor(dfdeer$id))

individs=levels(d$id) #individua
class(individs)

hist(dfdeer$rel.angle[dfdeer$id==individs[5]])

hist(dfdeer$rel.angle)

par(mfrow=c(2,2))
for(i in 1:n.vit){
    hist(dfdeer$rel.angle[dfdeer$id==individs[i]])
    plot(dfdeer$rel.angle[dfdeer$id==individs[i]])
    hist(dfdeer$abs.angle[dfdeer$id==individs[i]])
    plot(dfdeer$rel.angle[dfdeer$id==individs[i]])
    }

head(dfdeer)

d$rel.angle=dfdeer$rel.angle
d$dist=dfdeer$dist
d$R2n=dfdeer$R2n
d$dx=dfdeer$dx
d$dy=dfdeer$dy
d$dt=dfdeer$dt

for(i in 16:22){
    d[is.na(d[,i]),i]=2*pi
}
names(d)

#: the squared distance between the first relocation of the trajectory


# and the current relocation is often used to test some movements models
# (e.g. the correlated random walk, see the seminal paper of Kareiva and
# Shigesada, 1983).
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
#     d.tmp=d[d$id==individs[j],]
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
#     d.tmp=d[d$id==individs[j],]
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
#     birth.pred[i]=d$julian[d$id==individs[i]][birth.pred.indx[i]]
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
#     dates_min[i]=as.character(min(d$date_time_local[d$id==individs[i]]))
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

table(d$id)
individs=levels(d$id) #individua
n.obs=as.numeric(table(d$id))

step.dist=matrix(NA,nr=n.vit,nc=max(n.obs)-1)
vit.drop=matrix(NA,nr=n.vit,nc=max(n.obs)-1)
julian.out=matrix(NA,nr=n.vit,nc=max(n.obs)-1)
time.step=matrix(NA,nr=n.vit,nc=max(n.obs)-1)

for(i in 1:n.vit){
  data.tmp=d[d$id==individs[i],]
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
