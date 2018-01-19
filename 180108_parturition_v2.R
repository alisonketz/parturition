###
### 1/19/2018 
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



setwd("~/Documents/Parturition/180108_parturition")


###
### Load VIT data
###

d.vit=mdb.get('~/Documents/Data/SWDPPdeerDB.MDB',tables= "VIT")
names(d.vit)=tolower(gsub('[[:punct:]]',"",names(d.vit)))
names(d.vit)[10]="lowtag"
d.vit$datedropped=as.character(d.vit$datedropped)
d.vit$juliandropped=yday(mdy_hms(d.vit$datedropped))

n.vit=length(d.vit$lowtag)


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

d$id=as.factor(d$id)

#calculating julian day and omitting outside of parturition window
d$julian=yday(mdy_hms(d$date_time_gmt))

# d=d[d$julian>start & d$julian < end,]

#or keep all data, but exclude first couple hundred observations
#for plotting all data, not just window of parturition

d=data.frame(d %>% group_by(id) %>% slice(300:(n()-300)))
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


#all missing data points are missing for both lats and longs
which(is.na(d$longitude))==which(is.na(d$latitude))
d.spdf=d
d.spdf= SpatialPointsDataFrame(coords=as.data.frame(cbind(d$longitude,d$latitude)),data=d.spdf,proj4string =
                                   CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))

plot(d.spdf)
d.traj <- as.ltraj(coordinates(d.spdf), date=d.spdf$date_time_local,id=d.spdf$id)
head(d.spdf)


###
### Minimum convex polygon
###
d.mcp <- mcp(d.spdf[, "id"], percent = 100,unin = "m", unout = "km2")

# You can plot the home ranges - use the 'col' argument to specify colors
# for different individuals - there are 7 individuals tracked, so we can
# specify we want to use colors 1:7, based on standard R colors. You could
# customize this by specifying indivual colors (one per individual ID) add
# the argument 'axes=TRUE' to add x/y axes based on the data projection
library(RColorBrewer)


plot(d.mcp, col = brewer.pal(12,"Paired",color))
# Custom-specified colors with axes


###
### Kernal Density Estimate
###

d.Khref <- kernelUD(d.spdf[, "id"], h = "href")
# You can plot the KDE for each individual separately:
image(d.Khref)
#can output to raster



###
### Dealing with NA's
###

for(i in 2:records){
    if(is.na(d[i,5]))d[i,5]=d[i-1,5]
    if(is.na(d[i,6]))d[i,6]=d[i-1,6]
    if(is.na(d[i,7]))d[i,7]=d[i-1,7]
    }

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
summary(d.traj$id)
head(d.traj)

hist(d.traj[1], "dist", freq = TRUE)





#converts traj object to data frame
dfdeer <- ld(d.traj)
dfdeer$id=as.character(dfdeer$id)


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


d$rel.angle=dfdeer$rel.angle
d$abs.angle=dfdeer$abs.angle
d$dist=dfdeer$dist
d$R2n=dfdeer$R2n
d$dx=dfdeer$dx
d$dy=dfdeer$dy
d$dt=dfdeer$dt

for(i in 16:22){
    d[is.na(d[,i]),i]=0
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







###
### Looking at moving window correlation of turn angle and distance
###

n.obs=as.numeric(table(d$id))
n.obs


#relative angle      
angle.dist.corr=matrix(NA,nr=n.vit,nc=max(n.obs))
window=18
for(j in 1:n.vit){
    d.tmp=d[d$id==individs[j],]
    for(i in (1+window):dim(d.tmp)[1]){
        angle.dist.corr[j,i]=corr(cbind(d.tmp$rel.angle[(i-window):i],d.tmp$dist[i-window:i]))
    }    
}
angle.dist.corr=angle.dist.corr[,(window+1):dim(angle.dist.corr)[2]]

par(mfrow=c(4,3))
for(j in 1:n.vit){
    plot(angle.dist.corr[j,])
}

angle.dist.corr=angle.dist.corr[,-1]


#absolute angle      
absangle.dist.corr=matrix(NA,nr=n.vit,nc=max(n.obs))
window=18
for(j in 1:n.vit){
    d.tmp=d[d$id==individs[j],]
    for(i in (1+window):dim(d.tmp)[1]){
        absangle.dist.corr[j,i]=corr(cbind(d.tmp$abs.angle[(i-window):i],d.tmp$dist[i-window:i]))
    }    
}
absangle.dist.corr=absangle.dist.corr[,(window+2):dim(absangle.dist.corr)[2]]
par(mfrow=c(4,3))
for(j in 1:n.vit){
    plot(absangle.dist.corr[j,])
}

head(t(angle.dist.corr))

apply(angle.dist.corr,1,min,na.rm=TRUE)
apply(angle.dist.corr,1,min,na.rm=TRUE)





# ### brownian motion movement model
# 
# # The Brownian movement model requires input of some parameters.
# # • sig1: a first smoothing parameter related to the speed of animals;
# # • sig2: a second smoothing parameter related to imprecision of the location data
# 
# likDeer <- liker(d.traj, sig2 = 10000, rangesig1 = c(.001, 10))
# # This works on all individuals at once
# zeb.bb <- kernelbb(zeb.traj, sig1 = c(5.46, 5.16, 5.16, 6.35, 5.76, 6.05, 4.96),
#                    sig2 = 100)
# 
# ## This next line (commented out) works on the first individual (denoted by
# ## the [1]); a list of individuals could be specified using as '[c(1,3,5)]'
# ## zeb.bb.1 <- kernelbb(zeb.traj[1], sig1 = 5.45, sig2 = 100)
# # Again, you can plot the modeled home range as with the Kernel Density
# # Estimates.
# image(zeb.bb)
# 
