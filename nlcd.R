###
### National landcover joinin
###


###
### Preliminaries
###

rm(list=ls())


library(raster)
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
library(rgdal)


setwd("~/Documents/Parturition/180108_parturition")

#######################################################################################################
###
### Load VIT data
###
#######################################################################################################
d.vit=mdb.get('~/Documents/Data/SWDPPdeerDB.MDB',tables= "VIT")
names(d.vit)=tolower(gsub('[[:punct:]]',"",names(d.vit)))
names(d.vit)[10]="lowtag"
d.vit$datedropped=as.character(d.vit$datedropped)
d.vit$juliandropped=yday(mdy_hms(d.vit$datedropped))
d.vit$datefawnfound=as.character(d.vit$datefawnfound)
d.vit$julianfawnfound=yday(mdy_hms(d.vit$datefawnfound))
n.vit=length(d.vit$lowtag)

#reorder vit data by lowtag/lowtag
d.vit=d.vit[order(d.vit$lowtag),]

#extract the individual ids
individs=d.vit$lowtag

#manipulate the julian day for the doe 6865 with best guess based on movement rates
d.vit$juliandropped[9] = 147


#######################################################################################################
###
### Load data GPS location data
###
#######################################################################################################

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

#calculating julian day and omitting outside of parturition window
d$julian=yday(mdy_hms(d$date_time_local))

# increased fixes parturition window
start=yday(mdy("03/20/2017")) # May 6 beginning of parturition period
end=yday(mdy("07/07/2017")) #end of parturition period

# subset entire dataset to parturition window
d=d[d$julian>start & d$julian <= end,]

#removing last day of parturition of 5004, outlier movements
rm5004 = (dim(d[d$lowtag==5004,])[1]-15):dim(d[d$lowtag==5004,])[1]
d=d[-rm5004,]

d$dropped = 0
records=dim(d)[1]
for(i in 1:records){
    for(j in 1:n.vit){
        if(d$lowtag[i]==d.vit$lowtag[j]){
            if(d$julian[i]==d.vit$juliandropped[j])d$dropped[i]=1
        }
    }
}

# Converting date to POSIXct format
d$date_time_local=as.POSIXct(strptime(d$date_time_local,format="%m-%d-%Y %H:%M:%S"),tz="CST6CDT")
d$date_time_gmt=as.POSIXct(strptime(d$date_time_gmt,format="%m-%d-%Y %H:%M:%S"),tz="GMT")

#Create time lag between successive locations to censor data if needed.
time.diff <- diff(d$date_time_local)
time.diff
d=d[-1,]
d$timediff <-round(as.numeric(abs(time.diff)))
rm=which(d$timediff>10)
d=d[-rm,]
names(d)[1]="devname"

###
### Dealing with missing data
###

#impute with midpoint
d$missing=0
for( i in 1:dim(d)[1]){
    if(is.na(d$longitude[i]))d$missing[i]=1
}
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
d=data.frame(d.sp.proj)#20 and 21 are the coordinates in UTMs x y

###
### Animal paths
### Using the adelehabitatLT

d.traj <- as.ltraj(d[,20:21], date=d$date_time_local,id=d$lowtag)
dfdeer <- ld(d.traj) #converts traj object to data frame
dfdeer$id=as.character(dfdeer$id)
dfdeer=rbind(rep(NA,dim(dfdeer)[2]),dfdeer)
dfdeer=dfdeer[-dim(dfdeer)[1],]
d$rel.angle=dfdeer$rel.angle
d$dist.traj=dfdeer$dist
d$R2n=dfdeer$R2n
d$dx=dfdeer$dx
d$dy=dfdeer$dy
d$dt=dfdeer$dt

remove.indx=which(is.na(d$dist.traj))
d[remove.indx,]
d=d[-remove.indx,]

#######################################################################################################3
###
### Adding Landcover data
###
#######################################################################################################3

### Convert Data back to spatial points data frame

# load points data
class(d)
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

# plot(nlcd.img)
# crs(nlcd.img)

# load cover data
NLCD <- raster("~/Documents/Data/nlcd_2011_landcover_2011_edition_2014_10_10/nlcd_2011_landcover_2011_edition_2014_10_10.img")

# convert lat/lon to appropriate projection
crs_args <- NLCD@crs@projargs

d.sp.proj.nlcd=spTransform(d.sp.proj,crs_args)
    
#clip the NLCD to the d.sp.proj size +buffer
b <- bbox(d.sp.proj.nlcd)
b[,1]=b[,1]*.99
b[,2]=b[,2]*1.01

NLCD.crop=crop(NLCD,b)
plot(NLCD.crop)
plot(d.sp.proj.nlcd,add=T)    

#extract land cover data for each point, given buffer size
Landcover <- extract(NLCD.crop, d.sp.proj.nlcd)
num.codes <- unique(unlist(Landcover))
cover.names <- NLCD@data@attributes[[1]]$NLCD.2011.Land.Cover.Class[num.codes + 1]
levels(cover.names)[1] <- NA # first level is ""
conversions <- data.frame(num.codes, cover.names)
conversions <- na.omit(conversions)
conversions <- conversions[order(conversions$num.codes),]
conversions

# # summarize each site's data by proportion of each cover type
summ <- lapply(Landcover, function(x){prop.table(table(x))})

# generate land cover number to name conversions
num.codes <- unique(unlist(Landcover))
cover.names <- NLCD@data@attributes[[1]]$NLCD.2011.Land.Cover.Class[num.codes + 1]
levels(cover.names)[1] <- NA # first level is ""
conversions <- data.frame(num.codes, cover.names)
conversions <- na.omit(conversions)
conversions <- conversions[order(conversions$num.codes),]
# convert to data frame
nlcd.points = data.frame(d.sp.proj.nlcd,cover2 = names(unlist(summ)))

# create cover name column
nlcd.points$cover <- nlcd.points$cover2
levels(nlcd.points$cover) <- conversions$cover.names
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
head(d)
names(d)[rm.columns]
rm.columns=c(29:32,34:36)
d=d[,-rm.columns]


# ###
# ### Import Wiscland 2.0
# ###
# par(mfrow=c(1,1))
# wisc.level1 <- raster("~/Documents/Data/wiscland2.0/level1/wiscland2_level1.tif")
# plot(wisc.level1)
# summary(wisc.level1)
# wisc.level2 <- raster("~/Documents/Data/wiscland2.0/level2/wiscland2_level2.tif")
# plot(wisc.level2)
# summary(wisc.level2)
# 
# wisc.level3 <- raster("~/Documents/Data/wiscland2.0/level3/wiscland2_level3.tif")
# plot(wisc.level3)
# summary(wisc.level3)
# 
# wisc.level4 <- raster("~/Documents/Data/wiscland2.0/level4/wiscland2_level4.tif")
# plot(wisc.level4)
# summary(wisc.level1)