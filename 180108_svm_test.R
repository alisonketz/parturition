###
### 1/6/2018 
### Alison Ketz
### test for running machine learning of GPS data using svm
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

load("predict.out.Rdata")    

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
### Dealing with NA's
###
# 
# for(i in 2:records){
#     if(is.na(d[i,5]))d[i,5]=d[i-1,5]
#     if(is.na(d[i,6]))d[i,6]=d[i-1,6]
#     if(is.na(d[i,7]))d[i,7]=d[i-1,7]
#     }
# 
# which(is.na(d[,7]))
# d=d[-1,]


head(d[,5:7])
###
### separating individuals with 4 hour step rates
###

individs=levels(d$id) #individuals

individs_step4h = c("5726","7221","6897","5732")
rm.individs=c()
for(i in 1:4){
    rm.individs = c(rm.individs,which(individs==individs_step4h[i]))
}
individs = individs[-rm.individs]

d.sub=d
for(i in 1:4){
    d.sub = d.sub[d.sub$id!=individs_step4h[i],]
}
d.sub$id=as.character(d.sub$id)
table(d.sub$id)
d.sub$id=as.factor(d.sub$id)

n.obs=as.numeric(table(d.sub$id))
n.obs
n.vit=8


###
### Calculating step rates
###
#for all data since beginning of time
n.obs=as.numeric(table(d$id))

n.vit=12
step.rate=matrix(NA,nr=n.vit,nc=max(n.obs)-1)
vit.drop=matrix(NA,nr=n.vit,nc=max(n.obs)-1)
julian.out=matrix(NA,nr=n.vit,nc=max(n.obs)-1)
for(i in 1:n.vit){
    data.tmp=d[d$id==individs[i],]
    for(j in 2:n.obs[i]){
        step.rate[i,j-1]=distHaversine(data.tmp[j-1,6:5],data.tmp[j,6:5]) # Haversine distance in meters
        vit.drop[i,j-1]=data.tmp[j-1,15]
        julian.out[i,j-1]=data.tmp$julian[j-1]        
    }
}


move.rate=step.rate/4 #fixes occured every 2-4 hours
num.obs.day=6 #there are twelve fixes each day


#calculate a moving average for 3 day period. This window could be changed, and we should check for different windows
moverate.mean <- rollapply(move.rate, num.obs.day*3, mean, na.rm=T, by.column=T,partial=T)
dim(moverate.mean)
class(moverate.mean[1,])

vit.records = rep( list(list()), n.vit ) 
for(i in 1:n.vit){
    vit.records[[i]]=which(vit.drop[i,]==1)
}
vit.records





###
### for spring/summer subset
###



step.rate=matrix(NA,nr=n.vit,nc=max(n.obs)-1)
vit.drop=matrix(NA,nr=n.vit,nc=max(n.obs)-1)
julian.out=matrix(NA,nr=n.vit,nc=max(n.obs)-1)

for(i in 1:n.vit){
    data.tmp=d.sub[d.sub$id==individs[i],]
    for(j in 2:n.obs[i]){
        step.rate[i,j-1]=distHaversine(data.tmp[j-1,6:5],data.tmp[j,6:5]) # Haversine distance in meters
        vit.drop[i,j-1]=data.tmp[j-1,15]
        julian.out[i,j-1]=data.tmp$julian[j-1]        
    }
}
move.rate=step.rate/2 #fixes occured every 2 hours
num.obs.day=12 #there are twelve fixes each day


#calculate a moving average for 3 day period. This window could be changed, and we should check for different windows
moverate.mean <- rollapply(move.rate, num.obs.day*3, mean, na.rm=T, by.column=T,partial=T)
dim(moverate.mean)
class(moverate.mean[1,])

vit.records = rep( list(list()), n.vit ) 
for(i in 1:n.vit){
    vit.records[[i]]=which(vit.drop[i,]==1)
}
vit.records

###
### plots of movement rates
###

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#999999", "#E69F00", "#56B4E9", "#009E73")
length(cbPalette)

###
### from changepoint prediction for determining alternative response
###

predict.out$individs=as.character(predict.out$individs)
predict.out$dates_julian_min=as.numeric(predict.out$dates_julian_min)
remove=setdiff(predict.out$individs,individs)
predict.out.sub=predict.out
for(i in 1:4){
    predict.out.sub=predict.out.sub[-which(predict.out.sub$individ==remove[i]),]
}

predict.out.sub
individs_step4h
keep=c(4,12,11,5)

predict.out.star=predict.out
predict.out.star=predict.out.star[keep,]
predict.out.star



pred.out.records = rep(list(list()),8) 
for(i in 1:8){
    pred.out.records[[i]]=which(julian.out[i,]==predict.out.sub$birth.pred[i])
}
pred.out.records


# pdf("moverate1.pdf",width=6,height=8)
par(mfrow=c(4,1))
for(i in 1:4){
    plot(move.rate[i,],col=cbPalette[i],ylab="Movement rate (m/hr)",xlab="Record",main=individs[i])
    lines(moverate.mean[i,])
    abline(v=head(vit.records[[i]],1),col="darkred",lwd=2)
    abline(v=tail(vit.records[[i]],1),col="darkred",lwd=2)
    abline(v=head(pred.out.records[[i]],1),col="darkblue",lwd=2)
    abline(v=tail(pred.out.records[[i]],1),col="darkblue",lwd=2)
}
# dev.off()

# pdf("moverate2.pdf",width=6,height=8)
par(mfrow=c(4,1))
for(i in 5:8){
    plot(move.rate[i,],col=cbPalette[i],ylab="Movement rate (m/hr)",xlab="Record",main=individs[i])
    # lines(moverate.mean[i,])  
    abline(v=head(vit.records[[i]],1),col="darkred",lwd=2)
    abline(v=tail(vit.records[[i]],1),col="darkred",lwd=2)
    abline(v=head(pred.out.records[[i]],1),col="darkblue",lwd=2)
    abline(v=tail(pred.out.records[[i]],1),col="darkblue",lwd=2)
}
# dev.off()



###
### Calculating the step/movement/vit for the 4h vixed collars
###

individs_step4h = c("5726","7221","6897","5732")

d.star=d
dim(d)
for(i in 1:8){
    d.star = d.star[d.star$id!=individs[i],]
}
d.star$id=as.character(d.star$id)
d.star$id=as.factor(d.star$id)

n.obs.star=as.numeric(table(d.star$id))
n.vit.star=4

which(is.na(d.star[,5:7]))


###
### Calculating step rates
###

step.rate.star=matrix(NA,nr=n.vit.star,nc=max(n.obs.star)-1)
vit.drop.star=matrix(NA,nr=n.vit.star,nc=max(n.obs.star)-1)
julian.out.star=matrix(NA,nr=n.vit.star,nc=max(n.obs.star)-1)

for(i in 1:n.vit.star){
    data.tmp=d.star[d.star$id==individs_step4h[i],]
    for(j in 2:n.obs.star[i]){
        step.rate.star[i,j-1]=distHaversine(data.tmp[j-1,6:5],data.tmp[j,6:5]) # Haversine distance in meters
        vit.drop.star[i,j-1]=data.tmp[j-1,15]
        julian.out.star[i,j-1]=data.tmp$julian[j-1]        
    }
}
move.rate.star=step.rate.star/2 #fixes occured every 2 hours
num.obs.day.star=6 #there are twelve fixes each day
head(move.rate.star)

move.rate.star[which(is.na(move.rate.star))]




#calculate a moving average for 3 day period. This window could be changed, and we should check for different windows
moverate.mean.star <- rollapply(move.rate.star, num.obs.day.star*3, mean, na.rm=T, by.column=T,partial=T)

vit.records.star = rep( list(list()), n.vit.star )
pred.out.records.star = rep(list(list()),n.vit.star) 

for(i in 1:n.vit.star){
    vit.records.star[[i]]=which(vit.drop.star[i,]==1)
    pred.out.records.star[[i]]=which(julian.out.star[i,]==predict.out.star$birth.pred[i])
}
vit.records.star

pred.out.records.star




predict.out.star
as.numeric(predict.out.star$birth.pred[i])
min(d.star$julian)
min(predict.out.star$birth.pred)

which(d.star$julian==445,2)

max(d.star$julian)

###
### remove anamalous first observation
###
move.rate.star[4,]=c(move.rate.star[4,-1],NA)
n.obs.star[4]=n.obs.star[4]-1

# pdf("moverate3.pdf",width=6,height=8)
par(mfrow=c(4,1))
for(i in 1:4){
    plot(move.rate.star[i,],col=cbPalette[i+8],ylab="Movement rate (m/hr)",xlab="Record",main=individs_step4h[i])
    lines(moverate.mean.star[i,])
    abline(v=head(vit.records.star[[i]],1),col="darkred")
    abline(v=tail(vit.records.star[[i]],1),col="darkred")
    abline(v=head(pred.out.records.star[[i]],1),col="darkblue",lwd=2)
    abline(v=tail(pred.out.records.star[[i]],1),col="darkblue",lwd=2)
        }
# dev.off()




###
### ACF plots
###
par(mfrow=c(2,2))
for(i in 1:4){
    acf(move.rate[i,1:(n.obs[i]-1)],main=i)
}

par(mfrow=c(2,2))
for(i in 5:8){
    acf(move.rate[i,1:(n.obs[i]-1)],main=i)
}


par(mfrow=c(2,2))
for(i in 1:4){
    acf(move.rate.star[i,1:(n.obs.star[i]-1)],main=paste("Star ",i))
}

par(mfrow=c(2,2))
for(i in 1:4){
    acf(moverate.mean[i,],main=i)
}

par(mfrow=c(2,2))
for(i in 5:8){
    acf(moverate.mean[i,],main=i)
}


par(mfrow=c(2,2))
for(i in 1:4){
    acf(move.rate.star[i,],main=paste("star",i),na.rm=TRUE)
}
which(is.na(move.rate.star))

###
### mapping points using ggmap
###

###
### Overall map of all individuals

# Make a bounding box for data
box <- make_bbox(longitude,latitude,data = d)
calc_zoom(box)
map=get_map(location=c(mean(d$longitude),mean(d$latitude)),source="google",zoom="auto",maptype="satellite",crop=box)
#Map from URL : http://maps.googleapis.com/maps/api/staticmap?center=43.069755,-90.280971&zoom=14&size=640x640&scale=2&maptype=satellite&language=en-EN&sensor=false

pdf("Map_overall_Parturition")
ggmap(map)+geom_point(aes(x = longitude, y = latitude,colour=id),data = d, size = 3)
dev.off()
###
### mapping points for each individual
###

n.map=length(unique(d$devid))
k = unique(d$devid)
d$dropped=as.factor(d$dropped)

for(i in 1:n.map){
    d.temp=d[d$devid==k[i],]
    reflong=mean(d.temp$longitude)
    reflat=mean(d.temp$latitude)
    n.leg=dim(d.temp)[1]
    x=d.temp$longitude
    xend=c(d.temp$longitude[2:n.leg],d.temp$longitude[n.leg])
    y=d.temp$latitude
    yend=c(d.temp$latitude[2:n.leg],d.temp$latitude[n.leg])
    df.leg=data.frame(x,xend,y,yend)
    box <- make_bbox(longitude,latitude,data = d.temp)
    map=get_map(location=c(reflong,reflat),source="google",zoom=15,maptype="satellite",crop=box)
    pdf(file=paste("Map_",i,".pdf",sep=""))
    print(ggmap(map)+geom_point(aes(x = longitude, y = latitude,colour=dropped),data = d.temp, size = 2,alpha=.25)+
              scale_colour_manual(name = "", values = c(cbPalette[1],cbPalette[2]))+
              geom_point(aes(x = longitude, y = latitude),data = d.temp[d.temp$dropped=="1",], size = 2,alpha=1,color=cbPalette[2])+
              ggtitle(d.temp[1,1]))
    dev.off()
}

###
### Histogram of step rates 
###


step.dist=step.rate
dim(step.rate)
dim(step.dist)
step.event=matrix(NA,nr=n.vit,nc=50)
indx1.vec=c()
indx2.vec=c()
for(i in 1:n.vit){
    indx1=head(vit.records[[i]],1)
    indx2=tail(vit.records[[i]],1)
    step.event[i,]=step.dist[i,indx1:(indx1+49)]
    indx1.vec=c(indx1.vec,indx1)
    indx2.vec=c(indx2.vec,indx2)
}

step.event
apply(step.event,1,mean)
hist(step.event)

step.dist=step.rate

par(mfrow=c(1,1))
quantile(step.rate,c(.05,.1,.15,.2,.5,.99),na.rm = TRUE)

par(mfrow=c(1,2))
hist(step.rate[step.rate<851],breaks=30)
hist(step.event,breaks=30)


vit.drop

which(is.na(step.dist[3,]))
step.rate
chpt.out=list(rep(list(),n.vit))
par(mfrow=c(2,1))
for(i in 1:n.vit){
    chpt.out[[i]]=cpt.meanvar(na.trim(step.dist[i,]))
    plot(chpt.out[[i]],main=individs[i])
}
chpt.out=cpt.mean(step.dist)

plot(chpt.out[[1]])

###
### Preliminary exploratory analysis
###


out=list()
out
for(i in 1:8){
    out[[i]]=naiveBayes(y=vit.drop[i,],x=moverate.mean[i,])
}


print(out)


length(t(vit.drop))
length(t(moverate.mean))

svm(c(vit.drop)~c(step.rate))
vit.s=c(vit.drop)
move.s=c(moverate.mean)

svm.out=svm(c(vit.drop)~c(step.rate))


plot(predict(svm.out))

###
### Why are some of these plots of points ending so much less than the max
###

#answer: some of the interevals are 4 hour intervals, some are 2 hour intervals


d.mort=mdb.get('~/Documents/Data/SWDPPdeerDB.MDB',tables= "Mortalities")
names(d.mort)=tolower(gsub('[[:punct:]]',"",names(d.mort)))
for(i in 1:n.vit){
    if(d.mort$lowtag==individs[i])cat(individs[i],'\t')
}


class(individs)
class(d.mort$lowtag)="character"
class(d.mort$lowtag)

cat(which(d.mort$lowtag==individs[i]),'\t')


short1=d[d$id==5732,]#only monitoring every 4 hours
short1


