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

keep.window=NULL
d=d[d$julian>start & d$julian < end,]

### Adding Vit data to main GPS dataframe

d$dropped = 0

records=dim(d)[1]

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

for(i in 2:records){
    if(is.na(d[i,6]))d[i,6]=d[i-1,6]
    if(is.na(d[i,7]))d[i,7]=d[i-1,7]
    }

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


n.obs.star



cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#999999", "#E69F00", "#56B4E9", "#009E73")
length(cbPalette)

pdf("moverate1.pdf",width=6,height=8)
par(mfrow=c(4,1))
for(i in 1:4){
    plot(move.rate[i,],col=cbPalette[i],ylab="Movement rate (m/hr)",xlab="Record",main=individs[i])
    lines(moverate.mean[i,])
    abline(v=head(vit.records[[i]],1),col="darkred",lwd=2)
    abline(v=tail(vit.records[[i]],1),col="darkred",lwd=2)
}
dev.off()

pdf("moverate2.pdf",width=6,height=8)
par(mfrow=c(4,1))
for(i in 5:8){
    plot(move.rate[i,],col=cbPalette[i],ylab="Movement rate (m/hr)",xlab="Record",main=individs[i])
    lines(moverate.mean[i,])  
    abline(v=head(vit.records[[i]],1),col="darkred",lwd=2)
    abline(v=tail(vit.records[[i]],1),col="darkred",lwd=2)
}
dev.off()




###
### Calculating the step/movement/vit for the 4h vixed collars
###

individs_step4h = c("5726","7221","6897","5732")

d.star=d
for(i in 1:8){
    d.star = d.star[d.star$id!=individs[i],]
}
d.star$id=as.character(d.star$id)
d.star$id=as.factor(d.star$id)

n.obs.star=as.numeric(table(d.star$id))
n.vit.star=4


###
### Filling in NA's for missing lat/longs
###

for(i in 2:records){
    if(is.na(d.star[i,6]))d.star[i,6]=d.star[i-1,6]
    if(is.na(d.star[i,7]))d.star[i,7]=d.star[i-1,7]
}

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

#calculate a moving average for 3 day period. This window could be changed, and we should check for different windows
moverate.mean.star <- rollapply(move.rate.star, num.obs.day.star*3, mean, na.rm=T, by.column=T,partial=T)

vit.records.star = rep( list(list()), n.vit.star ) 
for(i in 1:n.vit.star){
    vit.records.star[[i]]=which(vit.drop.star[i,]==1)
}
vit.records.star





pdf("moverate3.pdf",width=6,height=8)
par(mfrow=c(4,1))
for(i in 1:4){
    plot(move.rate.star[i,],col=cbPalette[i+8],ylab="Movement rate (m/hr)",xlab="Record",main=individs_step4h[i])
    lines(moverate.mean.star[i,])
    abline(v=head(vit.records.star[[i]],1),col="darkred")
    abline(v=tail(vit.records.star[[i]],1),col="darkred")
}
dev.off()

###
### Preliminary analysis
###

svm(t(vit.drop)~t(moverate.mean))
vit.s=c(vit.drop)
move.s=c(moverate.mean)

svm.out=svm(vit.s~move.s)

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
