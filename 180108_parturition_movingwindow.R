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

