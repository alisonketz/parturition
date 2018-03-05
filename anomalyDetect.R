
anomalyDetect = function(n.vit,part.window=126,id,d,eps,covs.indx){
    individual=id
    #For individual case, run anomaly dection on single individual
    if(n.vit==1){
        cv.indx = covs.indx
        numcovs = length(cv.indx)
        epsilon=eps
        ind.mean=rep(NA,numcovs)
        ind.sd=rep(NA,numcovs)
        d.temp=d
        for(k in 1:numcovs){
            ind.mean[k]=mean(d.temp[d.temp$julian<part.window,cv.indx[k]],na.rm=TRUE)
            ind.sd[k]=sd(d.temp[d.temp$julian<part.window,cv.indx[k]],na.rm=TRUE)
        }
        
        d.temp=d[d$julian>part.window,]
        n.out=dim(d.temp)[1]
        n.out.max=n.out
        
        if(numcovs==1){
            ind.mean=c(ind.mean)
            ind.sd=c(ind.sd)
            source("anomalyDetectUniv.R")
            individual=id
            single=anomalyDetectUniv(n.out,n.vit=1,part.window=part.window,id=individual,d,eps=epsilon,im=ind.mean,is=ind.sd,cov.indx = cv.indx)
            alarm=single$alarm
            detect.density=single$detect.density
            threshold.density=single$threshold.density
            detect.quant=single$detect.quant
            }
        else{
            detect.quant=rep(NA,numcovs)
            detect.density=rep(NA,n.out)
            alarm=rep(list(),n.vit)
            d.temp=d[d$julian>part.window,]
            n.temp=dim(d.temp)[1]
            for(k in 1:numcovs){
                detect.quant[k]=quantile(d.temp[,cv.indx[k]],epsilon[k],na.rm=TRUE)
            }
            threshold.density = dmvnorm(detect.quant,ind.mean,diag(ind.sd))
            for(i in 1:n.temp){
                detect.density[i] = dmvnorm(d.temp[i,covs.indx],ind.mean,diag(ind.sd))
            }
            alarm= list(d.temp[detect.density<threshold.density,])
        }#end else numcovs
    }#endif n.vit=1
    else {# if n.vit>1
        cv.indx = covs.indx
        numcovs = length(cv.indx)
        epsilon=eps
        ind.mean=matrix(NA,nr=numcovs,nc=n.vit)
        ind.sd=matrix(NA,nr=numcovs,nc=n.vit)
        for(j in 1:n.vit){
            d.temp=d[d$lowtag==id[j],]
            for(k in 1:numcovs){
                ind.mean[k,j]=mean(d.temp[d.temp$julian<part.window,cv.indx[k]],na.rm=TRUE)
                ind.sd[k,j]=sd(d.temp[d.temp$julian<part.window,cv.indx[k]],na.rm=TRUE)
            }
        }
        
        n.out=c()
        for(j in 1:n.vit){
            d.temp=d[d$lowtag==id[j] & d$julian>part.window,]
            n.temp=dim(d.temp)[1]
            n.out=c(n.out,n.temp)
        }
        n.out.max=max(n.out)

        if(numcovs==1){
            ind.mean=c(ind.mean)
            ind.sd=c(ind.sd)
            epsilon=eps
            indiv=id
            nv=n.vit
            source("anomalyDetectUniv.R")
            single=anomalyDetectUniv(n.out.max,n.vit=nv,window=part.window,id=indiv,d,eps=epsilon,im=ind.mean,is=ind.sd,cov.indx=covs.indx)
            alarm=single$alarm
            detect.density=single$detect.density
            threshold.density=single$threshold.density
            detect.quant=single$detect.quant
            }
        
        else{
            detect.quant=matrix(NA,nr=numcovs,nc=n.vit)
            epsilon=eps
            threshold.density=rep(NA,n.vit)
            detect.density=matrix(NA,nr=n.out.max,nc=n.vit)
            alarm=rep(list(),n.vit)
            for(j in 1:n.vit){
                d.temp1=d[d$lowtag==id[j],]
                d.temp=d.temp1[d.temp1$julian>part.window,]
                n.temp=dim(d.temp)[1]
                for(k in 1:numcovs){
                    detect.quant[k,j]=quantile(d.temp[,cv.indx[k]],epsilon[k],na.rm=TRUE)
                }
                threshold.density[j] = dmvnorm(detect.quant[,j],ind.mean[,j],diag(ind.sd[,j]))
                for(i in 1:n.temp){
                    detect.density[i,j] = dmvnorm(d.temp[i,cv.indx],ind.mean[,j],diag(ind.sd[,j]))
                }
                alarm[[j]]= d.temp[is.finite(detect.density[,j])&detect.density[,j]<threshold.density[j],]
            }#endfor        
        }#end else
        
    }#end group n.vit>1 else

    return(list(alarm=alarm,detect.density=detect.density,threshold.density=threshold.density,detect.quant=detect.quant))
}
