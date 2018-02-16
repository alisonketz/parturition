anomalyDetectUniv = function(n.out.max,n.vit,window,id,d,eps,im,is,cov.indx){
    #can only be run on a single covariate
    if(n.vit==1){
        covariateindx=cov.indx
        epsilon=eps
        d.temp=d[d$julian>window,]
        n.temp=dim(d.temp)[1]
        detect.quant=quantile(d.temp[,covariateindx],epsilon,na.rm=TRUE)
        threshold.density = dnorm(detect.quant,im,is)
        for(i in 1:n.temp){
            detect.density[i] = dnorm(d.temp[i,covariateindx],im,is)
        }
        alarm= list(d.temp[detect.density<threshold.density,])
    }else{
        covariateindx=cov.indx
        detect.quant=rep(NA,n.vit)
        threshold.density=rep(NA,n.vit)
        detect.density=matrix(NA,nr=n.out.max,nc=n.vit)
        alarm=rep(list(),n.vit)
        epsilon=eps
        for(j in 1:n.vit){
            d.temp1=d[d$lowtag==id[j],]
            d.temp=d.temp1[d.temp1$julian>window,]
            n.temp=dim(d.temp)[1]
            detect.quant[j]=quantile(d.temp[,covariateindx],epsilon,na.rm=TRUE)
            threshold.density[j] = dnorm(detect.quant[j],im[j],is[j])
            for(i in 1:n.temp){
                detect.density[i,j] = dnorm(d.temp[i,covariateindx],im[j],is[j])
            }
            alarm[[j]]= d.temp[detect.density[,j]<threshold.density[j],]
        }
    }
    return(list(alarm=alarm,detect.density=detect.density,threshold.density=threshold.density,detect.quant=detect.quant))
}