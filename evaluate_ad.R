
evaluate.time.pop = function(alarm,possible.hits,n.vit=12,vitdropday){
    
    if(n.vit==1){
        alarm=alarm[[1]][alarm[[1]]$julian<=vitdropday,]
        
        # TRUE POSITIVE:  groundtruth data  says it's an anomaly and so algorithm does.
        tp = sum(alarm$dropped,na.rm=TRUE)
        
        # FALSE POSITIVE:  groundtruth data says it's not an anomaly, but algorithm says anomaly.
        fp = sum(alarm$dropped==0,na.rm=TRUE)
        
        # FALSE NEGATIVE: groundtruth data says it's an anomaly, but algorithm says not anomaly.
        fn = possible.hits - tp
        
        tp.pop = tp
        fp.pop = fp
        fn.pop = fn
        
        # precision and recall
        out.prec = tp.pop/(tp.pop+fp.pop)
        out.recall = tp.pop/(tp.pop+fn.pop)
        
        # F1 value
        out.F1 = (2*out.prec*out.recall)/(out.prec+out.recall)
    }
    else{
        for(j in 1:n.vit){
            alarm[[j]]=alarm[[j]][alarm[[j]]$julian<=vitdropday[j],]
        }
        
        tp = rep(NA,n.vit)
        fp = rep(NA,n.vit)
        fn = rep(NA,n.vit)
        
        for(j in 1:n.vit){
            
            # TRUE POSITIVE:  groundtruth data  says it's an anomaly and so algorithm does.
            tp[j] = sum(alarm[[j]]$dropped,na.rm=TRUE)
            
            # FALSE POSITIVE:  groundtruth data says it's not an anomaly, but algorithm says anomaly.
            fp[j] = sum(alarm[[j]]$dropped==0,na.rm=TRUE)
            
            # FALSE NEGATIVE: groundtruth data says it's an anomaly, but algorithm says not anomaly.
            fn[j] = possible.hits[j] - tp[j]
            
        }
        
        tp.pop = sum(tp)
        fp.pop = sum(fp)
        fn.pop = sum(fn)
        
        # precision and recall
        out.prec = tp.pop/(tp.pop+fp.pop)
        out.recall = tp.pop/(tp.pop+fn.pop)
        
        # F1 value
        out.F1 = (2*out.prec*out.recall)/(out.prec+out.recall)
    }#endelse
    
    return(list(tp=tp,fp=fp,fn=fn,out.prec=out.prec,out.recall=out.recall,out.F1=out.F1,tp.pop=tp.pop,fp.pop=fp.pop,fn.pop=fn.pop))
    
}
