print.NSM3Ch7MCp <-function(x,...){
    if(!is.null(x$ties)){cat("Ties are present, so p-values are based on conditional null distribution. \n")}
    
    if(is.null(x$trt)){
      cat(paste0("Number of blocks: n=", x$n,"\n"))
      cat(paste0("Number of treatments: k=", x$k, "\n"))
    } else{
      cat("Control group size:", x$trt, "Treatment group size(s):", x$n, "\n")
    }
    
    cat("Using the", x$method)
    if(x$method=="Monte Carlo"){
      cat(" (with", x$n.mc, "Iterations) ")
    }    
    cat("method: \n \n")
    for(i in 1:x$num.comp){
      round.p.val<-ifelse(round(x$p.val[i],4)==0,x$p.val[i],round(x$p.val[i],4))
      cat(paste0("For treatments ", x$labels[i], ", the ",x$stat.name," Statistic is ", round(x$obs.stat[i],4),". \n"))
      cat(paste0("The smallest experimentwise error rate leading to rejection is ",round.p.val ,".\n  \n"))
    }
      
    if(!is.null(x$extra)){
      cat(x$extra, "\n")
    }
}