print.NSM3Ch7p <-function(x,...){
  if(!is.null(x$ties)&&x$ties){cat("Ties are present, so p-values are based on conditional null distribution. \n")}
    
  if(is.null(x$trt)){
    cat(paste0("Number of blocks: n=", x$n,"\n"))
    cat(paste0("Number of treatments: k=", x$k, "\n"))
  } else{
    cat("Control group size:", x$trt, "Treatment group size(s):", x$n, "\n")
  }
  if(!is.null(x$ss)){
    cat(paste0("Number of treatments per block: s=", x$ss,"\n"))
  }
  if(!is.null(x$pp)){
    cat(paste0("Number of observations per treatment: p=", x$pp,"\n"))
  }
  if(!is.null(x$lambda)){
    cat(paste0("Number of times each pair of treatments occurs together within a block: lambda=", x$lambda,"\n"))
  }
  cat(x$stat.name, "Statistic:", round(x$obs.stat,4) , "\n")
  cat(x$method, "")
    
  if(x$method=="Monte Carlo"){
    cat("(Using ", x$n.mc, "Iterations) ")
  }
    
  cat("upper-tail probability:", ifelse(round(x$p.val,4)==0,x$p.val,round(x$p.val,4)), "\n")
    
  if(!is.null(x$extra)){
    cat(x$extra, "\n")
  }
    
}