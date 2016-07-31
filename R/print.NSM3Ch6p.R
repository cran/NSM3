print.NSM3Ch6p <-
function(x,...){
	if(!is.null(x$ties)&&x$ties){cat("Ties are present, so p-values are based on conditional null distribution. \n")}
  
  if(is.null(x$trt)){
    cat("Group sizes:", x$n, "\n")
  }
	
  if(!is.null(x$trt)){
	  cat("Control group size: ", x$trt, "Treatment group size(s): ", x$n, "\n")
  }
  
  
  cat(paste0(x$stat.name, " Statistic: ", round(x$obs.stat,4) , "\n"))
	cat(paste0(x$method, " "))
	
  if(x$method=="Monte Carlo"){
		cat(paste0("(Using ", x$n.mc, " Iterations) "))
	}
  
	cat(paste0("upper-tail probability: ", ifelse(round(x$p.val,4)==0,x$p.val,round(x$p.val,4)), "\n"))
	
  if(!is.null(x$extra)){
	  cat(x$extra, "\n")
	}
  
}
