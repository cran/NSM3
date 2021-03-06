print.NSM3Ch6c <-
function(x,...){
  if(x$method!="Exact"){
    cat(paste0("\n",x$method, " Approximation "))
    if(x$method=="Monte Carlo"){cat(paste0("(with ",x$n.mc, " Iterations) "))}
    cat("used: \n \n")
  }
  
  if(is.null(x$trt)){
    cat("Group sizes:", x$n, "\n")
  }
  
  if(!is.null(x$trt)){
    cat("Control group size:", x$trt, "Treatment group size(s): ", x$n, "\n")
  }
  
  if(x$method=="Asymptotic"){  
    cat(paste0("For the given alpha=", x$alpha, ", the approximate upper cutoff value is ",x$stat.name, "=",x$cutoff.U, ",\n"))
  } else{
    cat(paste0("For the given alpha=", x$alpha, ", the upper cutoff value is ",x$stat.name, "=" ,x$cutoff.U, ",\n", "with true alpha level=",round(x$true.alpha.U,4), "\n"))
  }
  if(!is.null(x$extra)){
    cat(x$extra, "\n")
  }
}
