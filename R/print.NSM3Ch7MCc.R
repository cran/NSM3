print.NSM3Ch7MCc <-
  function(x,...){
    if(x$method!="Exact"){
      cat(paste0("\n",x$method, " Approximation"))
      if(x$method=="Monte Carlo"){cat("(with ",x$n.mc, " Iterations)")}
      cat("used: \n \n")
    }
    
    if(is.null(x$trt)){
      cat(paste0("Number of blocks: n=", x$n,"\n"))
      cat(paste0("Number of treatments: k=", x$k, "\n"))
    } else{
      cat("Control group size:", x$trt, "Treatment group size(s):", x$n, "\n")
    }
    
    if(x$method!="Asymptotic"){  
      cat(paste0("For the given experimentwise alpha=", x$alpha, ", the upper cutoff value is ",x$stat.name, "=" ,x$cutoff.U, ", with true alpha level=",round(x$true.alpha.U,4), "\n") )
    }
    if(x$method=="Asymptotic"){
      cat(paste0("For the given experimentwise alpha=", x$alpha, ", the approximate upper cutoff value is ",x$stat.name, "=",x$cutoff.U, ",\n"))
    }
    if(!is.null(x$extra)){
      cat(x$extra, "\n")
    }
  }