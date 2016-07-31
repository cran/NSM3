print.NSM3Ch7c <-
  function(x,...){
    if(x$method!="Exact"){
      cat(paste0("\n",x$method, " Approximation "))
      if(x$method=="Monte Carlo"){cat("(with",x$n.mc, "Iterations) ")}
      cat("used: \n \n")
    }
    
    if(is.null(x$trt)){
      cat(paste0("Number of blocks: n=", x$n,"\n"))
      cat(paste0("Number of treatments: k=", x$k, "\n"))
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
    
    if(!is.null(x$trt)){
      cat("Control group size:", x$trt, "Treatment group size(s):", x$n, "\n")
    }
    
    if(x$method=="Asymptotic"){  
      cat(paste0("For the given alpha=", x$alpha, ", the approximate upper cutoff value is ",x$stat.name, "=",x$cutoff.U, ",\n"))
    } else{
      cat(paste0("For the given alpha=", x$alpha, ", the upper cutoff value is ",x$stat.name, "=" ,x$cutoff.U, ", with true alpha level=",round(x$true.alpha.U,4), "\n"))
    }
    if(!is.null(x$extra)){
      cat(x$extra, "\n")
    }
  }