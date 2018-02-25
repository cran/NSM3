print.NSM3Ch9ChickFn <-
  function(x,...){
    if (x$type == "t") {
      print.type <- " two-sided CI for beta:"
      null.dir <- " not equal to "
    } else if (x$type == "l") {
      print.type <- " lower bound for beta:"
      null.dir <- " less than " 
    } else if (x$type == "u"){
      print.type <- " upper bound for beta:"
      null.dir <- " greater than "
    }
      
    cat(paste("Alternative: beta", null.dir, x$beta.0, sep = "")) 
    cat("\n")
    cat(paste("C = ", 
              round(x$C.stat, x$r), 
              ", C.bar = ", 
              round(x$C.bar, x$r), 
              ", P = ", 
              round(x$p.val, x$r), sep = ""))
    
    cat("\n")
    cat(paste("beta.hat = ", round(x$beta.hat, x$r), sep = ""))
    cat("\n")
    cat(paste("alpha.hat = ", round(x$alpha.hat, x$r), sep = ""))
    cat("\n")
    
    if (x$slopes) {
      cat("\n")
      cat("All slopes:")
      cat("\n")
      print(x$slopes.table)
      cat("\n")
    }
    
    cat("\n")
    cat(paste("1 - alpha = ", 1 - x$alpha, print.type, sep = ""))
    cat("\n")
    cat(paste(x$L, ", ", x$U, sep = ""), "\n")
    cat("\n")
  }