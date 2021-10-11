# methods for 'vctest' objects 

#' @title Extract the Fisher Information Matrix
#' @param object an object of class vctest
#' @export fim.vctest
fim.vctest <- function(object)  {
  object$parameters$FIM
}

#' @title Print
#' @param x an object of class vctest
#' @param ... additional arguments
#' @export print.vctest
#' @export
print.vctest <- function(x, ...) 
{
  cat("Variance components testing in mixed effects models\n")
  cat("Testing that:\n",paste(paste(names(x$null.value),"is equal to 0"),collapse="\n "))
  cat("\nagainst the alternative that:\n",paste(x$alternative,collapse="\n "))
  
  lrt <- x$statistic
  pval_sample <- !is.na(x$p.value["pvalue.sample"])
  pval_weights <- !is.na(x$p.value["pvalue.weights"])
  
  cat("\n\n","Likelihood ratio test statistic:\n\tLRT = ",format(lrt,5))
  
  if (pval_weights | pval_sample){
    if (max(x$parameters$sdweights) == 0){
      cat("\n\n","exact p-value:",x$p.value["pvalue.weights"])
    }else{
      if (pval_weights) cat("\n\n","p-value from estimated weights:",x$p.value["pvalue.weights"])
      if (pval_sample) cat("\n","p-value from random sample:",x$p.value["pvalue.sample"])
    }
  }else{
    if (x$p.value["pvalue.lowerbound"] != x$p.value["pvalue.upperbound"]){
      cat("\n\n","bounds on p-value: lower ",format(x$p.value["pvalue.lowerbound"],5),
          "upper ",format(x$p.value["pvalue.upperbound"],5))  
      
      message("\n\nBounds based on the smallest and biggest degrees of freedom of the chi-bar-square distribution components. Re-run with option 'pval.comp=\"both\" or pval.comp=\"comp\" to approximate the weights of each chi-bar-square component and the p-value.")
    }else{
      cat("\n\n","exact p-value: ",format(x$p.value["pvalue.lowerbound"],5)) 
    }
  }
  cat("\n")
  invisible(x)
}

#' @title Summary
#' @param object an object of class vctest
#' @param ... additional arguments
#' @export summary.vctest
#' @export
summary.vctest <- function(object, ...)
{
  cat("Variance components testing in mixed effects models\n")
  cat("Testing that:\n",paste(paste(names(object$null.value),"is equal to 0"),collapse="\n "))
  cat("\nagainst the alternative that:\n",paste(object$alternative,collapse="\n "))
  
  lrt <- object$statistic
  w <- object$parameters$weights
  df <- object$parameters$df
  sdw <- object$parameters$sdweights
  pval_sample <- !is.na(object$p.value["pvalue.sample"])
  pval_weights <- !is.na(object$p.value["pvalue.weights"])
  
  cat("\n\n","Likelihood ratio test statistic:\n\tLRT = ",format(lrt,5))
  cat("\n\n","Limiting distribution:")
  if (length(df) > 1){
    cat("\n\tmixture of",length(df),"chi-bar-square distributions with degrees of freedom",df)
    if (!prod(is.na(w))){
      if (max(sdw) > 0) cat("\n\tassociated weights (and sd):",paste0(format(w,5)," (",sdw,")"))
      if (max(sdw) == 0) cat("\n\tassociated (exact) weights:",paste0(format(w,5)))
    }
    
    cat("\n\n","p-value of the test:")
    if (pval_weights){
      if (max(sdw) == 0){
        cat("\n\tfrom exact weights:",object$p.value["pvalue.weights"])
      }else{
        cat("\n\tfrom estimated weights:",object$p.value["pvalue.weights"])
      }
    }
    if (pval_sample) cat("\n\tfrom random sample:",object$p.value["pvalue.sample"])
    if (object$p.value["pvalue.lowerbound"] != object$p.value["pvalue.upperbound"]){
      cat("\n\tbounds on p-value: lower ",format(object$p.value["pvalue.lowerbound"],5),
          "upper ",format(object$p.value["pvalue.upperbound"],5)) 
      
      message("\n\nBounds based on the smallest and biggest degrees of freedom of the chi-bar-square distribution components. Re-run with option 'pval.comp=\"both\" or pval.comp=\"comp\" to approximate the weights of each chi-bar-square component and the p-value.")
    }
  }else{
    cat("chi-bar-square distributions with ",df," degree of freedom\n")
    cat("p-value: ",format(object$p.value["pvalue.weights"],5))
  }
  
  cat("\n\n")
  invisible(object)
}

