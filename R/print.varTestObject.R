##' @name print.varTestObject
##' @rdname print.varTestObject
##'
##' @title Print basic information about the variance components test
##'
##' @description Displays the likelihood ratio test statistics and the p-value of the test
##'
##' @param x a object of class \code{\link{varTestObject}} obtained from a call to function \code{\link{varTest}}
##'
NULL

print.varTestObject <- function(x){
  if (length(x@chibarsquare) == 0){
    if (length(x@pvalueWeights) == 0){
      uppboundpval <- (1/2)*sum(stats::pchisq(x@lrt,x@cbs@df[(length(x@cbs@df)-1):length(x@cbs@df)],lower.tail = F))
      lowboundpval <- (1/2)*sum(stats::pchisq(x@lrt,x@cbs@df[1:2],lower.tail = F))
      cat("Variance components testing\n LRT = ",x@lrt,
          paste("\n\nlower-bound for p-value:",format(lowboundpval,digits=5)," upper bound for p-value:",format(uppboundpval,digits=5)))
    }else{
      cat("Variance components testing\n LRT = ",x@lrt,
          "\n\np-value (from estimated weights):",format(x@pvalueWeights,digits = 5))
    }
  }else{
    cat("Variance components testing\n LRT = ",x@lrt,
        "\n\np-value (from estimated weights):",format(x@pvalueWeights,digits = 5),
        "\np-value (from simulated chi-bar-square distribution):",format(x@pvalueWeights,digits = 5))
  }
}
