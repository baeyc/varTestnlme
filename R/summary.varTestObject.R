#' @name summary.varTestObject
#' @rdname summary.varTestObject
#'
#' @title Summary information for the variance components test
#'
#' @description Displays the likelihood ratio test statistics, the limiting distribution and the p-value of the test
#'
#' @param x a object of class \code{\link{varTestObject}} obtained from a call to function \code{\link{varTest}}
#' @return a list containing the following elements:
#' \item{\code{lrt}}{the likelihood ratio test statistics}
#' \item{\code{df}}{the degrees of freedom of the chi-bar distributions involved in the chi-bar-square distribution}
#' \item{\code{weights}}{the weights of the limiting chi-bar-square distribution}
#' \item{\code{pvalWeights}}{the p-value of the test calculated using the cdf of the chi-bar-square based on (approximated) weights}
#' \item{\code{pvalMC}}{the Monte-Carlo estimate of the p-value of the test based on the simulated chi-bar-square distribution}
#'
##'
NULL

summary.varTestObject <- function(x){
  if (length(x@namesTestedParams$fixed)==1 & length(x@namesTestedParams$random)==1) cat(paste("Testing that the mean of",gsub("[()]","",x@namesTestedParams$fixed)," and the variance of",gsub("[()]","",x@namesTestedParams$random),"are null\n"))
  if (length(x@namesTestedParams$fixed)==1 & length(x@namesTestedParams$random)>1) cat(paste("Testing that the mean of",gsub("[()]","",x@namesTestedParams$fixed)," and the variance of",gsub("[()]","",paste(x@namesTestedParams$random,sep="",collapse = " and ")),"are null\n"))
  if (length(x@namesTestedParams$fixed)>1 & length(x@namesTestedParams$random)==1) cat(paste("Testing that the means of",gsub("[()]","",paste(x@namesTestedParams$fixed,sep="",collapse = " and "))," and the variance of",gsub("[()]","",x@namesTestedParams$random),"are null\n"))
  if (length(x@namesTestedParams$fixed)>1 & length(x@namesTestedParams$random)>1) cat(paste("Testing that the means of",gsub("[()]","",paste(x@namesTestedParams$fixed,sep="",collapse = " and "))," and the variance of",gsub("[()]","",paste(x@namesTestedParams$random,sep="",collapse = " and ")),"are null\n"))
  if (length(x@namesTestedParams$fixed)==0 & length(x@namesTestedParams$random)==1) cat(paste("Testing that the variance of",gsub("[()]","",x@namesTestedParams$random),"is null\n"))
  if (length(x@namesTestedParams$fixed)==0 & length(x@namesTestedParams$random)>1) cat(paste("Testing that the variances of",gsub("[()]","",paste(x@namesTestedParams$random,sep="",collapse = " and ")),"are null\n"))

  cat(paste("\nLikelihood ratio test statistics: \n LRT = ",format(x@lrt,digits=5),
            "\n\nLimiting distribution: \n",
            "mixture of",length(x@cbs@df),"chi-bar-square distributions",
            "with degrees of freedom",paste(x@cbs@df,sep="",collapse = ", "),"\n"))

  if (length(x@chibarsquare) == 0){
    if (length(x@pvalueWeights) == 0){
      uppboundpval <- (1/2)*sum(stats::pchisq(x@lrt,x@cbs@df[(length(x@cbs@df)-1):length(x@cbs@df)],lower.tail = F))
      lowboundpval <- (1/2)*sum(stats::pchisq(x@lrt,x@cbs@df[1:2],lower.tail = F))
      cat(paste("\nlower-bound for p-value:",format(lowboundpval,digits=5)," upper bound for p-value:",format(uppboundpval,digits=5)))
    }else{
      cat(" associated weights and sd:",paste(paste(round(x@weights,3)," (",round(x@sdWeights,3),")",sep=""),sep="",collapse = ", "),
          "\n\np-value (from estimated weights):",format(x@pvalueWeights,digits = 5))
    }
  }else{
    cat(" associated weights and sd:",paste(paste(round(x@weights,3)," (",round(x@sdWeights,3),")",sep=""),sep="",collapse = ", "),
        "\n\np-value (from estimated weights):",format(x@pvalueWeights,digits = 5),
        "\np-value (from simulated chi-bar-square distribution):",format(x@pvalueMC,digits = 5))
  }

  result <- list(lrt=x@lrt,df=x@cbs@df,weights=x@weights,pvalWeights=x@pvalueWeights,pvalMC=x@pvalueMC)
  invisible(result)
}

