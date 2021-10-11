#' Useful intern functions
#' 
#' @title null.desc
#' @description create null.value description
#' @param msdata a list containing the structure of the model and data, as an output from 
#' \code{extractStruct.<package_name>} functions
null.desc <- function(msdata){
  null.value <- rep(0,1+length(msdata$nameFixedTested))
  if (length(msdata$nameFixedTested)==0){
    if (length(msdata$nameVarTested)==1){
      names(null.value) <- paste("variance of the random effect associated to",msdata$nameVarTested)
    }else{
      names(null.value) <- paste("covariance matrix of",paste(msdata$nameVarTested,collapse = " and ")) 
    }
  }else{
    names(null.value) <- c(paste("mean of the random effect associated to",paste(msdata$nameFixedTested,collapse = " and ")),
                           paste(ifelse(length(msdata$nameVarTested)==1,"variance of ","covariance matrix of "),paste(msdata$nameVarTested,collapse = " and ")))
  }
  return(null.value)
}

#' @title alt.desc
#' @description create alternative description
#' @param msdata a list containing the structure of the model and data, as an output from 
#' \code{extractStruct.<package_name>} functions
alt.desc <- function(msdata){
  if (length(msdata$nameFixedTested)==0){
    if (length(msdata$nameVarTested)==1){
      alternative=c(paste("variance of the random effect associated to",msdata$nameVarTested,"> 0 "))
    }else{
      alternative <- paste("covariance matrix of",paste(msdata$nameVarTested,collapse = " and "),"> 0 ")
    }
  }else{
    alternative <- c(paste("mean of the random effect associated to",paste(msdata$nameFixedTested,collapse = " is different from 0 and ")),
                     paste(" and ",ifelse(length(msdata$nameVarTested)==1,"variance of ","covariance matrix of "),paste(msdata$nameVarTested,collapse = " and "),"> 0 "))
  }
  return(alternative)
}

#' @title print.desc.message
#' @description print a message to indicate the null and alternative hypotheses
#' @param msdata a list containing the structure of the model and data, as an output from 
#' \code{extractStruct.<package_name>} functions
print.desc.message <- function(msdata){
  if (length(msdata$nameFixedTested)==0){
    if (length(msdata$nameVarTested)==1){
      message(paste("Testing that the variance of the random effect associated to",msdata$nameVarTested,"is equal to 0"))
    }else if (length(msdata$nameVarTested) > 1){
      message(paste("Testing that the covariance matrix of",paste(msdata$nameVarTested,sep="",collapse = " and "),"is equal to 0\n"))
    }else{
      covTested <- msdata$detailStruct[msdata$detailStruct$tested,]
      namesToPrint <- paste(covTested$var1," and ",covTested$var2)
      message("Testing that covariances between the random effects ",paste(namesToPrint,collapse=", "),
              ifelse(length(namesToPrint)==1," is"," are")," equal to 0")
    }
  }else if (length(msdata$nameFixedTested)==1){
    message("Testing that the mean of the random effect associated to ",msdata$nameFixedTested," is equal to 0 and that")
    if (length(msdata$nameVarTested)==1){
      message(paste(" the variance of the random effect associated to",msdata$nameVarTested,"is equal to 0"))
    }else if (length(msdata$nameVarTested) > 1){
      message(paste("the covariance matrix of",paste(msdata$nameVarTested,sep="",collapse = " and "),"is equal to 0\n"))
    }else{
      covTested <- msdata$detailStruct[msdata$detailStruct$tested,]
      namesToPrint <- paste(covTested$var1," and ",covTested$var2)
      message("Testing that covariances between the random effects ",paste(namesToPrint,collapse=", "),
              ifelse(length(namesToPrint)==1," is"," are")," equal to 0")
    }
  }else{
    message("Testing that the means of the random effects associated to ",paste(msdata$nameFixedTested,sep="",collapse = " and ")," are equal to 0 and that")
    if (length(msdata$nameVarTested)==1){
      message(paste(" the variance of the random effect associated to",msdata$nameVarTested,"is equal to 0"))
    }else if (length(msdata$nameVarTested) > 1){
      message(paste("the covariance matrix of",paste(msdata$nameVarTested,sep="",collapse = " and "),"is equal to 0\n"))
    }else{
      covTested <- msdata$detailStruct[msdata$detailStruct$tested,]
      namesToPrint <- paste(covTested$var1," and ",covTested$var2)
      message("Testing that covariances between the random effects ",paste(namesToPrint,collapse=", "),
              ifelse(length(namesToPrint)==1," is"," are")," equal to 0")
    }
  }
}

#' @title print.res.message
#' @description print a message with the results
#' @param results an object of class vctest
print.res.message <- function(results){
  pval_sample <- !is.na(results$p.value["pvalue.sample"])
  pval_weights <- !is.na(results$p.value["pvalue.weights"])
  
  message("Likelihood ratio test statistic:\n\tLRT = ",format(results$statistic,5))
  
  if (pval_weights | pval_sample){
    if (results$p.value["pvalue.lowerbound"] == results$p.value["pvalue.upperbound"]){
      message("\np-value from exact weights: ",format(results$p.value["pvalue.weights"],5))
    }else{
      if (pval_weights) message("\np-value from estimated weights: ",format(results$p.value["pvalue.weights"],5))
      if (pval_sample) message("\np-value from random sample: ",format(results$p.value["pvalue.sample"],5))
      
      message("bounds on p-value: lower ",format(results$p.value["pvalue.lowerbound"],5),
              "\tupper ",format(results$p.value["pvalue.upperbound"],5))  
    }
  }else{
    if (results$p.value["pvalue.lowerbound"] != results$p.value["pvalue.upperbound"]){
      message("bounds on p-value: lower ",format(results$p.value["pvalue.lowerbound"],5),
          "\tupper ",format(results$p.value["pvalue.upperbound"],5))  
    }else{
      message("exact p-value: ",format(results$p.value["pvalue.lowerbound"],5)) 
    }
  }
  
  message("\n")
}