#' @rdname varCompTest
#' @importFrom stats formula pchisq
#' @rawNamespace export(varCompTest.lme)
#' @export
varCompTest.lme <- function(m1,m0,control = list(M=5000,parallel=T,nb_cores=1,B=1000),pval.comp = "bounds",fim = "extract"){
  
  # Specify default arguments in control
  if (!is.null(control)) {
    optionNames <- names(control)
    if (!"M" %in% optionNames) control$M=5000
    if (!"parallel" %in% optionNames) control$parallel=T
    if (!"nbcores" %in% optionNames) control$nbcores=1
    if (!"B" %in% optionNames) control$B = 1000
  }
  
  message("Variance components testing in mixed effects models")
  
  pkg <- pckName(m1)
  randm0 <- !(max(class(m0) %in% c("lm","glm","nls"))) # are there any random effect under H0?
  
  # Extract data structure
  msdata <- extractStruct(m1,m0,randm0)
  
  if (length(msdata$nameVarTested)==1){
    message(paste("Testing that the variance of",msdata$nameVarTested,"is null\n"))
  }else if (length(msdata$nameVarTested) > 1){
    message(paste("Testing that the variances of",paste(msdata$nameVarTested,sep="",collapse = " and "),"are null\n"))
  }else{
    message(paste("Testing that covariances ",paste0(msdata$detailStruct$names[msdata$detailStruct$tested],collapse=", "),"are null"))
  }
  
  # Compute LRT
  ll1 <- m1$logLik
  ll0 <- ifelse(randm0,m0$logLik,stats::logLik(m0))
  lrt <- -2*(ll0 - ll1)

  
  # Degrees of freedom of the chi-square components
  cbs.df.dims <- dfChiBarSquare(msdata)
  
  # FIM to compute the weights
  if (pval.comp != "bounds" & (length(cbs.df.dims$df) > 2)){
    if (fim == "extract"){
      invfim <- extractFIM.lme(m1,msdata$structGamma) # error message in case apVar is non positive definite
    }else if (fim == "compute"){
      message("Computing Fisher Information Matrix by bootstrap...\n")
      invfim <- bootinvFIM(m1, control$B)
    }else if (is.matrix(fim)){
      invfim <- chol2inv(fim)
    }else{
      stop("Unknown option for fim. Please use fim='extract' or fim='compute'")
    }
    
    # re-order FIM to match the definition of the convex cone in the chi-bar-square computation
    if (msdata$struct %in% c("diag","full")){
      neworder <- as.numeric(rownames(msdata$detailStruct[order(msdata$detailStruct$tested,msdata$detailStruct$names),]))
      invfim <- invfim[c(neworder,nrow(invfim)-msdata$dims$dimSigma+1),c(neworder,nrow(invfim)-msdata$dims$dimSigma+1)]
    }else{
      dd <- msdata$detailStruct[order(msdata$detailStruct$tested,msdata$detailStruct$block,msdata$detailStruct$covInBlock),]
      neworder <- as.numeric(rownames(dd))
      invfim <- invfim[c(neworder,nrow(invfim)-msdata$dims$dimSigma+1),c(neworder,nrow(invfim)-msdata$dims$dimSigma+1)]
    }
    
    fim <- chol2inv(chol(invfim))
  }else{
    invfim <- fim <- NA
  }

  # Compute chi-bar-square weights and p-value
  if (length(cbs.df.dims$df)>1){
    if (pval.comp %in% c("approx","both")){
      cbs.weights.sample <- weightsChiBarSquare(df=cbs.df.dims$df,
                                         V=invfim,
                                         dimsCone=cbs.df.dims$dimsCone,
                                         orthan=(msdata$structGamma == "diag"),
                                         control=control)
      
      pvalue1 <- sum(cbs.weights.sample$weights * stats::pchisq(lrt,df=cbs.df.dims$df,lower.tail = F))   # p-value from weights
      pvalue2 <- mean(cbs.weights.sample$randomCBS >= lrt)                                      # p-value from random sample
      
      if (min(cbs.weights.sample$weights)<0) warning("\nSome weights were estimated to be negative. Results can be improved by increasing the sampling size M.\n")
    }else{
      pvalue1 <- NA
      pvalue2 <- NA
      cbs.weights.sample <- list(weights=NA,sdWeights=NA,randomCBS=NA)
    }
    if (length(cbs.df.dims$df)==2){
      cbs.weights.sample <- list(weights=c(0.5,0.5),sdWeights=c(0,0),randomCBS=NA)
      pvalue1 <- sum(cbs.weights.sample$weights * stats::pchisq(lrt,df=cbs.df.dims$df,lower.tail = F))   # p-value from weights
      pvalue2 <- NA
    }
  }else{
    pvalue1 <- stats::pchisq(lrt,cbs.df.dims$df[1],lower.tail = F)
    pvalue2 <- NA
    cbs.weights.sample <- list(weights=NA,sdWeights=NA,randomCBS=NA)
  }
  
  # Bounds on p-value
  uppboundpval <- (1/2)*sum(stats::pchisq(lrt,cbs.df.dims$df[(length(cbs.df.dims$df)-1):length(cbs.df.dims$df)],lower.tail = F))
  lowboundpval <- (1/2)*sum(stats::pchisq(lrt,cbs.df.dims$df[1:2],lower.tail = F))
  
  # Print results
  message(paste("Likelihood ratio test statistic: \n LRT = ",format(lrt,digits=5),
                "\n\nLimiting distribution:"))
  if (length(cbs.df.dims$df) > 1){
    message(paste("mixture of",length(cbs.df.dims$df),"chi-bar-square distributions with degrees of freedom",paste(cbs.df.dims$df,sep="",collapse = ", "),"\n"))
    
    if (pval.comp %in% c("approx","both")){
      message(paste(" associated weights and sd: ",paste(paste(round(cbs.weights.sample$weights,3)," (",round(cbs.weights.sample$sdWeights,3),")",sep=""),sep="",collapse = ", "),
                    "\n\np-value (from estimated weights):",format(pvalue1,digits = 5)))
      
      if (length(cbs.weights.sample$randomCBS) > 0) message(paste("\np-value (from simulated chi-bar-square distribution):",format(pvalue2,digits=5),"\n"))
    }
    
    if (pval.comp %in% c("bounds","both")) message(paste("lower-bound for p-value:",format(lowboundpval,digits=5)," upper bound for p-value:",format(uppboundpval,digits=5)))
  }else{
    message(paste0("  chi-bar-square distributions with ",cbs.df.dims$df," degree of freedom\n"))
    message(paste0("  p-value: ",format(pvalue1,digits=5)))
  }
  
  # create results, object of class htest
  null.value <- rep(0,length(msdata$nameVarTested)+length(msdata$nameFixedTested))
  if (length(msdata$nameFixedTested)==0){
    names(null.value) <- paste("variance of",msdata$nameVarTested)
    alternative=c(paste("variance of",msdata$nameVarTested,"> 0 "))
  }else{
    names(null.value) <- c(paste("variance of",msdata$nameVarTested),paste(" mean of",msdata$nameFixedTested)) 
    alternative=c(paste("variance of",msdata$nameVarTested,"> 0 "),paste(" mean of",msdata$nameFixedTested," different from 0 "))
  }

  results <- list(statistic=c(LRT=lrt),
                  null.value=null.value,
                  alternative=alternative,
                  parameters=list(df=cbs.df.dims$df,weights=cbs.weights.sample$weights,FIM=fim),
                  method="Likelihood ratio test for variance components in mixed effects models",
                  p.value=c(pvalue.weights=pvalue1,pvalue.sample=pvalue2,pvalue.lowerbound=lowboundpval,pvalue.upperbound=uppboundpval))
  class(results) <- "htest"
  
  invisible(results)
}
