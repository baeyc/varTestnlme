#' @rdname varCompTest
#' @importFrom stats formula pchisq
#' @rawNamespace export(varCompTest.SaemixObject)
#' @export
varCompTest.SaemixObject <- function(m1,m0,control = list(M=5000,parallel=FALSE,nb_cores=1,B=1000),pval.comp = "bounds",fim = "extract"){
  
  # Specify default arguments in control
  if (!is.null(control)) {
    optionNames <- names(control)
    if (!"M" %in% optionNames) control$M=5000
    if (!"parallel" %in% optionNames) control$parallel=FALSE
    if (!"nbcores" %in% optionNames) control$nbcores=1
    if (!"B" %in% optionNames) control$B = 1000
  }
  
  message("Variance components testing in mixed effects models")
  
  # Identify the packages from which m0 and m1 come from
  randm0 <- !(class(m0) %in% c("lm","glm","nls")) # are there any random effect under H0?
  
  # Extract data structure
  msdata <- extractStruct(m1,m0,randm0)

  # Print message
  print.desc.message(msdata)
  
  # Compute LRT
  lrt <- -2*(saemix::logLik.SaemixObject(m0) - saemix::logLik.SaemixObject(m1))
  
  
  # Degrees of freedom of the chi-square components
  cbs.df.dims <- dfChiBarSquare(msdata)
  
  # FIM to compute the weights
  if (pval.comp != "bounds" & (length(cbs.df.dims$df) > 2)){
    if (fim == "extract"){
      invfim <- chol2inv(chol(m1@results@fim))
    }else if (fim == "compute"){
      invfim <- bootinvFIM(m1,control$B)
    }else if (is.matrix(fim)){
      invfim <- chol2inv(fim)
    }else{
      stop("Unknown option for fim. Please use fim='extract' or fim='compute'")
    }
    
    # re-order FIM
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
    cbs.weights.sample <- list(weights=1,sdWeights=0,randomCBS=NA)
  }
  
  # Bounds on p-value
  uppboundpval <- (1/2)*sum(stats::pchisq(lrt,cbs.df.dims$df[(length(cbs.df.dims$df)-1):length(cbs.df.dims$df)],lower.tail = F))
  lowboundpval <- (1/2)*sum(stats::pchisq(lrt,cbs.df.dims$df[1:2],lower.tail = F))
  
  # create results, object of class htest
  null.value <- null.desc(msdata)
  alternative <- alt.desc(msdata)
  
  results <- list(statistic=c(LRT=lrt),
                  null.value=null.value,
                  alternative=alternative,
                  parameters=list(df=cbs.df.dims$df,weights=cbs.weights.sample$weights,sdweights=cbs.weights.sample$sdWeights,FIM=fim),
                  method="Likelihood ratio test for variance components in mixed effects models",
                  p.value=c(pvalue.weights=pvalue1,pvalue.sample=pvalue2,pvalue.lowerbound=lowboundpval,pvalue.upperbound=uppboundpval))
  class(results) <- c("vctest","htest")
  
  print.res.message(results)
  
  invisible(results)
}
