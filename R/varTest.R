#' Variance component testing
#'
#' Perform a likelihood ratio test to test whether a subset of the variances of the random effects
#' are equal to zero. The test is defined by two hypotheses, H0 and H1, and the model under H0 is
#' assumed to be nested within the model under H1.
#'
#' It is possible to tests if any subset of the variances are equal to zero. However, the function does not
#' currently support nested random effects, and assumes that the random effects are Gaussian.
#' 
#' The asymptotic distribution of the likelihood ratio test is a chi-bar-square, with weights that need to be
#' approximated by Monte Carlo methods, apart from some specific cases where they are available explicitly. 
#' Therefore, the p-value of the test is not exact but approximated. This computation can be time-consuming, so
#' the default behaviour of the function is to provide bounds on the exact p-value, which can be enough in practice
#' to decide whether to reject or not the null hypothesis. This is triggered by the option \code{pval.comp="bounds"}.
#' To compute an approximation of the exact p-value, one should use the option \code{pval.comp="approx"} or \code{pval.comp="both"}.
#' 
#' The \code{control} argument controls the options for chi-bar-square weights computation. It is a list with the
#' following elements: \code{M} the size of the Monte Carlo simulation, \code{parallel} a boolean for parallel computing
#' and \code{nbcores} the number of cores to be used in case of parallel computing. Default is \code{M=5000}, \code{parallel=FALSE}
#' and \code{nbcores=1}.
#'
#' @name varTest
#' @aliases test lrt likelihood
#'
#' @param  m1 a fit of the model under H1, obtained from \code{\link[nlme]{nlme}}, \code{\link[lme4]{lme4}}
#' or \code{\link[saemix]{saemix}}
#' @param m0 a fit of the model under H0, obtained from the same package as \code{m0}
#' @param control (optional) a list of control options for the computation of the chi-bar-weights
#' @param pval.comp (optional) the method to be used to compute the p-value, one of: \code{"bounds"} (the default),
#' \code{"approx"} or \code{"both"} (see Details section)
#' @param fim (optional) the method to compute the Fisher Information Matrix. Currently, only \code{fim="extract"} is supported.
#' @return A list with the following components:
#' \item{\code{lrt}}{the likelihood ratio test statistics}
#' \item{\code{ddl}}{the degrees of freedom of the chi-bar distributions involved in the chi-bar-square distribution}
#' \item{\code{weights}}{the weights of the limiting chi-bar-square distribution}
#' \item{\code{pval}}{the p-value of the test}
#'
#' @author Charlotte Baey <\email{charlotte.baey@univ-lille.fr}>
#'
#' @examples 
#' # load nlme package and example dataset
#' library(nlme)
#' data(Orthodont)
#' 
#' # fit the two models under H1 and H0
#' lm1.h1.nlme <- lme(distance ~ 1 + Sex + age + age*Sex, random = ~ 1 + age | Subject, 
#' data = Orthodont, method = "ML")
#' lm1.h0.nlme <- lme(distance ~ 1 + Sex + age + age*Sex, random = ~ 1 | Subject, 
#' data = Orthodont, method = "ML")
#' 
#' # compare them (order is important: m1 comes first)
#' varTest(lm1.h1.nlme,lm1.h0.nlme)
#'
#' @references Baey C, Cournède P-H, Kuhn E, 2019. Asymptotic distribution of likelihood ratio test
#' statistics for variance components in nonlinear mixed effects models. \emph{Computational
#' Statistics and Data Analysis} 135:107-122.
#'
#' Silvapulle  MJ, Sen PK, 2011. Constrained statistical inference: order, inequality and shape constraints.
#' @export varTest
#' @importFrom stats formula pchisq
varTest <- function(m1,m0,control = list(M=5000),
                    pval.comp = "bounds",
                    fim = "extract"){
  
  # Specify default arguments in control
  if (!is.null(control)) {
    optionNames <- names(control)
    if (!"M" %in% optionNames) control$M=5000
    if (!"parallel" %in% optionNames) control$parallel=F
    if (!"nbcores" %in% optionNames) control$nbcores=1
  }
  
  message("Variance components testing in mixed-effects models\n")
  
  # Identify the packages from which m0 and m1 come from
  cl0 <- class(m0)
  cl1 <- class(m1)
  
  linmodel <- ("lme" %in% cl1 || cl1 %in% "lmerMod")
  if (cl1[1] %in% "nlmerMod") stop("Error: the Fisher information matrix is not available for nonlinear mixed effect models fitted with nlmer() of package lme4. Please use nlme or saemix packages.")
  pkg <- pckName(m1)  
  
  # Identify whether there is any random effects under m0
  randm0 <- TRUE # boolean to indicate whether there are any random effects under H0. This constant is modified if m0 is fitted using lm, glm or nls
  if (cl0[1] %in% c("lm","glm","nls")){ # if there is no random term in m0 it means that we are testing that ALL the variances are null
    
    if (cl0[1] == "lm" & !linmodel) stop("Models should be nested, but m1 is nonlinear while m0 is linear. Please check the formulation.")
    if ((cl0[1] %in% c("glm") || cl0[1] == "nls") & linmodel) stop("Models should be nested, but m1 is linear while m0 is not. Please check the formulation.")
    if (cl1[1] == "nlmerMod" & cl0[1] != "nls") stop("Models should be nested, but m1 is nonlinear while m0 is not. Please check the formulation.")
    if (cl1[1] == "glmerMod" & !(cl0[1] %in% c("glm"))) stop("Models should be nested, but m1 is generalized linear while m0 is not. Please check the formulation.")
    if (pkg == "saemix") warning("We cannot check that the model formulation under H0 and H1 is identical. Please make sure it is the case.")
    
    randm0 <- FALSE
  }else if (cl0[1] != cl1[1]){ # if m0 and m1 both have random effects but are fitted from different package we throw an error
    stop("Error: models m0 and m1 should be fitted with the same function from the same package")
  }
  
  # Check for non nested models in nlme package
  if (cl1[1] == "nlme"){
    if (cl0[1] == "nlme") if((m1$call$model != m0$call$model)) stop("Models should be nested. Check the nonlinear model")
    if (cl0[1] == "nls") if((m1$call$model != m0$call$formula)) stop("Models should be nested. Check the nonlinear model")
  }
  # TODO : the same for saemix and lme4
  
  if (pkg=="nlme") msdata <- modelStructnlme(m1,m0,randm0)
  if (pkg=="lme4") msdata <- modelStructlme4(m1,m0,linmodel,randm0)
  if (pkg=="saemix") msdata <- modelStructsaemix(m1,m0,randm0)
  
  message("(models fitted using the",pkg,"package)\n\n")
  
  # LRT and Fisher Information Matrix
  # LRT
  if (randm0){
    if (pkg=="nlme") lrt <- -2*(m0$logLik - m1$logLik)
    if (pkg=="lme4") lrt <- -2*(stats::logLik(m0) - stats::logLik(m1))
    if (pkg=="saemix") lrt <- -2*(m0@results@ll.is - m1@results@ll.is)
  }else{
    ll0 <- stats::logLik(m0)
    if (pkg=="nlme") lrt <- -2*(ll0 - m1$logLik)
    if (pkg=="lme4") lrt <- -2*(ll0 - stats::logLik(m1))
    if (pkg=="saemix") lrt <- -2*(ll0 - m1@results@ll.is)
  }
  
  # FIM
  if (pkg=="nlme") invfim <- as.matrix(Matrix::bdiag(m1$varFix,m1$apVar)) # error message in case apVar is non positive definite
  if (pkg=="lme4"){
    if (linmodel){
      invfim <- as.matrix(merDeriv::vcov.lmerMod(m1,full=T))
    }else{
      invfim <- as.matrix(merDeriv::vcov.glmerMod(m1,full=T))
    }
  }
  if (pkg=="saemix") invfim <- chol2inv(chol(m1@results@fim))
  
  
  # Constructing chi-bar-square object
  cbs <- methods::new("chiBarSquareObject")
  cbs@orthan <- msdata$structGamma$diag
  
  if (msdata$structGamma$diag || msdata$structGamma$full){
    # reorder FIM so that the tested variances are at the end, before the residual covariance structure
    neworder <- as.numeric(rownames(msdata$detailStruct[order(msdata$detailStruct$tested,msdata$detailStruct$names),]))
    invfim <- invfim[c(neworder,nrow(invfim)-msdata$dims$dimSigma+1),c(neworder,nrow(invfim)-msdata$dims$dimSigma+1)]
    
    if (msdata$structGamma$diag){
      cbs@dims <- list(dimBeta=list(dim0=msdata$dims$nbFE0,
                                    dimR=msdata$dims$nbFE1-msdata$dims$nbFE0),
                       dimGamma=list(dim0=msdata$dims$nbRE0,
                                     dimR=0,
                                     dimSplus=msdata$dims$nbRE1-msdata$dims$nbRE0),
                       dimSigma=msdata$dims$dimSigma)
    }else{
      nbCovTestedAlone <- (msdata$dims$nbRE0==msdata$dims$nbRE1)*(msdata$dims$nbCov1-msdata$dims$nbCov0) # nb of covariances tested without the corresponding variances being tested
      cbs@dims <- list(dimBeta=list(dim0=msdata$dims$nbFE0,
                                    dimR=msdata$dims$nbFE1-msdata$dims$nbFE0),
                       dimGamma=list(dim0=msdata$dims$nbRE0*(msdata$dims$nbRE0+1)/2,
                                     dimR=(msdata$dims$nbRE0)*(msdata$dims$nbRE1-msdata$dims$nbRE0)+nbCovTestedAlone,
                                     dimSplus=msdata$dims$nbRE1-msdata$dims$nbRE0),
                       dimSigma=msdata$dims$dimSigma)
    }
    
  }else if (msdata$structGamma$blockDiag){
    # get indices for the re-ordered parameters : first the fixed effects, then the variances and covariances which are NOT tested,
    # then the covariances tested that are not part of a tested block of sub-matrix, and then the blocks of subsets of variances tested
    dd <- msdata$detailStruct[order(msdata$detailStruct$tested,msdata$detailStruct$block,msdata$detailStruct$covInBlock),]
    neworder <- as.numeric(rownames(dd))
    invfim <- invfim[c(neworder,nrow(invfim)-msdata$dims$dimSigma+1),c(neworder,nrow(invfim)-msdata$dims$dimSigma+1)]
    
    dim0 <- nrow(dd[!dd$tested,])-msdata$dims$nbFE1
    #dimRplus <- length(dd$names[dd$tested & dd$diagBlock])
    dimR <- length(dd$names[dd$covInBlock!=1 & dd$type=="co" & dd$tested])
    ddSp <- dd[dd$tested & dd$covInBlock,]
    dimSplus <- floor(sqrt(2*as.vector(table(ddSp$block))))
    rm(ddSp)
    
    cbs@dims <- list(dimBeta=list(dim0=msdata$dims$nbFE0,
                                  dimR=msdata$dims$nbFE1-msdata$dims$nbFE0),
                     dimGamma=list(dim0=dim0,
                                   dimR=dimR,
                                   dimSplus=dimSplus),
                     dimSigma=msdata$dims$dimSigma)
  }
  cbs@V <- invfim
  cbs@invV <- chol2inv(chol(invfim))
  
  if (length(msdata$nameVarTested)==1){
    message(paste("Testing that the variance of",msdata$nameVarTested,"is null\n"))
  }else{
    message(paste("Testing that the variances of",paste(msdata$nameVarTested,sep="",collapse = " and "),"are null\n"))
  }
  
  
  ## TO COMPLETE
  if (pkg=="lme4"){
    message(paste("model under H1:",deparse(formula(m1),width.cutoff=500),"\n"))
    message(paste("model under H0:",deparse(formula(m0),width.cutoff=500),"\n"))
  }
  if (pkg=="nlme"){
    message(paste("model under H1:",deparse(m1$call$fixed,width.cutoff=500),"(fixed effects)",",",deparse(m1$call$random,width.cutoff=500),"(random effects)\n"))
    if (randm0) message(paste("model under H0:",deparse(m0$call$fixed,width.cutoff=500),"(fixed effects)",",",deparse(m0$call$random,width.cutoff=500),"(random effects)\n"))
    if (!randm0) message(paste("model under H0:",deparse(formula(m0),width.cutoff=500),"(fixed effects)\n"))
  }
  
  # Identify the components of the mixture
  # -> some weights are null, depending on whether the cone includes, or is contained in a linear space
  cbs@df <- dfchisqbar(cbs)
  
  uppboundpval <- (1/2)*sum(stats::pchisq(lrt,cbs@df[(length(cbs@df)-1):length(cbs@df)],lower.tail = F))
  lowboundpval <- (1/2)*sum(stats::pchisq(lrt,cbs@df[1:2],lower.tail = F))
  
  # Compute FIM only if necessary (i.e. when weights of the chi-bar-square need to be computed)
  if (pval.comp != "bounds" & length(cbs@df) > 2){
    if (fim == "extract"){
      message("\nExtracting Fisher Information matrix...")
      if (pkg=="nlme"){
        if (msdata$structGamma$diag) struct <- "diag"
        if (msdata$structGamma$full) struct <- "full"
        if (msdata$structGamma$blockDiag) struct <- "blockDiag"
        apVarTheta <- extractFIMnlme(m1,struct)
        invfim <- as.matrix(Matrix::bdiag(m1$varFix,apVarTheta)) # error message in case apVar is non positive definite
        colnames(invfim) <- rownames(invfim) <- c(names(m1$coefficients$fixed),colnames(apVarTheta))
      }
      if (pkg=="lme4"){
        if (linmodel){
          invfim <- as.matrix(merDeriv::vcov.lmerMod(m1,full=T))
        }else{
          invfim <- as.matrix(merDeriv::vcov.glmerMod(m1,full=T))
        }
      }
      if (pkg=="saemix") invfim <- chol2inv(chol(m1@results@fim))
    }else if (fim == "compute"){
      message("Computation of the FIM is not yet implemented")
    }else{
      stop("Unknown option for fim. Please use fim='extract' or fim='compute'")
    }  
  }
  
  # Compute chi-bar-square weights
  if (length(cbs@df)>1){
    if (pval.comp %in% c("approx","both")){
      wcbs <- weightsChiBarSquare(cbs,control)
      cbs@weights <- wcbs$weights
      # Two ways to compute the pvalue : using the estimated weights (pvalue1) or computing the empirical pvalue (pvalue2)
      pvalue1 <- pchisqbar(lrt[1],cbs,lower.tail = F)
      pvalue2 <- mean(wcbs$randomCBS>=lrt)
      if (min(wcbs$w)<0) warning("\nSome weights were estimated to be negative. Results can be improved by increasing the sampling size M.\n")
    }
  }else{
    pvalue1 <- pchisq(lrt[1],cbs@df[1],lower.tail = F)
  }
  
  # print results
  message(paste("\nLikelihood ratio test statistics: \n LRT = ",format(lrt,digits=5),
            "\n\nLimiting distribution: \n"))
  if (length(cbs@df)>1){
    message(paste("mixture of",length(cbs@df),"chi-bar-square distributions with degrees of freedom",paste(cbs@df,sep="",collapse = ", "),"\n"))
    if (pval.comp %in% c("approx","both")){
      message(paste(" associated weights and sd: ",paste(paste(round(wcbs$w,3)," (",round(wcbs$sdWeights,3),")",sep=""),sep="",collapse = ", "),
                "\n\np-value (from estimated weights):",format(pvalue1,digits = 5)))
      if (length(wcbs$randomCBS)>0) message(paste("\np-value (from simulated chi-bar-square distribution):",format(pvalue2,digits=5),"\n"))
    }
    if (pval.comp %in% c("bounds","both")) message(paste("\nlower-bound for p-value:",format(lowboundpval,digits=5)," upper bound for p-value:",format(uppboundpval,digits=5)))
  }else{
    message(paste0("  chi-bar-square distributions with ",cbs@df," degree of freedom\n"))
    message(paste0("  p-value: ",format(pvalue1,digits=5)))
  }
  
  
  # Creating varTestObject
  vtobj <- varTestObject()
  vtobj@cbs <- cbs
  vtobj@lrt <- lrt[1] # to get rid of the logLik object class
  vtobj@namesTestedParams <- list(fixed=msdata$nameFixedTested,random=msdata$nameVarTested)
  if (pval.comp %in% c("approx","both")){
    if (length(cbs@df)>1) {
      vtobj@chibarsquare <- wcbs$randomCBS
      vtobj@weights <- as.vector(wcbs$weights)
      vtobj@sdWeights <- wcbs$sdWeights
      vtobj@pvalueWeights <- pvalue1
      vtobj@pvalueMC <- pvalue2
    }else{
      vtobj@weights <- cbs@df
      vtobj@sdWeights <- 0
      vtobj@pvalueWeights <- pvalue1
    }
  }
  
  invisible(vtobj)
}
