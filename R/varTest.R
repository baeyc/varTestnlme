#' Variance component testing
#'
#' Perform a likelihood ratio test to test whether a subset of the variances of the random effects
#' are equal to zero. The test is defined by two hypotheses, H0 and H1, and the model under H0 is
#' assumed to be nested within the model under H1.
#'
#' It is possible to tests if any subset of the variances are equal to zero. However, the function does not
#' currently support nested random effects, and assumes that the random effects are Gaussian.
#'
#' @name varTest
#' @aliases test lrt likelihood
#'
#' @param  m1 a fit of the model under H1, obtained from \code{\link[nlme]{nlme}}, \code{\link[lme4]{lme4}}
#' or \code{\link[saemix]{saemix}}
#' @param m0 a fit of the model under H0, obtained from the same package as \code{m0}
#' @param control (optional) a list of control options for the computation of the chi-bar-weights
#' @return A list with the following components:
#' \item{\code{lrt}}{the likelihood ratio test statistics}
#' \item{\code{ddl}}{the degrees of freedom of the chi-bar distributions involved in the chi-bar-square distribution}
#' \item{\code{weights}}{the weights of the limiting chi-bar-square distribution}
#' \item{\code{pval}}{the p-value of the test}
#'
#' @author Charlotte Baey <\email{charlotte.baey@univ-lille.fr}>
#'
#' @references Baey C, Cournède P-H, Kuhn E, 2019. Asymptotic distribution of likelihood ratio test
#' statistics for variance components in nonlinear mixed effects models. \emph{Computational
#' Statistics and Data Analysis} (in press).
#'
#' Silvapulle  MJ, Sen PK, 2011. Constrained statistical inference: order, inequality and shape constraints.
#' @export varTest
varTest <- function(m1,m0,control = list(N=5000)){ # ajouter paramètres pour option de calcul

  cat("Variance components testing in mixed-effects models\n")

  # Identify the packages from which m0 and m1 come from
  cl0 <- class(m0)
  cl1 <- class(m1)
  if (cl0[1] != cl1[1]) stop("Error: models m0 and m1 should be fitted with the same function from the same package")
  ## NB : HANDLE THE CASE WHERE ALL THE RANDOM EFFECTS ARE TESTED (m0 would be fitted using lm or gls, ...)

  linmodel <- ("lme" %in% cl1 || cl1 %in% "lmerMod")
  if (cl1[1] %in% "nlmerMod") stop("Error: the Fisher information matrix is not available for nonlinear mixed effect models fitted with nlmer() of package lme4. Please use nlme or saemix packages.")
  pkg <- pckName(m1)

  if (pkg=="nlme") msdata <- modelStructnlme(m1,m0)
  if (pkg=="lme4") msdata <- modelStructlme4(m1,m0,linmodel)
  if (pkg=="saemix") msdata <- modelStructsaemix(m1,m0)

  cat("(models fitted using the",pkg,"package)\n\n")

  # LRT and Fisher Information Matrix
  if (pkg=="nlme") lrt <- -2*(m0$logLik - m1$logLik)
  if (pkg=="lme4") lrt <- -2*(stats::logLik(m0) - stats::logLik(m1))
  if (pkg=="saemix") lrt <- -2*(m0@results@ll.is - m1@results@ll.is)

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
  cbs <- methods::new("chiBarObject")
  cbs@orthan <- msdata$structGamma$diag

  if (msdata$structGamma$diag || msdata$structGamma$full){
    # reorder FIM so that the tested variances are at the end, before the residual covariance structure
    neworder <- as.numeric(rownames(msdata$detailStruct[order(msdata$detailStruct$tested,msdata$detailStruct$names),]))
    invfim <- invfim[c(neworder,nrow(invfim)-msdata$dims$dimSigma+1),c(neworder,nrow(invfim)-msdata$dims$dimSigma+1)]

    if (msdata$structGamma$diag){
      cbs@dims <- list(dimBeta=msdata$dims$nbFE1,
                       dimGamma=list(dim0=msdata$dims$nbRE0,
                                     dimR=0,
                                     dimSplus=msdata$dims$nbRE1-msdata$dims$nbRE0),
                       dimSigma=msdata$dims$dimSigma)
    }else{
      cbs@dims <- list(dimBeta=msdata$dims$nbFE1,
                       dimGamma=list(dim0=msdata$dims$nbRE0*(msdata$dims$nbRE0+1)/2,
                                     dimR=(msdata$dims$nbRE0)*(msdata$dims$nbRE1-msdata$dims$nbRE0),
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

    cbs@dims <- list(dimBeta=msdata$dims$nbFE1,
                     dimGamma=list(dim0=dim0,
                                   dimR=dimR,
                                   dimSplus=dimSplus),
                     dimSigma=msdata$dims$dimSigma)
  }
  cbs@V <- invfim
  cbs@invV <- chol2inv(chol(invfim))

  if (length(msdata$nameVarTested)==1){
    cat(paste("Testing that the variance of",msdata$nameVarTested,"is null\n"))
  }else{
    cat(paste("Testing that the variances of",paste(msdata$nameVarTested,sep="",collapse = " and "),"are null\n"))
  }


  ## TO COMPLETE
  if (pkg=="lme4"){
    cat(paste("model under H1:",deparse(formula(m1),width.cutoff=500),"\n"))
    cat(paste("model under H0:",deparse(formula(m0),width.cutoff=500),"\n"))
  }
  if (pkg=="nlme"){
    text1 <- deparse(m1$call$random,width.cutoff=500)
    text0 <- deparse(m0$call$random,width.cutoff=500)
    if (!is.null(m1$call$model)){
      cat(paste("model under H1:",deparse(m1$call$model,width.cutoff=500),"(model),",gsub(".*formula = structure[(]list*|, class*.*$", "", text1),"(random effects)\n"))
      cat(paste("model under H0:",deparse(m0$call$model,width.cutoff=500),"(model),",gsub(".*formula = structure[(]list*|, class*.*$", "", text0),"(random effects)\n"))
    }else{
      cat(paste("model under H1:",deparse(m1$call$fixed,width.cutoff=500),"(fixed effects)",",",gsub(".*formula = structure[(]list*|, class*.*$", "", text1),"(random effects)\n"))
      cat(paste("model under H0:",deparse(m0$call$fixed,width.cutoff=500),"(fixed effects)",",",gsub(".*formula = structure[(]list*|, class*.*$", "", text0),"(random effects)\n"))
    }
  }

  # Compute chi-bar-square weights
  wcbs <- weightsChiBarSquare(cbs,control)
  pvalue <- sum(wcbs$weights * stats::pchisq(lrt,df=wcbs$df,lower.tail = F))

  # print results
  cat("\nLikelihood ratio test statistics: \n LRT = ",lrt,
      "\n\nLimiting distribution: \n",
      " mixture of",length(wcbs$df),"chi-bar-square distributions",
      "with degrees of freedom",paste(wcbs$df,sep="",collapse = ", "),"\n",
      " associated weights",paste(round(wcbs$w,3),sep="",collapse = ", "),
      "\n\np-value:",pvalue)

  #return(list(lrt=lrt,ddl=wcbs$df,weights=wcbs$w,pvalue))
}
