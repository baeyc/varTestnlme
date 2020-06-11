#' Extraction of the Fisher Information Matrix for nlme package
#'
#' Extraction of the Fisher Information Matix for variance components fitted with nlme using Delta method
#'
#' This function extract the FIM computed by the nlme for the transformed variance components, and uses
#' the Delta method to compute the FIM for the natural variance components (i.e. variances and covariances)
#'
#' @name extractFIMnlme
#' @aliases FIM nlme
#'
#' @param m a model fitted using nlme
#' @param struct a string giving the structure of the covariance matrix: either \code{diag} for a diagonal
#' matrix, \code{blockDiag} for a block diagonal matrix of \code{full} for a matrix with non-zero components
#' @return the FIM matrix for the variance components
#' @author Charlotte Baey <\email{charlotte.baey@univ-lille.fr}>
#'
#' @export extractFIMnlme
#' @importFrom stats as.formula
extractFIMnlme <- function(m,struct){
  nameGroup <- names(m$groups)
  namesFE <- names(m$coefficients$fixed)
  vc <- eval(parse(text=paste("m$modelStruct$reStruct$",nameGroup,sep="")))
  namesRE <- attr(vc,"Dimnames")[[1]]
  
  apVarNu <- m$apVar # covariance matrix of TRANSFORMED variance components
  if (is.character(apVarNu)) stop(paste0("Error in nlme package for apVar in m1: ",apVarNu,"\n Try to re-run with fim='compute'"))
  # variances are log-transformed in nu and covariances are logit-transformed (see Pinheiro and Bates)
  meanNu <- attr(apVarNu,"Pars")
  nre <- length(namesRE)
  if (struct == "diag"){
    formulaDeltaMethod <- paste0("~exp(x",1:(nre+1),")") # nb of random effects + residual
    vecAsFor <- Vectorize(as.formula,"object")
    apVarTheta <- msm::deltamethod(lapply(1:length(formulaDeltaMethod),FUN = function(i){as.formula(formulaDeltaMethod[i])}), mean=meanNu, cov=apVarNu, ses=F)
    colnames(apVarTheta) <- rownames(apVarTheta) <- c(paste0("sd(",namesRE,")"),"residual")
  }else if (struct == "full"){
    # loop over columns of apVar -> first line is diagonal elements and the others are covariance param
    formulaDeltaMethod <- c(paste0("~exp(x",1:nre,")"))
    nameParams <- paste0("sd(",namesRE,")")
    for (k in (2:nre)){
      minInd <- (k-1)*nre - sum(0:max(0,k-2)) + 1
      maxInd <- minInd + (nre - (k-1)) - 1
      formulaDeltaMethod <- c(formulaDeltaMethod,paste0("~sqrt(exp(x",k-1,")*exp(x",k:nre,"))*(exp(x",minInd:maxInd,")-1)/(exp(x",minInd:maxInd,")+1)"))
      nameParams <- c(nameParams,paste0("cor(",namesRE[k-1],",",namesRE[k:nre],")"))
    }
    formulaDeltaMethod <- c(formulaDeltaMethod,paste0("~exp(x",length(meanNu),")"))
    nameParams <- c(nameParams,"residual")
    vecAsFor <- Vectorize(as.formula,"object")
    apVarTheta <- msm::deltamethod(lapply(1:nrow(apVarNu),FUN = function(i){as.formula(formulaDeltaMethod[i])}), mean=meanNu, cov=apVarNu, ses=F)
    colnames(apVarTheta) <- rownames(apVarTheta) <- nameParams
  }else{ # block diag
    sizeBlocks <- attr(vc,"plen")
    formulaDeltaMethod <- character()
    nameParams <- character()
    indb <- 0
    for (b in 1:length(sizeBlocks)){
      nre <- floor(sqrt(2*sizeBlocks[b]))
      nameREinblock <- attr(vc[[b]],"Dimnames")[[1]]
      formulaDeltaMethod <- c(formulaDeltaMethod,paste0("~exp(x",(indb + 1):(indb + nre),")"))
      nameParams <- c(nameParams,paste0("sd(",nameREinblock,")"))
      if (nre > 1){
        for (k in (2:nre)){
          minInd <- indb + (k-1)*nre - sum(0:max(0,k-2)) + 1
          maxInd <- indb + minInd + (nre - (k-1)) - 1
          formulaDeltaMethod <- c(formulaDeltaMethod,paste0("~sqrt(exp(x",k-1,")*exp(x",k:nre,"))*(exp(x",minInd:maxInd,")-1)/(exp(x",minInd:maxInd,")+1)"))
          nameParams <- c(nameParams,paste0("cor(",nameREinblock[k-1],",",nameREinblock[k:nre],")"))
        }
      }
      indb <- indb + length(formulaDeltaMethod)
    }
    formulaDeltaMethod <- c(formulaDeltaMethod,paste0("~exp(x",length(meanNu),")"))
    nameParams <- c(nameParams,"residual")
    vecAsFor <- Vectorize(as.formula,"object")
    apVarTheta <- msm::deltamethod(lapply(1:length(formulaDeltaMethod),FUN = function(i){as.formula(formulaDeltaMethod[i])}), mean=meanNu, cov=apVarNu, ses=F)
  }
  return(apVarTheta)
}