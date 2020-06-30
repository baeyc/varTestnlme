#' @name getVarCovnlme
#' @rdname getVarCovnlme
#'
#' @title Extract covariance matrix 
#'
#' @description Extract covariance matrix of the random effects for a model fitted with nlme.
#'
#' @param m a fit from nlme package (either linear or nonlinear)
getVarCovnlme <- function(m){
  
  pkg <- class(m)
  if (! ("nlme" %in% pkg)){
    v <- nlme::getVarCov(m)
    v <- matrix(as.numeric(v),ncol=ncol(v),nrow=nrow(v))
  }else{
    varcorr <- nlme::VarCorr(m)
    coln <- colnames(varcorr)
    if ("Corr" %in% coln){
      stdev <- as.numeric(varcorr[-nrow(varcorr),2])
      corrmat <- diag(1,length(stdev))
      # NA appear during the conversion from characters to numbers, due to the fact that only the lower
      # triangular of the correlation matrix is printed -> NA can be removed safely and correspond to blanks
      corrmat[lower.tri(corrmat,diag=FALSE)] <- na.omit(as.numeric(varcorr[-c(1,nrow(varcorr)),-c(1:2)]))
      corrmat[upper.tri(corrmat,diag=FALSE)] <- na.omit(as.numeric(varcorr[-c(1,nrow(varcorr)),-c(1:2)]))
      v <-  diag(stdev)%*%corrmat%*%diag(stdev)
      colnames(v) <- rownames(varcorr)[1:length(stdev)]
      rownames(v) <- rownames(varcorr)[1:length(stdev)]
    }else{
      v <- diag(as.numeric(varcorr[-nrow(varcorr),1]),nrow=nrow(varcorr)-1)
    }
    
    v <- matrix(as.numeric(v),ncol=ncol(v),nrow=nrow(v))
    colnames(v) <- rownames(varcorr)[1:ncol(v)]
    rownames(v) <- rownames(varcorr)[1:nrow(v)]
  }
  return(v)
}


getVarCovlme4 <- function(m){
  varcorr <- as.data.frame(lme4::VarCorr(m))
  varcorr <- varcorr[varcorr$grp!="Residual",]
  indCov <- which(!is.na(varcorr$var2))
  if (length(indCov)>0){
    stdev <- varcorr$sdcor[-indCov]
    corrmat <- diag(1,length(stdev))
    corrmat[lower.tri(corrmat,diag=FALSE)] <- varcorr$sdcor[indCov]
    corrmat[upper.tri(corrmat,diag=FALSE)] <- varcorr$sdcor[indCov]
    v <-  diag(stdev)%*%corrmat%*%diag(stdev)
    colnames(v) <- varcorr$var1[-indCov]
    rownames(v) <- varcorr$var1[-indCov]
  }else{
    v <- diag(varcorr$vcov,nrow=nrow(varcorr))
    colnames(v) <- varcorr$var1
    rownames(v) <- varcorr$var1
  }
  
  return(v)
}