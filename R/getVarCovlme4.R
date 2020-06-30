#' @name getVarCovlme4
#' @rdname getVarCovlme4
#'
#' @title Extract covariance matrix 
#'
#' @description Extract covariance matrix of the random effects for a model fitted with lme4.
#'
#' @param m a fit from lme4 package (either linear or nonlinear)
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