#' Extract package name from a fitted mixed-effects model
#' 
#' @param m a model with random effects fitted with \code{nlme}, \code{lme4} or \code{saemix}
#' @return a string giving the name of the package
pckName <- function(m){
  if (inherits(m,"lme")) {
    pkg <- "nlme"
  }else if(inherits(m,c("lmerMod","glmerMod","nlmerMod","lmerModLmerTest"))){
    pkg <- "lme4"
  }else if(inherits(m,"SaemixObject")){
    pkg <- "saemix"
  }else{
    pkg <- class(m)
  }
  return(pkg)
}
