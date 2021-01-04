#' Extract package name from a fitted mixed-effects model
#' 
#' @param m a model with random effects fitted with \code{nlme}, \code{lme4} or \code{saemix}
#' @return a string giving the name of the package
pckName <- function(m){
  cl <- class(m)
  if ("lme" %in% cl) {
    pkg <- "nlme"
  }else if(max(c("lmerMod","glmerMod","nlmerMod") %in% cl)){
    pkg <- "lme4"
  }else if("SaemixObject" %in% cl){
    pkg <- "saemix"
  }else{
    pkg <- cl
  }
  return(pkg)
}
