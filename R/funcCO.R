##' @name funcCO
##' @rdname funcCO
##'
##' @title Internal functions for constrained minimization
##'
##' @description Groups of functions used for the constrained minimization problem arising in the computation of the
##' likelihood ratio test statistics.
##'
##' @param x A vector
##' @param cst A list of constants to be passed to the optimisation function
##' @return value of the objective function, its gradient, and the set of inequality and euqality constraints
##'
NULL

##' @rdname funcCO
##' @export
objFunction <- function(x,cst){
  return(t(cst$Z-x)%*%cst$invV%*%(cst$Z-x))
}


##' @rdname funcCO
##' @export
gradObjFunction <- function(x,cst){
  return(-2*t(cst$Z-x)%*%cst$invV)
}

# Function creating a symmetric matrix from its unique elements stored in a vector.
#
# @param v A vector
# @return A matrix containing the elements of \code{v}.
# @example
# symMatrixFromVect(c(1,2,3,4,5,6))
#
symMatrixFromVect <- function(v){
  n <- length(v)

  if (n==1){
    return(v)
  }else{
    p <- floor(sqrt(2*n))# nb of rows/columns

    m <- matrix(0,nrow=p,ncol=p)
    m[lower.tri(m,diag=TRUE)] <- v

    return(m+t(m)-diag(diag(m)))
  }
}


##' @rdname funcCO
##' @export
ineqCstr <- function(x,cst){
  r <- cst$cbs@dims$dimGamma$dimSplus
  n0R <- cst$cbs@dims$dimGamma$dim0 + cst$cbs@dims$dimGamma$dimR
  dimMats <- r*(r+1)/2
  n <- length(x)
  ds <- cst$cbs@dims$dimSigma

  if (sum(r) == 1){
    constr <- x[n-ds]
  }else{
    constr <- numeric()
    for (i in 1:length(r)){
      m <- symMatrixFromVect(x[(n0R+1):(n0R+dimMats[i])])
      n0R <- n0R + dimMats[i]
      if (dimMats[i]==1){
        constr <- c(constr,m)
      }else{
        constr <- c(constr,det(m))
      }
    }
  }

  return(constr)
}



##' @rdname funcCO
##' @export
ineqCstrDiag <- function(x,cst){
  r <- cst$cbs@dims$dimGamma$dimSplus
  n0R <- cst$cbs@dims$dimGamma$dim0 + cst$cbs@dims$dimGamma$dimR
  n <- length(x)
  ds <- cst$cbs@dims$dimSigma

  if (r == 1){
    constr <- x[n-ds]
  }else{
    constr <- x[(n0R+1):(n0R+r)]
  }

  return(constr)
}


##' @rdname funcCO
##' @export
eqCstr <- function(x,cst){
  n0 <- cst$cbs@dims$dimGamma$dim0
  n <- length(x)
  ds <- cst$cbs@dims$dimSigma

  if (n0==0)
  {
   constr <- x[(n-ds+1):n]
  }else{
   constr <- x[c(1:n0,(n-ds+1):n)]
  }
  return(constr)
}




