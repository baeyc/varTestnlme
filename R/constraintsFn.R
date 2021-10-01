#' Internal functions for constrained minimization
#'
#' Groups of functions used for the constrained minimization problem arising in the computation of the
#' likelihood ratio test statistics. 
#'
#' @describeIn objFunction objective function to be optimized
#'
#' @param x A vector
#' @param cst A list of constants to be passed to the optimisation function
#' @return value of the objective function, its gradient, and the set of inequality and equality constraints
#'
objFunction <- function(x,cst){
  return(t(cst$Z-x)%*%cst$invV%*%(cst$Z-x))
}


#' @describeIn objFunction gradient of the objective function
gradObjFunction <- function(x,cst){
  return(-2*t(cst$Z-x)%*%cst$invV)
}

#' @describeIn objFunction function creating a symmetric matrix from its unique elements stored in a vector
symMatrixFromVect <- function(x){
  n <- length(x)

  if (n==1){
    return(x)
  }else{
    p <- floor(sqrt(2*n)) 

    m <- matrix(0,nrow=p,ncol=p)
    m[lower.tri(m,diag=TRUE)] <- x

    return(m+t(m)-diag(diag(m)))
  }
}


#' @describeIn objFunction set of inequality constraints
ineqCstr <- function(x,cst){
  r <- cst$dimsCone$dimGamma$dimSplus # nb of variances tested in each block
  # n0R : nb of components in the cone corresponding to {0} or R
  n0R <- cst$dimsCone$dimBeta$dim0 + cst$dimsCone$dimBeta$dimR + cst$dimsCone$dimGamma$dim0 + cst$dimsCone$dimGamma$dimR
  dimMats <- r*(r+1)/2 # dimensions of the blocks of variances tested
  n <- length(x) # total dimension of parameter space
  ds <- cst$dimsCone$dimSigma # dimension of residual parameter

  if (sum(r) == 1){
    constr <- x[n-ds]
  }else{
    if (cst$orthan){
      constr <- x[(n0R+1):(n0R+sum(r))]
    }else{
      constr <- numeric()
      for (i in 1:length(r)){
        m <- symMatrixFromVect(x[(n0R+1):(n0R+dimMats[i])]) # reconstruction of the tested block from its unique elements in a vectorized form
        n0R <- n0R + dimMats[i]
        if (dimMats[i]==1){
          constr <- c(constr,m)
        }else{
          constr <- c(constr,det(m)) # add the constraint that the determinant should be positive
        }
      }
    }
  }
  return(constr)
}


#' @describeIn objFunction jacobian of the inequality constraints
jacobianIneqCstr <- function(x,cst){
  r <- cst$dimsCone$dimGamma$dimSplus
  n0R <- cst$dimsCone$dimBeta$dim0 + cst$dimsCone$dimBeta$dimR + cst$dimsCone$dimGamma$dim0 + cst$dimsCone$dimGamma$dimR # dimensions of spaces {0} and R (corresponding to non-constrained elements)
  dimMats <- r*(r+1)/2
  n <- length(x)
  ds <- cst$dimsCone$dimSigma

  # Jacobian of the inequality constrains
  jacobian <- matrix(0,nrow=length(dimMats),ncol=n)

  if (cst$orthan){
    jacobian <- matrix(0,nrow=r,ncol=n)
    jacobian[1:r,(n0R+1):(n0R+r)] <- diag(r)
  }else{
    if (sum(r) == 1){
      # If r=1 there is only one non-null term in the Jacobian
      jacobian[1,n-ds] <- 1
    }else{
      # If r>1 the gradient is computed for each block of the matrix
      for (i in 1:length(r)){
        m <- symMatrixFromVect(x[(n0R+1):(n0R+dimMats[i])])
        if (dimMats[i]==1){
          jacobian[i,n0R+1] <- 1
        }else{
          derivMatrixForm <- det(m) * (2*solve(m) - solve(m)*diag(r[i]))
          jacobian[i,(n0R+1):(n0R+dimMats[i])] <- derivMatrixForm[lower.tri(derivMatrixForm,diag = T)]
        }
        n0R <- n0R + dimMats[i]
      }
    }
  }

  return(jacobian)
}


#' @describeIn objFunction set of equality constraints
eqCstr <- function(x,cst){
  # equality constraints come from the fixed effects, the untested blocks, the untested covariances in partially tested blocks, and the residual
  nontestedFix <- cst$dimsCone$dimBeta$dim0
  nbFix <- cst$dimsCone$dimBeta$dim0 + cst$dimsCone$dimBeta$dimR
  n0 <- nbFix + cst$dimsCone$dimGamma$dim0
  n <- length(x)
  ds <- cst$dimsCone$dimSigma

  if (nontestedFix < nbFix){ # if some fixed effects are tested
    constr <- x[c(1:nontestedFix,(nbFix+1):n0,(n-ds+1):n)]
  }else{
    if (ds>=1) constr <- x[c(1:n0,(n-ds+1):n)]
    if (ds==0) constr <- x[1:n0]
  }

  return(constr)
}


#' @describeIn objFunction jacobian of the inequality constraints
jacobianEqCstr <- function(x,cst){
  n0 <- cst$dimsCone$dimBeta$dim0 + cst$dimsCone$dimBeta$dimR + cst$dimsCone$dimGamma$dim0
  n <- length(x)
  ds <- cst$dimsCone$dimSigma

  jacobian <- matrix(0,nrow=n0+ds,ncol=n)

  jacobian[1:n0,1:n0] <- diag(n0)
  if (ds>1) jacobian[(n0+1):(n0+ds),(n-ds+1):n] <- diag(ds)

  return(jacobian)
}


