#' Chi-bar-square weights computation
#'
#' Computation of the chi-bar-square weigths via Monte Carlo approximation.
#'
#' The function computes an approximation of the weights of the chi-bar-square distribution
#' \eqn{\bar{\chi}^2(I,C)} arising as the limiting distribution of the likelihood ratio test
#' statistics under the null hypothesis. More details can be found in the references listed below
#'
#' @name weightsChiBarSquare
#' @aliases weights chibarsquare limit
#'
#' @param cbs an object of class \code{\link{chibarObject}}, containing the parameters of the chi-bar-square distribution
#' @param control (optional) a list of control options for the computation of the chi-bar-weights
#' @return A list containing the degrees of freedom of the chi-bar distributions involved in the chi-bar-square, along with
#' the associated weights.
#' @author Charlotte Baey <\email{charlotte.baey@univ-lille.fr}>
#' @references Baey C, Cournède P-H, Kuhn E, 2019. Asymptotic distribution of likelihood ratio test
#' statistics for variance components in nonlinear mixed effects models. \emph{Computational
#' Statistics and Data Analysis} 135:107-122.
#'
#' Silvapulle  MJ, Sen PK, 2011. Constrained statistical inference: order, inequality and shape constraints.
#' @export weightsChiBarSquare
weightsChiBarSquare <- function(cbs,control = list(N=5000)){

  # some weights are null, depending on whether the cone includes, or is contained in a linear space
  dimLSincluded <- cbs@dims$dimGamma$dimR # weights of chi-bar with df=0, ..., dimLSincluded-1 are null
  if (cbs@orthan) dimLScontaining <- cbs@dims$dimGamma$dimR+sum(cbs@dims$dimGamma$dimSplus) # weights of chi-bar with df=dimLScontaining+1, ..., q are null
  if (!cbs@orthan) dimLScontaining <- cbs@dims$dimGamma$dimR+sum(cbs@dims$dimGamma$dimSplus*(cbs@dims$dimGamma$dimSplus+1)/2)
  dimTot <- ncol(cbs@V)
  df <- seq(dimLSincluded,dimLScontaining,1)

  w <- rep(0,length(df))

  b <- cbs@dims$dimBeta
  r <- sum(cbs@dims$dimGamma$dimSplus)
  ds <- cbs@dims$dimSigma
  p <- dimTot-b-ds
  #p <- cbs@dims$dimGamma$dim0 + cbs@dims$dimGamma$dimR + sum(cbs@dims$dimGamma$dimSplus)

  # Exact weights when testing that r variances are null with r<=3, and for the diagonal case
  # Otherwise, simulation of weights
  if(cbs@orthan){
    R <- cbind(matrix(0,ncol=cbs@dims$dimBeta+cbs@dims$dimGamma$dim0,nrow=r),
               diag(r),
               matrix(0,nrow=r,ncol=cbs@dims$dimSigma))
    W <- R%*%cbs@V%*%t(R)
    invW <- chol2inv(chol(W))
    C <- stats::cov2cor(W)

    if (r == 1){
      w <- c(0.5,0.5)
    }else if (r == 2){
      w[1] <- acos(C[1,2])/(2*pi)
      w[2] <- 0.5
      w[3] <- 0.5-w[1]
    }else if (r == 3){
      pC <- corpcor::cor2pcor(C)
      w[4] <- (3*pi - acos(pC[1,2]) - acos(pC[1,3]) - acos(pC[2,3]))/(4*pi)
      w[3] <- (2*pi - acos(C[1,2]) - acos(C[1,3]) - acos(C[2,3]))/(4*pi)
      w[2] <- 0.5 - w[4]
      w[1] <- 0.5 - w[3]
    }else{
      print("Simulating chi-bar-square weights ...")
      pb = utils::txtProgressBar(min = 0, max = control$N, initial = 0)
      Z <- mvtnorm::rmvnorm(control$N,mean=rep(0,nrow(W)),sigma=W)
      chibarsquare <- t(sapply(1:control$N,FUN = function(i){
                                           projZ <- quadprog::solve.QP(invW, Z[i,]%*%invW, diag(ncol(Z)), rep(0,ncol(Z)), meq=0, factorized=FALSE)$solution
                                           utils::setTxtProgressBar(pb,i)}))
      nbPosComp <- apply(chibarsquare,1,FUN = function(x){length(x[x>1e-6])})
      w <- summary(as.factor(nbPosComp),maxsum=length(df))/control$N
    }
  }else if (r == 1){
    w <- c(0.5,0.5)
  }else{
    invI <- cbs@V[(b+1):dimTot,(b+1):dimTot]
    I <- cbs@invV[(b+1):dimTot,(b+1):dimTot]
    Z <- mvtnorm::rmvnorm(control$N,mean=rep(0,nrow(invI)),sigma=invI)
    chibarsquare <- rep(0,control$N)

    cat("\nSimulating chi-bar-square weights ...\n")
    pb = utils::txtProgressBar(min = 0, max = control$N, initial = 0)
    if (cbs@orthan){
      chibarsquare <- t(sapply(1:control$N,FUN = function(i){
                                            x0 <- Z[i,]
                                            constantes <- list(Z=Z[i,],invV=I,cbs=cbs)
                                            projZ <- alabama::auglag(par=x0, fn = objFunction, gr = gradObjFunction, hin = ineqCstrDiag, heq = eqCstr, cst=constantes, control.outer=list(trace=F))
                                            utils::setTxtProgressBar(pb,i)
                                            return(Z[i,]%*%I%*%Z[i,] - projZ$value)}))
    }else{
      chibarsquare <- t(sapply(1:control$N,FUN = function(i){
                                            x0 <- Z[i,]
                                            constantes <- list(Z=Z[i,],invV=I,cbs=cbs)
                                            projZ <- alabama::auglag(par=x0, fn = objFunction, gr = gradObjFunction, hin = ineqCstr, heq = eqCstr, cst=constantes, control.outer=list(trace=F))
                                            utils::setTxtProgressBar(pb,i)
                                            return(projZ$value)}))
    }

    c <- seq(min(chibarsquare*1.1),max(chibarsquare)*0.8,length.out = length(df)-2)
    aij <- sapply(1:length(df),FUN = function(i){stats::pchisq(c,df=i-1)})
    pj <- sapply(1:length(c),FUN = function(i){mean(chibarsquare<=c[i])})

    # We add constraints on the weights: sum of weights is equal to 1 and sum of even and odd weights are equal to 1/2
    aij <- rbind(aij,rep(1,length(df)),rep(c(1,0),length.out=length(df)))
    pj <- c(pj,1,0.5)

    w <- solve(aij,pj)
  }

  return(list(df=df,weights=w))
}
