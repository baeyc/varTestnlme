#' Chi-bar-square weights computation
#'
#' Computation of the chi-bar-square weigths.
#'
#' The function computes an approximation of the weights of the chi-bar-square distribution
#' \eqn{\bar{\chi}^2(I,C)} arising as the limiting distribution of the likelihood ratio test
#' statistics under the null hypothesis. More details can be found in the references listed below
#'
#' @name weightsChiBarSquare
#' 
#' @param cbs an object of class \code{\link{chiBarSquareObject}}, containing the parameters of the chi-bar-square distribution
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
#' @importFrom foreach %dopar%
#' @importFrom stats na.omit
weightsChiBarSquare <- function(cbs,control){

  # global variables
  df <- i <- NULL
  
  # Initialize vector of weights
  w <- rep(0,length(cbs@df))

  # define some quantities : b=nb of fixed effects, r=nb of variances being tested equal to 0, ds=dim of residual, p=nb of variance components
  b <- sum(unlist(cbs@dims$dimBeta)) # size of fixed effects
  rf <- cbs@dims$dimBeta$dimR # size of tested fixed effects
  r <- sum(cbs@dims$dimGamma$dimSplus) # number of variances tested
  ds <- cbs@dims$dimSigma # size of residual
  #p <- ncol(cbs@V)-b-ds #
  #p <- cbs@dims$dimGamma$dim0 + cbs@dims$dimGamma$dimR + sum(cbs@dims$dimGamma$dimSplus)

  # Initialize vector of simulated chi-bar-square
  chibarsquare <- numeric()

  # Exact weights when testing that r variances are null with r<=3, and for the diagonal case
  # Otherwise, simulation of weights
  if(cbs@orthan){
    R <- cbind(matrix(0,ncol=b+cbs@dims$dimGamma$dim0,nrow=r),
               diag(r),
               matrix(0,nrow=r,ncol=cbs@dims$dimSigma))
    W <- R%*%cbs@V%*%t(R)
    invW <- chol2inv(chol(W))

    if (r == 1){
      w <- c(0.5,0.5)
      sdw <- rep(0,2)
    }else if (r == 2){
      C <- stats::cov2cor(W)
      w[1] <- acos(C[1,2])/(2*pi)
      w[2] <- 0.5
      w[3] <- 0.5-w[1]
      sdw <- rep(0,3)
    }else if (r == 3){
      C <- stats::cov2cor(W)
      pC <- corpcor::cor2pcor(C)
      w[4] <- (3*pi - acos(pC[1,2]) - acos(pC[1,3]) - acos(pC[2,3]))/(4*pi)
      w[3] <- (2*pi - acos(C[1,2]) - acos(C[1,3]) - acos(C[2,3]))/(4*pi)
      w[2] <- 0.5 - w[4]
      w[1] <- 0.5 - w[3]
      sdw <- rep(0,4)
    }else{
      message("Simulating chi-bar-square weights ...")
      Z <- mvtnorm::rmvnorm(control$M,mean=rep(0,nrow(W)),sigma=W)
      projZ <- t(sapply(1:control$M,FUN = function(i){
        quadprog::solve.QP(invW, Z[i,]%*%invW, diag(ncol(Z)), rep(0,ncol(Z)), meq=0, factorized=FALSE)$solution}))
      nbPosComp <- apply(projZ,1,FUN = function(x){length(x[x>1e-6])})
      w <- summary(as.factor(nbPosComp),maxsum=length(df))/control$M
      sdw <- sqrt(w*(1-w)/control$M)
    }
  }else if (r == 1){
    w <- c(0.5,0.5)
    sdw <- rep(0,2)
  }else{
    invI <- cbs@V
    I <- cbs@invV
    Z <- mvtnorm::rmvnorm(control$M,mean=rep(0,nrow(invI)),sigma=invI)
    chibarsquare <- rep(0,control$M)

    cat("\nSimulating chi-bar-square weights ...\n")
    no_cores <- max(1,parallel::detectCores() - 1)
    # Initiate cluster
    doParallel::registerDoParallel(no_cores)

    chibarsquare <- foreach::foreach(i=1:control$M, 
                                     .packages='varTestnlme', 
                                     .combine=c) %dopar% {
          x0 <- Z[i,]
          constantes <- list(Z = Z[i,], invV = I, cbs = cbs)
          projZ <- alabama::auglag(par=x0,
                                   fn = objFunction,
                                   gr = gradObjFunction,
                                   hin = ineqCstr,
                                   heq = eqCstr,
                                   hin.jac = jacobianIneqCstr,
                                   heq.jac = jacobianEqCstr,
                                   cst=constantes,
                                   control.outer=list(trace=F,method="BFGS",eps=1e-5))
          # due to the tolerance threshold in auglag, we end up with equality constraints verified up to 1e-5
          # we set those values to 0 otherwise it results in incorrect projected Z
          projZzero <- projZ$par
          projZzero[abs(projZzero)<1e-5] <- 0
          Z[i,]%*%I%*%Z[i,] - objFunction(projZzero,constantes)
        }
    chibarsquare <- chibarsquare[chibarsquare>-1e-04]
    doParallel::stopImplicitCluster()
    
    # Approximation of chi-bar-square weights
    empQuant <- seq(0.001,1,0.001)
    w_covw <- lapply(empQuant, FUN = function(q){try(approxWeights(chibarsquare,cbs@df,q))})
    w <- na.omit(t(sapply(1:length(empQuant), FUN = function(i){if (is.null(names(w_covw[[i]]))) {NA*df} else {w_covw[[i]]$w}})))
    covw <- na.omit(t(sapply(1:length(empQuant), FUN = function(i){if (is.null(names(w_covw[[i]]))) {NA*df} else {sqrt(diag(w_covw[[i]]$covw))}})))

    # Take into account the constraints on the weights : they should be between 0 and 1/2
    minW <- apply(w,1,min)
    maxW <- apply(w,1,max)
    admissibleWeights <- w[minW>0 & maxW<0.5,]
    sdAdmissibleWeights <- covw[minW>0 & maxW<0.5,]

    # NCOL and NROW can be used on vector
    eps <- 0
    while (min(NROW(admissibleWeights),NCOL(admissibleWeights))==0){
      # if no combination is admissible, relax the constraints
      admissibleWeights <- w[minW>-(0.05+eps) & maxW<(0.55+eps),]
      sdAdmissibleWeights <- covw[minW>-(0.05+eps) & maxW<(0.55+eps),]
      eps <- eps+0.01
    }

    if (min(NROW(admissibleWeights),NCOL(admissibleWeights)) > 1){
      # if several weights are admissible, we keep those whose variances are minimal
      maxvar <- apply(sdAdmissibleWeights,1,max)
      #meanvar <- apply(sdAdmissibleWeights,1,mean)
      admissibleWeights <- admissibleWeights[which.min(maxvar),]
      sdAdmissibleWeights <- sdAdmissibleWeights[which.min(maxvar),]
    }

    w <- admissibleWeights
    sdw <- sdAdmissibleWeights
  }

  return(list(df=cbs@df,weights=w,sdWeights=sdw,randomCBS=chibarsquare))
}

