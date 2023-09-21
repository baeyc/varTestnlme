#' @title Chi-bar-square degrees of freedom computation
#'
#' @description Computation of the degrees of freedom of the chi-bar-square 
#'
#' @param msdata a list containing the structure of the model and data, as an output from 
#' \code{extractStruct.<package_name>} functions
#' @return a list containing the vector of the degrees of freedom of the chi-bar-square and the dimensions
#' of the cone of the chi-bar-square distribution
#' 
#' @name dfChiBarSquare
# #' @export dfChiBarSquare
dfChiBarSquare <- function(msdata){
  
  orthan <- (msdata$struct == "diag")
  
  if (msdata$struct == "diag"){
    dims <- list(dimBeta=list(dim0=msdata$dims$nbFE0,
                              dimR=msdata$dims$nbFE1-msdata$dims$nbFE0),
                 dimGamma=list(dim0=msdata$dims$nbRE0,
                               dimR=0,
                               dimSplus=msdata$dims$nbRE1-msdata$dims$nbRE0),
                 dimSigma=msdata$dims$dimSigma)
  }else if (msdata$struct == "full"){
    nbCovTestedAlone <- (msdata$dims$nbRE0==msdata$dims$nbRE1)*(msdata$dims$nbCov1-msdata$dims$nbCov0) # nb of covariances tested without the corresponding variances being tested
    dims <- list(dimBeta=list(dim0=msdata$dims$nbFE0,
                              dimR=msdata$dims$nbFE1-msdata$dims$nbFE0),
                 dimGamma=list(dim0=msdata$dims$nbRE0*(msdata$dims$nbRE0+1)/2,
                               dimR=(msdata$dims$nbRE0)*(msdata$dims$nbRE1-msdata$dims$nbRE0)+nbCovTestedAlone,
                               dimSplus=msdata$dims$nbRE1-msdata$dims$nbRE0),
                 dimSigma=msdata$dims$dimSigma)
  }else{
    dd <- msdata$detailStruct
    dim0 <- nrow(dd[!dd$tested,])-msdata$dims$nbFE1
    #dimRplus <- length(dd$names[dd$tested & dd$diagBlock])
    dimR <- length(dd$names[dd$covInBlock!=1 & dd$type=="co" & dd$tested])
    ddSp <- dd[dd$tested & dd$covInBlock,]
    dimSplus <- floor(sqrt(2*as.vector(table(ddSp$block))))
    rm(ddSp)
    
    dims <- list(dimBeta=list(dim0=msdata$dims$nbFE0,
                              dimR=msdata$dims$nbFE1-msdata$dims$nbFE0),
                 dimGamma=list(dim0=dim0,
                               dimR=dimR,
                               dimSplus=dimSplus),
                 dimSigma=msdata$dims$dimSigma)
  }
  
  # Identify the components of the mixture
  # dimLSincluded : dimension of the buggest linear space included in the cone
  # dimLScontaining : dimension of the smallest linear space containing the cone
  q <- sum(unlist(dims)) # total dimension
  dimLSincluded <- dims$dimBeta$dimR + dims$dimGamma$dimR # weights of chi-bar with df=0, ..., dimLSincluded-1 are null
  if (orthan) dimLScontaining <- dims$dimBeta$dimR + dims$dimGamma$dimR + sum(dims$dimGamma$dimSplus) # weights of chi-bar with df=dimLScontaining+1, ..., q are null
  if (!orthan) dimLScontaining <- dims$dimBeta$dimR + dims$dimGamma$dimR + sum(dims$dimGamma$dimSplus*(dims$dimGamma$dimSplus+1)/2)
  
  return(list(df=seq(dimLSincluded,dimLScontaining,1),dimsCone=dims))
}

#' @title Monte Carlo approximation of chi-bar-square weights
#'
#' @description The function provides a method to approximate the weights of the mixture components, 
#' when the number of components is known as well as the degrees of freedom of each chi-square distribution 
#' in the mixture, and given a vector of simulated values from the target \eqn{\bar{\chi}^2(V,C)} 
#' distribution. Note that the estimation is based on (pseudo)-random Monte Carlo samples. For reproducible
#' results, one should fix the seed of the (pseudo)-random number generator.
#' 
#' @name weightsChiBarSquare
# #' @export weightsChiBarSquare
#' 
#' @param df a vector with the degrees of freedom of the chi-square components of the chi-bar-square distribution
#' @param V a positive semi-definite matrix 
#' @param dimsCone a list with the dimensions of the cone C, expressed on the parameter space scale
#' @param orthan a boolean specifying whether the cone is an orthan
#' @param control (optional) a list of control options for the computation of the chi-bar-weights, containing 
#' two elements: \code{parallel} a boolean indicating whether computation should be done in parallel (FALSE 
#' by default), \code{nb_cores} the number of cores for parallel computing (if \code{parallel=TRUE} but no value is given
#' for \code{nb_cores}, it is set to number of detected cores minus 1), and \code{M} the Monte Carlo sample
#' size for the computation of the weights.
#' 
#' @return A list containing the estimated weights, the standard deviations of the estimated weights and the
#' random sample of \code{M} realizations from the chi-bar-square distribution
#'
#' @importFrom foreach %dopar%
#' @importFrom stats na.omit
weightsChiBarSquare <- function(df,V,dimsCone,orthan,control){
  
  # Initialize vector of weights
  w <- rep(0,length(df))
  
  b <- sum(unlist(dimsCone$dimBeta))      # nb of fixed effects
  rf <- dimsCone$dimBeta$dimR             # nb of tested fixed effects
  r <- sum(dimsCone$dimGamma$dimSplus)    # nb of variances tested
  ds <- dimsCone$dimSigma                 # size of residual
  
  # Initialize vector of simulated chi-bar-square
  chibarsquare <- numeric()
  
  if (length(df) == 2){
    w <- c(0.5,0.5)
    sdw <- rep(0,2)
  }else{
    if(orthan){
      R <- cbind(matrix(0,ncol=b+dimsCone$dimGamma$dim0,nrow=r),
                 diag(r),
                 matrix(0,nrow=r,ncol=dimsCone$dimSigma))
      W <- R%*%V%*%t(R)
      invW <- chol2inv(chol(W))
      
      if (r == 2){
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
    }else{
      invV <- chol2inv(chol(V))
      Z <- mvtnorm::rmvnorm(control$M,mean=rep(0,nrow(V)),sigma=V)
      chibarsquare <- rep(0,control$M)
      
      message("Simulating chi-bar-square weights ...")
      nb_cores <- ifelse(control$parallel,max(1,parallel::detectCores() - 1),1)
      
      # Initiate cluster
      doParallel::registerDoParallel(nb_cores)
      
      i <- 0
      chibarsquare <- foreach::foreach(i=1:control$M, 
                                       .packages='varTestnlme', 
                                       .combine=c) %dopar% {
                                         x0 <- Z[i,]
                                         constantes <- list(Z = Z[i,], invV = invV, orthan = orthan, dimsCone = dimsCone)
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
                                         Z[i,]%*%invV%*%Z[i,] - objFunction(projZzero,constantes)
                                       }
      chibarsquare <- chibarsquare[chibarsquare>-1e-04]
      doParallel::stopImplicitCluster()
      
      # Approximation of chi-bar-square weights
      empQuant <- seq(0.001,1,0.001)
      w_covw <- lapply(empQuant, FUN = function(q){try(approxWeights(chibarsquare,df,q))})
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
      
      w <- apply(admissibleWeights,2,mean)
      sdw <- sqrt(apply(sdAdmissibleWeights,2,mean)/nrow(admissibleWeights))
      
    }
  }
  
  return(list(weights=w,sdWeights=sdw,randomCBS=chibarsquare))
}

#' @title Monte Carlo approximation of chi-bar-square weights
#'
#' @description The chi-bar-square distribution \eqn{\bar{\chi}^2(I,C)} is a mixture of chi-square distributions. The function provides
#' a method to approximate the weights of the mixture components, when the number of components is known as well as the
#' degrees of freedom of each chi-square distribution in the mixture, and given a vector of simulated values from the target
#' \eqn{\bar{\chi}^2(I,C)} distribution. Note that the estimation is based on (pseudo)-random Monte Carlo samples. For reproducible
#' results, one should fix the seed of the (pseudo)-random number generator.
#' 
#' @details Let us assume that there are \eqn{p} components in the mixture, with degrees of
#' freedom between \eqn{n_1} and \eqn{n_p}. By definition of a mixture distribution, we have :
#' \deqn{ P(\bar{\chi}^2(I,C) \leq c) = \sum_{i=n_1}^{n_p} w_i P(\chi^2_{i} \leq c)}
#' Choosing \eqn{p-2} values \eqn{c_1, \dots, c_{p-2}}, the function will generate a system of \eqn{p-2} equations
#' according to the above relationship, and add two additional relationships stating that the sum of all the weights is
#' equal to 1, and that the sum of odd weights and of even weights is equal to 1/2, so that we end up with a system a \eqn{p}
#' equations with \eqn{p} variables.
#'
#' @name approxWeights
#' @export approxWeights
#'
#' @param x a vector of i.i.d. random realizations of the target chi-bar-square distribution
#' @param df a vector containing the degrees of freedom of the chi-squared components
#' @param q the empirical quantile of \code{x} used to choose the \eqn{p-2} values \eqn{c_1, \dots, c_{p-2}} (see Details)
#' @return A vector containing the estimated weights, as well as their covariance matrix.
#' @author Charlotte Baey <\email{charlotte.baey@univ-lille.fr}>
#'
#' @importFrom stats quantile pchisq qchisq cov
approxWeights <- function(x,df,q){
  maxcbs <- max(0,quantile(x,q))
  epsilon <- pchisq(maxcbs,df=max(df))
  c <- numeric()
  
  if(length(df)>2){
    c <- sapply(df[-c(1,2)],FUN = function(i){qchisq(epsilon,i,lower.tail = T)})
  }else{
    c <- sapply(df,FUN = function(i){qchisq(epsilon,i,lower.tail = T)})
  }
  
  aij <- sapply(df,FUN = function(i){pchisq(c,df=i)})
  phatcbs <- sapply(1:length(c),FUN = function(i){x<=c[i]})
  pj <- apply(phatcbs,2,mean)
  covpj <- cov(phatcbs)/length(x)
  
  # We add constraints on the weights: sum of weights is equal to 1 and sums of even and of odd weights are equal to 1/2
  aij <- rbind(aij,rep(1,length(df)),rep(c(1,0),length.out=length(df)))
  pj <- c(pj,1,0.5)
  
  invAij <- solve(aij)
  w <- as.vector(invAij%*%pj)
  
  covpj2 <- matrix(0,nrow=length(df),ncol=length(df))
  covpj2[1:(length(df)-2),1:(length(df)-2)] <- covpj
  covw <- invAij %*% covpj2 %*% t(invAij)
  
  return(list(w=w,covw=covw))
}
