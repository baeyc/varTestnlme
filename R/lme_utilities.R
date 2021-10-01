# Extract FIM
extractFIM.lme <- function(m,struct){
  nameGroup <- names(m$groups)
  namesFE <- names(m$coefficients$fixed)
  vc <- eval(parse(text=paste("m$modelStruct$reStruct$",nameGroup,sep="")))
  namesRE <- attr(vc,"Dimnames")[[1]]
  
  apVarNu <- m$apVar # covariance matrix of TRANSFORMED variance components
  if (is.character(apVarNu)) stop(paste0("in nlme package for apVar in m1: ",apVarNu,". Try to re-run with fim='compute'"))
  # variances are log-transformed in nu and covariances are logit-transformed (see Pinheiro and Bates)
  meanNu <- attr(apVarNu,"Pars")
  nre <- length(namesRE)
  if (struct == "diag"){
    formulaDeltaMethod <- paste0("~exp(x",1:(nre+1),")^2") # nb of random effects + residual
    vecAsFor <- Vectorize(stats::as.formula,"object")
    apVarTheta <- msm::deltamethod(lapply(1:length(formulaDeltaMethod),FUN = function(i){stats::as.formula(formulaDeltaMethod[i])}), mean=meanNu, cov=apVarNu, ses=F)
    colnames(apVarTheta) <- rownames(apVarTheta) <- c(paste0("sd(",namesRE,")"),"residual")
  }else if (struct == "full"){
    # loop over columns of apVar -> first line is diagonal elements and the others are covariance param
    formulaDeltaMethod <- c(paste0("~exp(x",1:nre,")^2"))
    nameParams <- paste0("sd(",namesRE,")")
    for (k in (2:nre)){
      minInd <- (k-1)*nre - sum(0:max(0,k-2)) + 1
      maxInd <- minInd + (nre - (k-1)) - 1
      formulaDeltaMethod <- c(formulaDeltaMethod,paste0("~sqrt(exp(x",k-1,")*exp(x",k:nre,"))*(exp(x",minInd:maxInd,")-1)/(exp(x",minInd:maxInd,")+1)"))
      nameParams <- c(nameParams,paste0("cor(",namesRE[k-1],",",namesRE[k:nre],")"))
    }
    formulaDeltaMethod <- c(formulaDeltaMethod,paste0("~exp(x",length(meanNu),")"))
    nameParams <- c(nameParams,"residual")
    vecAsFor <- Vectorize(stats::as.formula,"object")
    apVarTheta <- msm::deltamethod(lapply(1:nrow(apVarNu),FUN = function(i){stats::as.formula(formulaDeltaMethod[i])}), mean=meanNu, cov=apVarNu, ses=F)
    colnames(apVarTheta) <- rownames(apVarTheta) <- nameParams
  }else{ # block diag
    sizeBlocks <- attr(vc,"plen")
    formulaDeltaMethod <- character()
    nameParams <- character()
    indb <- 0
    for (b in 1:length(sizeBlocks)){
      nre <- floor(sqrt(2*sizeBlocks[b]))
      nameREinblock <- attr(vc[[b]],"Dimnames")[[1]]
      formulaDeltaMethod <- c(formulaDeltaMethod,paste0("~exp(x",(indb + 1):(indb + nre),")^2"))
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
    formulaDeltaMethod <- c(formulaDeltaMethod,paste0("~exp(x",length(meanNu),")^2"))
    nameParams <- c(nameParams,"residual")
    vecAsFor <- Vectorize(stats::as.formula,"object")
    apVarTheta <- msm::deltamethod(lapply(1:length(formulaDeltaMethod),FUN = function(i){stats::as.formula(formulaDeltaMethod[i])}), mean=meanNu, cov=apVarNu, ses=F)
  }
  
  invfim <- as.matrix(Matrix::bdiag(m$varFix,apVarTheta))
  colnames(invfim) <- rownames(invfim) <- c(colnames(m$varFix),colnames(apVarTheta))
  return(invfim)
}


#' @name extractStruct.lme
#' @rdname extractStruct.lme
#'
#' @title Extract model structure
#'
#' @param m1 the fit under H1
#' @param m0 the fit under H0
#' @param randm0 a boolean indicating whether random effects are present in m0
extractStruct.lme <- function(m1,m0,randm0){
  # name of the random grouping variable
  nameRE <- names(m1$groups)
  
  # get the residual variance structure
  varStructm1 <- m1$modelStruct$varStruct
  varStructm0 <- m0$modelStruct$varStruct
  if (length(varStructm0) != length(varStructm1)){
    stop("the residual variance model should be the same under both hypotheses")
  }else if (length(varStructm0) > 0 & length(varStructm1) > 0){
    if (varStructm0 != varStructm1){
      stop("the residual variance model should be the same under both hypotheses")
    }
  }
  
  # get the structure of the random effects in m0 and m1
  if (randm0){
    vc0 <- eval(parse(text=paste("m0$modelStruc$reStruct$",nameRE,sep="")))
  }else{
    vc0 <- NULL
  }  
  vc1 <- eval(parse(text=paste("m1$modelStruc$reStruct$",nameRE,sep="")))
  
  # Structure of the covariance matrix
  if (randm0) covStruct0 <- class(vc0)[1]
  covStruct1 <- class(vc1)[1]
  
  # names of fixed and random effects
  if (randm0){
    namesRE0 <- colnames(m0$coefficients$random[[1]])
    namesFE0 <- names(m0$coefficients$fixed)
  }else{
    namesRE0 <- NULL
    if (inherits(m0,"lm")) namesFE0 <- names(stats::coefficients(m0))
    if (inherits(m0,"nls")) namesFE0 <- names(m0$m$getPars())
  }
  namesRE1 <- colnames(m1$coefficients$random[[1]])
  nameVarTested <- namesRE1[!(namesRE1%in%namesRE0)]
  
  namesFE1 <- names(m1$coefficients$fixed)
  nameFixedTested <- namesFE1[!(namesFE1%in%namesFE0)]
  
  # Throwing errors for cases not covered by the package
  #if ( !!!!!!! ) stop("Error: the current version of the package does not support more than 1 level of random effects")
  if (!prod(namesFE0 %in% namesFE1)) stop("the models should be nested, but it seems that some fixed effects are in m0 but not in m1")
  if (!prod(namesRE0 %in% namesRE1)) stop("the models should be nested, but it seems that some random effects are in m0 but not in m1")
  
  # dimension of the parameters
  nbFixEff0 <- length(namesFE0)
  nbFixEff1 <- length(namesFE1)
  nbRanEff0 <- length(namesRE0)
  nbRanEff1 <- length(namesRE1)
  nbRanEffTest <- nbRanEff1-nbRanEff0
  
  # retrieve names of parameters to identify their order in the FIM and in the cone
  #nbTotParam <- ncol(invfim)
  # get the covariance parameters
  Gamma1 <- extractVarCov(m1)
  covNames1 <- NULL
  covNames0 <- NULL
  if (length(Gamma1)>1 & !Matrix::isDiagonal(Gamma1)){
    lowDiag <- Gamma1
    lowDiag[lower.tri(lowDiag,diag = T)] <- NA
    posNonZeros <- which(lowDiag!=0,arr.ind = TRUE)
    rowNonZeros <- posNonZeros[,"row"]
    colNonZeros <- posNonZeros[,"col"]
    covNames1 <- sapply(1:nrow(posNonZeros),FUN=function(i){paste0("cor(",namesRE1[rowNonZeros[i]],",",namesRE1[colNonZeros[i]],")")})
  }
  if (randm0) {
    Gamma0 <- extractVarCov(m0)
    if (length(Gamma0)>1 & !Matrix::isDiagonal(Gamma0)){
      lowDiag <- Gamma0
      lowDiag[lower.tri(lowDiag,diag = T)] <- NA
      posNonZeros <- which(lowDiag!=0,arr.ind = TRUE)
      rowNonZeros <- posNonZeros[,"row"]
      colNonZeros <- posNonZeros[,"col"]
      covNames0 <- sapply(1:nrow(posNonZeros),FUN=function(i){paste0("cor(",namesRE0[rowNonZeros[i]],",",namesRE0[colNonZeros[i]],")")})
    }
  }
  
  if (!prod(covNames0%in%covNames1)) stop("the models should be nested but there are some covariances in m0 which are not in m1")
  
  nameParams0 <- c(namesFE0,paste0("sd(",namesRE0,")"),covNames0)
  nameParams1 <- c(namesFE1,paste0("sd(",namesRE1,")"),covNames1)
  paramTested <- !(nameParams1 %in% nameParams0)
  dd <- data.frame(names=nameParams1,tested=paramTested)
  dd$type <- c(rep("beta",nbFixEff1),substr(dd$names[(nbFixEff1+1):nrow(dd)],1,2))
  
  if (nrow(dd[dd$type=="co",])>0){
    dd$var1 <- dd$var2 <- dd$covTested <- dd$covInBlock <- NA
    ddco <- dd[dd$type=="co",]
    for (i in 1:nrow(ddco)){
      tmp<-gsub("cor","",ddco$names[i]);
      tmp<-gsub("[(]","",tmp);
      tmp<-gsub(")","",tmp);
      tmp<-strsplit(tmp,",");
      dd[dd$type=="co",][i,]$var1 <- tmp[[1]][1];
      dd[dd$type=="co",][i,]$var2 <- tmp[[1]][2];
      dd[dd$type=="co",][i,]$covTested <- dd[dd$type=="co",][i,]$tested*prod(!(unlist(tmp) %in% nameVarTested)) # identify covariances that are tested without testing the associate variances
      dd[dd$type=="co",][i,]$covInBlock <- dd[dd$type=="co",][i,]$tested*prod((unlist(tmp) %in% nameVarTested)) # identify covariances that are tested as part of a block of the covariance matrix tested
    }
    dd$covInBlock[dd$type=="sd"] <- 1
  }
  
  # Throwing errors for cases not covered by the package
  #if ( !!!!!!! ) stop("Error: the current version of the package does not support more than 1 level of random effects")
  #if (nbFixEff1 != nbFixEff0) stop("Error: the current version of the package does not support simultaneously testing means and variances. Models should have the same fixed effects")
  if (!prod(namesRE0 %in% namesRE1)) stop("the models should be nested, but it seems that some random effects are in m0 but not in m1")
  
  # get the dimension of the residual variance
  # it should be identical under H0 and H1
  if(is.null(m1$modelStruct$corStruct)){
    dimSigma=1
    if (randm0 & !is.null(m0$modelStruct$corStruct)) stop("The same model should be used under H0 and H1 for the residual covariance matrix")
  }else{
    # if no random effects under H0 then the model is fitted with lm which only allow for residual of dimension 1
    if (!randm0) stop("The same model should be used under H0 and H1 for the residual covariance matrix")
    dimSigma=1+length(m1$modelStruct$corStruct) # to be refined
  }
  nbCov1 <- length(nameParams1) - nbFixEff1 - nbRanEff1# - dimSigma
  nbCov0 <- length(nameParams0) - nbFixEff0 - nbRanEff0# - dimSigma
  
  diag <- (covStruct1 == "pdDiag" || covStruct1 == "pdIdent" || !("co" %in% dd$type))
  blockDiag <- (covStruct1 == "pdBlocked")
  full <- !diag & !blockDiag #!!! compound symmetry !!!
  struct <- c("diag","full","blockDiag")[c(diag,full,blockDiag)]
  
  if (blockDiag){
    dimBlock1 <- lengths(vc1) # gives the number of elements in each block
    dimBlock0 <- lengths(vc0)
    
    dd$block <- 0
    for(i in 1:length(dimBlock1)){
      nameREinBlock <- attr(vc1[[i]],"Dimnames")[[1]]
      loc <- sapply(1:length(nameREinBlock),FUN = function(x){grep(nameREinBlock[x],dd$names)})
      dd$block[loc] <- i
    }
    dd$block[dd$type=="beta"] <- 0
    dd$diagBlock <- NA
    dd$testInBlock <- 0
    for (i in 1:length(dimBlock1)){
      dd$diagBlock[dd$block==i] <- !("co" %in% dd$type[dd$block==i])
      dd$testInBlock[dd$block==i] <- max(dd$tested[dd$block==i])
    }
  }
  
  nameVarTested <- gsub("[(]","",nameVarTested)
  nameVarTested <- gsub("[)]","",nameVarTested)
  
  return(list(detailStruct=dd,
              nameVarTested=nameVarTested,
              nameFixedTested=nameFixedTested,
              dims=list(nbFE1=nbFixEff1,nbFE0=nbFixEff0,nbRE1=nbRanEff1,nbRE0=nbRanEff0,nbCov1=nbCov1,nbCov0=nbCov0,dimSigma=dimSigma),
              structGamma=struct))
}


#' @name extractVarCov.lme
#' @rdname extractVarCov.lme
#'
#' @title Extract covariance matrix 
#'
#' @description Extract covariance matrix of the random effects for a model fitted with nlme.
#'
#' @param m a fit from nlme package (either linear or nonlinear)
#' @export extractVarCov.lme
#' @export
extractVarCov.lme <- function(m){
  if (! (inherits(m,"nlme"))){
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


#' @name bootinvFIM.lme
#' @rdname bootinvFIM.lme
#'
#' @title Compute the inverse of the Fisher Information Matrix using parametric bootstrap
#'
#' @param m the model under H1
#' @param B the bootstrap sample size
#' @export bootinvFIM.lme
#' @export
bootinvFIM.lme <- function(m, B=1000){
  
  mySumm <- function(m,diagSigma=F) {
    beta <- nlme::fixef(m)
    resStd <- stats::sigma(m)
    Sigma <- extractVarCov(m)
    if(diagSigma){
      theta <- c(beta,diag(Sigma),resStd)
    }else{
      theta <- c(beta,Sigma[lower.tri(Sigma,diag = T)],resStd)
    }
    return(theta)
  }
  nonlin <- inherits(m,"nlme")
  
  if (!nonlin){
    bootstrap <- lmeresampler::bootstrap(m, mySumm, B = B, type = "parametric")
    bootstrap <- bootstrap$t[, colSums(bootstrap$t != 0) > 0]
    invfim <- cov(bootstrap)
    
    Gamma1 <- extractVarCov(m)
    namesRE <- colnames(m$coefficients$random[[1]])
    if (length(Gamma1)>1 & !Matrix::isDiagonal(Gamma1)){
      lowDiag <- Gamma1
      lowDiag[lower.tri(lowDiag,diag = F)] <- NA
      posNonZeros <- which(lowDiag!=0,arr.ind = TRUE)
      rowNonZeros <- posNonZeros[,"row"]
      colNonZeros <- posNonZeros[,"col"]
      covNames1 <- sapply(1:nrow(posNonZeros),FUN=function(i){
        if (rowNonZeros[i]==colNonZeros[i]){
          nme <- paste0("var(",namesRE[rowNonZeros[i]],")")
        }else{
          nme <- paste0("cov(",namesRE[rowNonZeros[i]],",",namesRE[colNonZeros[i]],")")
        }
        return(nme)
      })
    }else{
      if (length(Gamma1) == 1){
        covNames1 <- namesRE[1]
      }else{
        covNames1 <- paste0("var(",namesRE,")")
      }
    }
    namesParams <- c(names(m$coefficients$fixed),covNames1,"sd_residual")
    colnames(invfim) <- rownames(invfim) <- namesParams
  }else{
    beta <- nlme::fixef(m) # fixed effects
    resStd <- stats::sigma(m)
    Sigma <- extractVarCov(m)
    
    nind <- nrow(m$coefficients$random[[1]]) # only one level of random effects
    nrandEfft <- nrow(Sigma)
    namesRE <- colnames(m$coefficients$random[[1]])
    grpFactor <- names(m$groups)
    
    diagSigma <- Matrix::isDiagonal(Sigma)
    if (diagSigma){
      namesParams <- c(names(m$coefficients$fixed),paste0("var(",namesRE,")"),"sd_residual")
    }else{
      Sigma2 <- Sigma
      Sigma2[lower.tri(Sigma2)] <- NA
      nonZeros <- which(Sigma2!=0, arr.ind = T)
      
      namesCovParams <- character()
      for (l in 1:nrow(nonZeros)){
        if (nonZeros[l,1] == nonZeros[l,2]) namesCovParams <- c(namesCovParams,paste0("var(",colnames(Sigma)[nonZeros[l,1]],")"))
        if (nonZeros[l,1] != nonZeros[l,2]) namesCovParams <- c(namesCovParams,paste0("cov(",colnames(Sigma)[nonZeros[l,1]],",",colnames(Sigma)[nonZeros[l,2]],")"))
      }
      
      namesParams <- c(names(m$coefficients$fixed),namesCovParams,"sd_residual")
    }
    
    
    message(paste0("\t ...generating the B=",B," bootstrap samples ...\n"))
    thetaBoot <- numeric()
    
    b <- 1
    tbar <- utils::txtProgressBar(min=1,max=B,char = ".", style = 3)
    grpVar <- m$groups[,1]
    while (b <= B){  
      utils::setTxtProgressBar(tbar,b)
      phi <- t(chol(Sigma)%*%matrix(stats::rnorm(nrow(Sigma)*nind,0,1),ncol=nind))
      betaAll <- as.data.frame(matrix(rep(beta,nind),nrow=nind,byrow = T))
      betaAll <- cbind(betaAll,grp=unique(grpVar))
      names(betaAll) <- c(names(beta),grpFactor)
      pos <- names(betaAll)%in%namesRE # TRUE/FALSE to identify where are the random effects
      c <- 1
      
      ## !!! pb de comparaison de facteurs qd grpVar n'est pas forcément numérique ... il faudrait "enlever" le type
      
      nmeInd <- unique(grpVar)
      for (i in 1:length(pos)){
        if (pos[i]){
          betaAll[,i] <- betaAll[,i] + phi[,c]#c(sapply(nmeInd, FUN = function(j){betaAll[grpVar==j,i] + phi[j,c]}))
          c <- c+1
        }
      }
      
      # get name of response variable
      responseVar <- unlist(strsplit(as.character(nlme::getResponseFormula(m)),"~"))[[2]]
      # get name of internal nonlinear function
      nlmod <- unlist(strsplit(as.character(nlme::getCovariateFormula(m))[2],"[(]"))[1]
      # get names of all variables, and identify the covariables as the complement of response variable and grouping factor
      data.m <- nlme::getData(m)
      namesAllVar <- names(data.m)
      namesCov <- namesAllVar[! (namesAllVar %in% c(responseVar,grpFactor))]
      
      # build the part of the formula with the covariables, to be found in m@frame
      modAndCov <- paste(paste0("data.m$",namesCov),collapse=",")
      # build the part of the formula with the parameters
      posParamInBeta <- (1:ncol(betaAll))[names(betaAll) %in% names(beta)]
      modAndParam <- paste(paste0("betaAll[,",posParamInBeta,"]"),collapse=",")
      
      namesFE <- names(m$coefficients$fixed)
      simuResp <- lapply(1:nrow(betaAll),FUN = function(i){
        eval(parse(text=paste0(namesFE,"=",betaAll[i,posParamInBeta],sep=";")))
        with(data.m[eval(parse(text=paste0("data.m$",grpFactor,"==nmeInd[i]"))),],
             eval(parse(text=as.character(getCovariateFormula(m))[2])))
      })
      simuResp <- do.call("c",simuResp)
      
      data.m[,responseVar] <- simuResp  + stats::rnorm(nrow(betaAll),0,resStd)
      
      #############################
      fitInd <- suppressWarnings(try({setTimeLimit(10)
        stats::update(m, data=data.m)},silent=TRUE))
      
      if (!inherits(fitInd,"try-error")){
        thetaHatBoot <- mySumm(fitInd,diagSigma)
        names(thetaHatBoot) <- namesParams
        thetaBoot <- rbind(thetaBoot,thetaHatBoot)
        b <- b + 1
      }
    }
    
    invfim <- cov(thetaBoot) # empirical covariance matrix of the bootstrap samples as an estimator of the inverse of the FIM
    colnames(invfim) <- namesParams
    rownames(invfim) <- namesParams
  }
  
  return(invfim)
}

