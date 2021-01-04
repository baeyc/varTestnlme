#' @name extractStruct.merMod
#' @rdname extractStruct.merMod
#'
#' @title Extract model structure
#'
#' @param m1 the fit under H1
#' @param m0 the fit under H0
#' @param randm0 a boolean indicating whether random effects are present in m0
extractStruct.merMod <- function(m1,m0,randm0){
  # get package used to fit m0
  pkgm0 <- class(m0)[1]
  
  # name of the grouping factor
  nameRE <- names(m1@flist)
  if (length(nameRE)>1) stop("Error: the package does not currently support more than one level of random effects")
  
  # dimension of the parameters
  nbFixEff1 <- lme4::getME(m1,"p")
  # names of fixed and random effects
  if (randm0){
    namesRE0 <- unlist(m0@cnms)
    namesFE0 <- names(lme4::getME(m0,"fixef"))
    nbFixEff0 <- lme4::getME(m0,"p")
  }else{
    namesRE0 <- NULL
    if (pkgm0 %in% c("lm","glm")) namesFE0 <- names(stats::coefficients(m0))
    if (pkgm0 == "nls") namesFE0 <- names(m0$m$getPars())
    nbFixEff0 <- length(namesFE0)
  }
  namesFE1 <- names(lme4::getME(m1,"fixef"))
  namesRE1 <- unlist(m1@cnms)
  nameVarTested <- namesRE1[!(namesRE1%in%namesRE0)]
  nameFixedTested <- namesFE1[!(namesFE1%in%namesFE0)]
  nbRanEff0 <- length(namesRE0)
  nbRanEff1 <- length(namesRE1)
  nbRanEffTest <- nbRanEff1-nbRanEff0
  
  #if (nbRanEff0==nbRanEff1) stop("Error: there are the same number of random effects in models m0 and m1. Please check the models' formulation.")
  
  # Structure of the covariance matrix (diagonal, blocked or full)
  nbCompVar1 <- lme4::getME(m1,"devcomp")$dims["nth"]
  if (randm0){
    nbCompVar0 <- lme4::getME(m0,"devcomp")$dims["nth"]
  }else{
    nbCompVar0 <- 0
  }
  diag <- (nbCompVar1==nbRanEff1)
  full <- (nbCompVar1==(nbRanEff1*(nbRanEff1+1)/2))
  if (nbCompVar1==1) full <- FALSE
  blockDiag <- !diag & !full
  struct <- c("diag","full","blockDiag")[c(diag,full,blockDiag)]
  nbCov1 <- nbCompVar1 - nbRanEff1
  nbCov0 <- nbCompVar0 - nbRanEff0
  
  if (nbCov1 < nbCov0) stop("Error: the models should be nested but there are some covariances in m0 which are not in m1")
  
  # CHECK IF ML WAS USED AND NOT REML
  if (lme4::isREML(m1)) stop("Error: the models should be fitted using Maximum Likelihood (ML) instead of Restricted ML (REML)")
  if (randm0) if (lme4::isREML(m0)) stop("Error: the models should be fitted using Maximum Likelihood (ML) instead of Restricted ML (REML)")
  
  # retrieve names of parameters to identify their order in the FIM and in the cone
  #nbTotParam <- ncol(invfim1)
  if (randm0) nameParams0 <- c(namesFE0,paste0("cov_",names(lme4::getME(m0,"theta"))),"residual")
  if (!randm0) nameParams0 <- c(namesFE0,"residual")
  nameParams1 <- c(namesFE1,paste0("cov_",names(lme4::getME(m1,"theta"))),"residual")
  paramTested <- !(nameParams1 %in% nameParams0)
  
  # create a dataset with the list of parameters and: their tyoe (fixed, variance or correlation)
  # whether they are tested equal to 0 or not, and if they are tested, if it's as a subpart of a block
  # or in a block which is fully tested
  dd <- data.frame(names=nameParams1,tested=paramTested)
  dd$type <- c(rep("beta",nbFixEff1),substr(dd$names[(nbFixEff1+1):nrow(dd)],1,2))
  indRes <- as.numeric(rownames(dd[dd$type=="re",])) # get indices of residual parameters
  dd <- dd[-indRes,]
  ddco <- dd[dd$type=="co",]
  
  if (nrow(ddco)>0){
    dd$var1 <- dd$var2 <- dd$covTested <- dd$covInBlock <- NA
    isVariance <- numeric()
    for (i in 1:nrow(ddco)){
      tmp<-gsub(paste("cov_",nameRE,".",sep=""),"",ddco$names[i]);
      tmp<-strsplit(tmp,"[.]");
      dd[dd$type=="co",][i,]$var1 <- tmp[[1]][1];
      dd[dd$type=="co",][i,]$var2 <- tmp[[1]][2];
      if (length(tmp[[1]])==1) isVariance <- c(isVariance,i);
      dd[dd$type=="co",][i,]$covTested <- dd[dd$type=="co",][i,]$tested*prod(!(unlist(tmp) %in% nameVarTested)) # identify covariances that are tested without testing the associate variances
      dd[dd$type=="co",][i,]$covInBlock <- dd[dd$type=="co",][i,]$tested*prod((unlist(tmp) %in% nameVarTested)) # identify covariances that are tested as part of a block of the covariance matrix tested
    }
    dd[nrow(dd)-nrow(ddco)+isVariance,]$type <- "sd"
    dd$covInBlock[dd$type=="sd"] <- 1
  }else{
    dd$covInBlock[dd$type=="beta"] <- NA
    dd$covInBlock[dd$type!="beta"] <- 1
  }
  
  if(blockDiag){
    dimBlock1 <- lengths(m1@cnms) # nb of ramdom effects per block
    dimBlock0 <- lengths(m0@cnms)
    
    dd$block <- 0
    for(i in 1:length(dimBlock1)){
      nameREinBlock <- m1@cnms[[i]]
      loc <- sapply(1:length(nameREinBlock),FUN = function(x){grep(nameREinBlock[x],dd$names)})
      dd$block[unlist(loc)] <- i
    }
    dd$block[dd$type=="beta"] <- 0
  }
  
  return(list(detailStruct=dd,
              nameVarTested=nameVarTested,
              nameFixedTested=nameFixedTested,
              dims=list(nbFE1=nbFixEff1,nbFE0=nbFixEff0,nbRE1=nbRanEff1,nbRE0=nbRanEff0,nbCov1=nbCov1,nbCov0=nbCov0,dimSigma=1*!(pkgm0%in%c("glm"))),
              structGamma=struct))
}

#' @name extractVarCov.merMod
#' @rdname extractVarCov.merMod
#'
#' @title Extract covariance matrix 
#'
#' @description Extract covariance matrix of the random effects for a model fitted with lme4.
#'
#' @param m a fit from lme4 package (either linear or nonlinear)
extractVarCov.merMod <- function(m){
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


#' @name bootinvFIM.merMod
#' @rdname bootinvFIM.merMod
#'
#' @title Compute the inverse of the Fisher Information Matrix using parametric bootstrap
#'
#' @param m the model under H1
#' @param B the bootstrap sample size
bootinvFIM.merMod <- function(m, B=1000){
  
  mySumm <- function(m,diagSigma=F) {
    beta <- lme4::fixef(m)
    resStd <- stats::sigma(m)
    grpFactor <- names(lme4::getME(m,"cnms"))
    vc <- lme4::VarCorr(m)
    Sigma <- as.matrix(Matrix::bdiag(as.matrix(vc)))
    if(diagSigma){
      theta <- c(beta,diag(Sigma),resStd)
    }else{
      theta <- c(beta,Sigma[lower.tri(Sigma,diag = T)],resStd)
    }
    return(theta)
  }
  
  nonlin <- (class(m) %in% "nlmerMod")
  
  # Use bootMer functions if linear or generalized linear, otherwise code our own bootstrap
  if (!nonlin){
    bootstrap <- lme4::bootMer(m, mySumm, use.u = F, type = "parametric", nsim = B)
    bootstrap <- bootstrap$t[, colSums(bootstrap$t != 0) > 0]
    invfim <- cov(bootstrap)
    
    namesParams <- c(names(lme4::getME(m,"fixef")),paste0("cov_",names(lme4::getME(m,"theta"))),"residual")
    colnames(invfim) <- rownames(invfim) <- namesParams
  }else{
    beta <- lme4::fixef(m) # fixed effects
    resStd <- stats::sigma(m)
    grpFactor <- unique(names(lme4::getME(m,"cnms")))
    vc <- lme4::VarCorr(m)
    Sigma <- extractVarCov(m)
    diagSigma <- Matrix::isDiagonal(Sigma)
    
    # Generate B bootstrap samples
    nind <- lme4::getME(m,"l_i")
    nrandEfft <- nrow(Sigma)
    namesRE <- lme4::getME(m,"cnms")
    
    namesParams <- c(names(lme4::fixef(m)),paste0("cov_",names(lme4::getME(m,"theta"))),"residual")
    
    message(paste0("\t ...generating the B=",B," bootstrap samples ...\n"))
    #no_cores <- max(1,parallel::detectCores() - 1)
    # Initiate cluster
    #doParallel::registerDoParallel(no_cores)
    
    thetaBoot <- numeric()
    
    b <- 1
    tbar <- utils::txtProgressBar(min=1,max=B,char = ".", style = 3)
    grpVar <- m@frame[,grpFactor]
    nmeInd <- unique(grpVar)
    while (b <= B){  
      utils::setTxtProgressBar(tbar,b)      
      phi <- t(t(chol(Sigma))%*%matrix(stats::rnorm(nrandEfft*nind,0,1),ncol=nind))
      betaAll <- as.data.frame(matrix(rep(beta,nind),nrow=nind,byrow = TRUE))
      betaAll <- cbind(betaAll,grp=grpVar)
      names(betaAll) <- c(names(beta),grpFactor)
      pos <- names(betaAll)%in%unlist(namesRE) # TRUE/FALSE to identify where are the random effects
      c <- 1
      for (i in 1:length(pos)){
        if (pos[i]){
          betaAll[,i] <- c(sapply(nmeInd, FUN = function(j){betaAll[grpVar==j,i] + phi[j,c]}))
          c <- c+1
        }
      }
      
      # get name of response variable
      responseVar <- unlist(strsplit(as.character(m@call)[2],"[~]"))[1]
      responseVar <- gsub(" ","",responseVar)
      # get name of internal nonlinear function
      nlmod <- unlist(strsplit(as.character(m@resp$nlmod)[2],"[(]"))[1]
      # get names of all variables, and identify the covariables as the complement of response variable and grouping factor
      namesAllVar <- names(m@frame)
      namesCov <- namesAllVar[! (namesAllVar %in% c(responseVar,grpFactor))]
      
      # build the part of the formula with the covariables, to be found in m@frame
      modAndCov <- paste(paste0("m@frame$",namesCov),collapse=",")
      # build the part of the formula with the parameters
      posParamInBeta <- (1:ncol(betaAll))[names(betaAll) %in% names(beta)]
      modAndParam <- paste(paste0("betaAll[,",posParamInBeta,"]"),collapse=",")
      
      simuResp <- eval(parse(text=paste0(nlmod,"(",modAndCov,",",modAndParam,")"))) 
      
      d <- m@frame
      d[,responseVar] <- simuResp + stats::rnorm(nrow(betaAll),0,resStd)
      
      fitInd <- suppressWarnings(try({setTimeLimit(100)
        stats::update(m, data=d, start=beta)},silent=TRUE))
      #fitInd <- update(m, data=d)
      
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

