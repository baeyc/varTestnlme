#' Approximation of the (inverse of the) Fisher Information Matrix
#'
#' Approximation of the inverse of the Fisher Information Matix via parametric bootstrap
#'
#' When the FIM is not available, this function provides an approximation of the FIM based on an estimate
#' of the covariancole matrix of the model's parameters obtained via parametric bootstrap.
#'
#' @name bootinvFIM
#' @aliases invFIM bootstrap
#'
#' @param m a fitted model that will be used as the basis of the parametric bootstrap (providing the initial maximum
#' likelihood estimate of the parameters and the modelling framework)
#' @param B the size of the bootstrap sample
#' @return the empirical covariancole matrix of the parameter estimates obtained on the bootstrap sample
#' @author Charlotte Baey <\email{charlotte.baey@univ-lille.fr}>
#'
#' @export bootinvFIM
#' @importFrom foreach %dopar%
bootinvFIM <- function(m, B=1000){

  mySummlme4 <- function(m,diagSigma=F) {
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
  mySummnlme <- function(m,diagSigma=F) {
    beta <- nlme::fixef(m)
    resStd <- stats::sigma(m)
    Sigma <- getVarCovnlme(m)
    if(diagSigma){
      theta <- c(beta,diag(Sigma),resStd)
    }else{
      theta <- c(beta,Sigma[lower.tri(Sigma,diag = T)],resStd)
    }
    return(theta)
  }

  pkg <- pckName(m)

  nonlin <- (class(m)[1] == "nlme" || class(m) %in% "nlmerMod")

  # Use bootMer functions if linear or generalized linear, otherwise code our own bootstrap
  if (!nonlin){
    if (pkg == "lme4"){
      bootstrap <- lme4::bootMer(m, mySummlme4, use.u = F, type = "parametric", nsim = B)
      bootstrap <- bootstrap$t[, colSums(bootstrap$t != 0) > 0]
      invfim <- cov(bootstrap)
      
      namesParams <- c(names(lme4::getME(m,"fixef")),paste0("cov_",names(lme4::getME(m,"theta"))),"residual")
      colnames(invfim) <- rownames(invfim) <- namesParams
    }
    if (pkg == "nlme"){
      bootstrap <- lmeresampler::parametric_bootstrap(m, mySummnlme, B = B)
      bootstrap <- bootstrap$t[, colSums(bootstrap$t != 0) > 0]
      invfim <- cov(bootstrap)
      
      Gamma1 <- getVarCovnlme(m)
      namesRE <- colnames(m$coefficients$random[[1]])
      if (length(Gamma1)>1 & !Matrix::isDiagonal(Gamma1)){
        print("truc")
        lowDiag <- Gamma1
        lowDiag[lower.tri(lowDiag,diag = F)] <- NA
        posNonZeros <- which(lowDiag!=0,arr.ind = TRUE)
        rowNonZeros <- posNonZeros[,"row"]
        colNonZeros <- posNonZeros[,"col"]
        covNames1 <- sapply(1:nrow(posNonZeros),FUN=function(i){
          if (rowNonZeros[i]==colNonZeros[i]){
            nme <- paste0("sd(",namesRE[rowNonZeros[i]],")")
          }else{
            nme <- paste0("cor(",namesRE[rowNonZeros[i]],",",namesRE[colNonZeros[i]],")")
          }
        return(nme)
        })
      }else{
        if (length(Gamma1) == 1){
          covNames1 <- namesRE[1]
        }else{
          covNames1 <- paste0("sd(",namesRE,")")
        }
      }
      namesParams <- c(names(m$coefficients$fixed),covNames1,"residual")
      colnames(invfim) <- rownames(invfim) <- namesParams
    }
  }else{
    if (pkg == "lme4"){
      beta <- lme4::fixef(m) # fixed effects
      resStd <- stats::sigma(m)
      grpFactor <- unique(names(lme4::getME(m,"cnms")))
      vc <- lme4::VarCorr(m)
      Sigma <- getVarCovlme4(m)
      diagSigma <- Matrix::isDiagonal(Sigma)
      
      # Generate B bootstrap samples
      nind <- lme4::getME(m,"l_i")
      nrandEfft <- nrow(Sigma)
      namesRE <- lme4::getME(m,"cnms")
      
      namesParams <- c(names(lme4::fixef(m)),paste0("cov_",names(lme4::getME(m,"theta"))),"residual")

      message(paste0("\t ...generating the B=",B," bootstrap samples ...\n"))
      no_cores <- max(1,parallel::detectCores() - 1)
      # Initiate cluster
      doParallel::registerDoParallel(no_cores)
      
      thetaBoot <- numeric()
      
      # thetaBoot <- foreach::foreach(i=1:B, 
      #                                  .packages='varTestnlme', 
      #                                  .combine=rbind) %dopar% {      
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

        fitInd <- suppressWarnings(try({setTimeLimit(2)
                                        stats::update(m, data=d)},silent=TRUE))
        #fitInd <- update(m, data=d)
        
        if (!inherits(fitInd,"try-error")){
          thetaHatBoot <- mySummlme4(fitInd,diagSigma)
          names(thetaHatBoot) <- namesParams
          thetaBoot <- rbind(thetaBoot,thetaHatBoot)
          b <- b + 1
        }
        #else{
        #  print("fail")
        #}
      }
    }

    if (pkg == "nlme"){
      beta <- nlme::fixef(m) # fixed effects
      resStd <- stats::sigma(m)
      Sigma <- getVarCovnlme(m)

      nind <- nrow(m$coefficients$random[[1]]) # only one level of random effects
      nrandEfft <- nrow(Sigma)
      namesRE <- colnames(m$coefficients$random[[1]])
      grpFactor <- names(m$groups)

      namesParams <- c(names(m$coefficients$fixed),paste0("sd(",namesRE,")"),"residual")
      diagSigma <- Matrix::isDiagonal(Sigma)
      
      message(paste0("\t ...generating the B=",B," bootstrap samples ...\n"))
      no_cores <- max(1,parallel::detectCores() - 1)
      # Initiate cluster
      doParallel::registerDoParallel(no_cores)
    
      thetaBoot <- numeric()
        
      b <- 1
      tbar <- utils::txtProgressBar(min=1,max=B,char = ".", style = 3)
      grpVar <- m$groups[,1]
      while (b <= B){  
        utils::setTxtProgressBar(tbar,b)
        phi <- t(chol(Sigma)%*%matrix(stats::rnorm(nrow(Sigma)*nind,0,1),ncol=nind))
        betaAll <- as.data.frame(matrix(rep(beta,nind),nrow=nind,byrow = T))
        betaAll <- cbind(betaAll,grp=grpVar)
        names(betaAll) <- c(names(beta),grpFactor)
        pos <- names(betaAll)%in%namesRE # TRUE/FALSE to identify where are the random effects
        c <- 1
        
        ## !!! pb de comparaison de facteurs qd grpVar n'est pas forcément numérique ... il faudrait "enlever" le type
        
        nmeInd <- unique(grpVar)
        for (i in 1:length(pos)){
          if (pos[i]){
            betaAll[,i] <- c(sapply(nmeInd, FUN = function(j){betaAll[grpVar==j,i] + phi[j,c]}))
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
        eval(parse(text=paste0(namesFE,"=",beta,sep=";")))
        simuResp <- with(data.m,
                         eval(parse(text=as.character(getCovariateFormula(m))[2])))

        data.m[,responseVar] <- simuResp  + stats::rnorm(nrow(betaAll),0,resStd)

        #############################
        fitInd <- suppressWarnings(try({setTimeLimit(2)
                                        stats::update(m, data=data.m)},silent=TRUE))
                        
        if (!inherits(fitInd,"try-error")){
          thetaHatBoot <- mySummnlme(fitInd,diagSigma)
          names(thetaHatBoot) <- namesParams
          thetaBoot <- rbind(thetaBoot,thetaHatBoot)
          b <- b + 1
        }
      }
      
      message("\n\n")
    }

    if (pkg == "saemix"){
      # newdata <- saemix::simul.saemix(m, nsim=B, predictions = T, res.var = T)
      # 
      # no_cores <- max(1,parallel::detectCores() - 1)
      # # Initiate cluster
      # doParallel::registerDoParallel(no_cores)
      # 
      # thetaBoot <- foreach::foreach(i=1:B, 
      #                               .packages='varTestnlme', 
      #                               .combine=rbind) %dopar% {   
      #   dataB <- cbind(m@data@data,newdata@sim.data@datasim[newdata@sim.data@datasim$irep==i,])
      #   saemixDataB <- saemixData(dataB, header = T, name.group = m@data@name.group, name.predictors = m@data@name.predictors, name.response = "ysim")
      #   newoptions <- m1@options
      #   newoptions$displayProgress <- F
      #   newoptions$fim <- F
      #   newoptions$ll.is <- F
      #   newoptions$print <- F
      #   newoptions$nbdisplay <- 1e+06
      #   saemix::saemix(model = m@model, data = saemixDataB, control = newoptions)
      #   mb <- update(m,data=newdata)
      #   tmp <- as.data.frame(VarCorr(mb))
      #   c(mb@beta,tmp$vcov)
      #   # TO FINISH
      #}
    }

    invfim <- cov(thetaBoot) # empirical covariance matrix of the bootstrap samples as an estimator of the inverse of the FIM
    colnames(invfim) <- namesParams
    rownames(invfim) <- namesParams
  }

  return(invfim)
}

