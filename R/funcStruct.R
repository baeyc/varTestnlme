#' @name funcStruct
#' @rdname funcStruct
#'
#' @title Extracting models structures
#'
#' @description Functions extracting the structure of the models for each package \code{nlme}, \code{lme4}
#' and \code{saemix}.
#'
#'
#' @param m1 the model under H1
#' @param m0 the model under H0
#' @param randm0 a boolean stating whether the model under H0 contains any random effect
#' @param linmodel (only for \code{modelStructlme4}) a boolean to specify whether the model is linear or not
#'
#' @return A list with the following components:
#' \item{\code{detailStruct}}{a data frame containing 8 variables: \code{name} with the name of the model parameters, \code{var1} and
#' \code{var2} with the names of the two variances associated with each covariance parameter, \code{type} giving the type of parameter
#' (\code{beta} for fixed effects, \code{sd} for variances and \code{co} for covariances), \code{tested} equal to \code{TRUE} if the parameter is
#' tested and \code{FALSE} otherwise, \code{block} giving the block number to which the variance component parameter belongs (equal 0 for
#' fixed effects), \code{covTested} indicating whether a covariance is tested without the associated variances being tested,
#' and \code{covInBlock} indicating whether a covariance is tested within a block of the complete covariance matrix}
#' \item{\code{dims}}{a list with the dimensions of the models (\code{nbFE1} and \code{nbFE0} the number of fixed effects in m1 and m0,
#' \code{nbRE1} and \code{nbRE0} the number of random effects in m1 and m0 and \code{dimSigma} the number of residual error parameters)}
#' \item{\code{structGamma}}{the structure of the covariance matrix of the random effects as a list of three logical elements: \code{diag},
#' \code{full} and \code{blockDiag}, equal to \code{TRUE} if the matrix is diagonal, full or block-diagonal respectively.}
#' 
NULL

#' @rdname funcStruct
modelStructnlme <- function(m1,m0,randm0){
  # get package used to fit m0
  pkgm0 <- class(m0)[1]
  
  # name of the random grouping variable
  nameRE <- names(m1$groups)
  
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
    if (pkgm0 == "lm") namesFE0 <- names(stats::coefficients(m0))
    if (pkgm0 == "nls") namesFE0 <- names(m0$m$getPars())
  }
  namesRE1 <- colnames(m1$coefficients$random[[1]])
  nameVarTested <- namesRE1[!(namesRE1%in%namesRE0)]

  namesFE1 <- names(m1$coefficients$fixed)
  nameFixedTested <- namesFE1[!(namesFE1%in%namesFE0)]
  
  # Throwing errors for cases not covered by the package
  #if ( !!!!!!! ) stop("Error: the current version of the package does not support more than 1 level of random effects")
  if (!prod(namesFE0 %in% namesFE1)) stop("Error: the models should be nested, but it seems that some fixed effects are in m0 but not in m1")
  if (!prod(namesRE0 %in% namesRE1)) stop("Error: the models should be nested, but it seems that some random effects are in m0 but not in m1")
  
  # dimension of the parameters
  nbFixEff0 <- length(namesFE0)
  nbFixEff1 <- length(namesFE1)
  nbRanEff0 <- length(namesRE0)
  nbRanEff1 <- length(namesRE1)
  nbRanEffTest <- nbRanEff1-nbRanEff0

  # retrieve names of parameters to identify their order in the FIM and in the cone
  #nbTotParam <- ncol(invfim)
  # get the covariance parameters
  Gamma1 <- getVarCovnlme(m1)
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
    Gamma0 <- getVarCovnlme(m0)
    if (length(Gamma0)>1 & !Matrix::isDiagonal(Gamma0)){
      lowDiag <- Gamma0
      lowDiag[lower.tri(lowDiag,diag = T)] <- NA
      posNonZeros <- which(lowDiag!=0,arr.ind = TRUE)
      rowNonZeros <- posNonZeros[,"row"]
      colNonZeros <- posNonZeros[,"col"]
      covNames0 <- sapply(1:nrow(posNonZeros),FUN=function(i){paste0("cor(",namesRE0[rowNonZeros[i]],",",namesRE0[colNonZeros[i]],")")})
    }
  }
  
  if (!prod(covNames0%in%covNames1)) stop("Error: the models should be nested but there are some covariances in m0 which are not in m1")
  
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
  if (!prod(namesRE0 %in% namesRE1)) stop("Error: the models should be nested, but it seems that some random effects are in m0 but not in m1")

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

  return(list(detailStruct=dd,
              nameVarTested=nameVarTested,
              nameFixedTested=nameFixedTested,
              dims=list(nbFE1=nbFixEff1,nbFE0=nbFixEff0,nbRE1=nbRanEff1,nbRE0=nbRanEff0,nbCov1=nbCov1,nbCov0=nbCov0,dimSigma=dimSigma),
              structGamma=list(diag=diag,full=full,blockDiag=blockDiag)))
}


#' @rdname funcStruct
modelStructlme4 <- function(m1,m0,linmodel,randm0){
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
  blockDiag <- !diag & !full
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
              structGamma=list(diag=diag,full=full,blockDiag=blockDiag)))
}


#' @rdname funcStruct
modelStructsaemix <- function(m1,m0,randm0){
  # get package used to fit m0
  pkgm0 <- class(m0)[1]
  
  # dimension of the parameters
  nbFixEff0 <- sum(m0@model@fixed.estim>0)
  nbFixEff1 <- sum(m1@model@fixed.estim>0)
  # names of fixed and random effects
  if (randm0){
    namesRE0 <- m0@model@name.random
    namesFE0 <- m0@model@name.fixed[m0@model@fixed.estim>0]
  }else{
    namesRE0 <- NULL
    if (pkgm0 == "lm" || pkgm0 == "glm") namesFE0 <- names(stats::coefficients(m0))
    if (pkgm0 == "nls") namesFE0 <- names(m0$m$getPars())
  }
  namesRE1 <- m1@model@name.random
  namesFE1 <- m1@model@name.fixed[m1@model@fixed.estim>0]
  nameFixedTested <- namesFE1[!(namesFE1%in%namesFE0)]
  nameVarTested <- namesRE1[!(namesRE1%in%namesRE0)]
  #nameVarTested <- gsub("omega2.","",nameVarTested)
  nbRanEff0 <- length(namesRE0)
  nbRanEff1 <- length(namesRE1)
  nbRanEffTest <- nbRanEff1-nbRanEff0

  covStruct1 <- m1@model@covariance.model
  if (randm0) covStruct0 <- m0@model@covariance.model

  diag <- (sum(covStruct1)==sum(diag(covStruct1)))
  full <- (min(covStruct1)==1)
  blockDiag <- !diag & !full
  
  nbCompVar1 <- sum(m1@model@covariance.model)
  nbCompVar0 <- sum(m0@model@covariance.model)
  nbCov1 <- (nbCompVar1 - nbRanEff1)/2
  nbCov0 <- (nbCompVar0 - nbRanEff0)/2

  # dimension of residual error
  if (m1@model@error.model == "combined"){
    dimSigma <- 2
  }else{
    dimSigma <- 1
  }

  if(nbRanEff0==nbRanEff1) stop("Error: there are the same number of random effects in models m0 and m1. Please check the models' formulation.")
  if(randm0) if (min(covStruct1-covStruct0) < 0) stop("Error: the models should be nested, but it seems that some random effects are in m0 but not in m1")
  
  # create a dataset with the list of parameters and: their tyoe (fixed, variance or correlation)
  # whether they are tested equal to 0 or not, and if they are tested, if it's as a subpart of a block
  # or in a block which is fully tested  
  dd <- expand.grid(rownames(covStruct1),rownames(covStruct1))
  names(dd) <- c("var1","var2")
  dd$names <- paste(dd[,1],dd[,2],sep=".")
  dd$lowertri <- as.vector(lower.tri(covStruct1,diag = T))
  dd$zero1 <- as.vector(covStruct1==0)
  dd$diag <- as.vector(lower.tri(covStruct1,diag = T)) - as.vector(lower.tri(covStruct1,diag = F))
  dd$type <- ifelse(dd$diag==1,"sd","co")
  if (randm0){
    dd$tested <- as.vector(covStruct0==0)
  }else{
    dd$tested <- TRUE
  }
  dd <- dd[dd$lowertri & !dd$zero1,]
  dd$zero1 <- dd$lowertri <- dd$diag <- NULL

  dd <- rbind(data.frame(names=namesFE1,var1=rep(NA,length(namesFE1)),var2=rep(NA,length(namesFE1)),type=rep("beta",length(namesFE1)),tested=rep(FALSE,length(namesFE1))),dd)
  ddco <- dd[dd$type=="co",]
  
  if (nrow(ddco)>0){
    dd$covTested <- dd$covInBlock <- NA
    for (i in 1:nrow(ddco)){
      dd[dd$type=="co",][i,]$covTested <- dd[dd$type=="co",][i,]$tested*(!(dd$var1[dd$type=="co"][i] %in% nameVarTested) & !(dd$var2[dd$type=="co"][i] %in% nameVarTested))
      dd[dd$type=="co",][i,]$covInBlock <- dd[dd$type=="co",][i,]$tested*(dd$var1[dd$type=="co"][i] %in% nameVarTested) & (dd$var2[dd$type=="co"][i] %in% nameVarTested) # identify covariances that are tested as part of a block of the covariance matrix tested
    }
    dd$covInBlock[dd$type=="sd"] <- 1
  }else{
    dd$covInBlock[dd$type=="beta"] <- NA
    dd$covInBlock[dd$type!="beta"] <- 1
  }
  
  if (blockDiag){
    # identify block structure using spectral clustering
    D <- apply(covStruct1,1,sum)
    L <- diag(nrow(covStruct1)) - diag(1/sqrt(D))%*% covStruct1 %*% diag(1/sqrt(D))
    ev <- eigen(L,symmetric=T)$values
    k <- length(ev[ev<1e-10])
    blocks <- anocva::spectralClustering(covStruct1,k)
    dd$block <- 0
    dd$block[dd$type=="sd"] <- blocks

    ddco <- dd[dd$type=="co",]
    for (i in 1:k){
      nameREinBlock <- dd$var1[dd$block==i]
      loc <- sapply(1:length(nameREinBlock),FUN = function(x){grep(nameREinBlock[x],dd$names)})
      dd$block[unlist(loc)] <- i
    }
    dd$block[dd$type=="beta"] <- 0
  }
  rownames(dd) <- 1:nrow(dd) # reinitialize rownames (used to be indices of elements in the covariance matrix)
  if (full){
    dd$block <- 0
    dd$block[dd$tested] <- 1
  }
  if (diag){
    dd$block <- 1:nrow(dd)
    dd$block[!dd$tested] <- 0
  }

  return(list(detailStruct=dd,
              nameVarTested=nameVarTested,
              nameFixedTested=nameFixedTested,
              dims=list(nbFE1=nbFixEff1,nbFE0=nbFixEff0,nbRE1=nbRanEff1,nbRE0=nbRanEff0,nbCov1=nbCov1,nbCov0=nbCov0,dimSigma=dimSigma),
              structGamma=list(diag=diag,full=full,blockDiag=blockDiag)))
}


pckName <- function(m1){
  cl1 <- class(m1)
  if ("lme" %in% cl1) {
    pkg <- "nlme"
  }else if(cl1 %in% c("lmerMod","glmerMod","nlmerMod")){
    pkg <- "lme4"
  }else if(cl1 %in% "SaemixObject"){
    pkg <- "saemix"
  }else{
    stop("Error: the current version of the package only supports models fitted by nlme, lme4 or saemix packages")
  }
  return(pkg)
}

