##' @name funcStruct
##' @rdname funcStruct
##'
##' @title Extracting models structures
##'
##' @description Functions extracting the structure of the models for each package \code{\link[nlme]{nlme}}, \code{\link[lme4]{lme4}}
##' and \code{\link[saemix]{saemix}}
##'
##'
##' @param m1 the model under H1
##' @param m0 the model under H0
##' @param linmodel (only for \code{modelStructlme4}) a boolean to specify whether the model is linear or not
##'
##' @return A list with the following components:
##' \item{\code{detailStruct}}{a data frame containing 8 variables: \code{name} with the name of the model parameters, \code{var1} and
##' \code{var2} with the names of the two variances associated with each covariance parameter, \code{type} giving the type of parameter
##' (\code{beta} for fixed effects, \code{sd} for variances and \code{co} for covariances), \code{tested} equal to \code{TRUE} if the parameter is
##' tested and \code{FALSE} otherwise, \code{block} giving the block number to which the variance component parameter belongs (equal 0 for
##' fixed effects), \code{covTested} indicating whether a covariance is tested without the associated variances being tested,
##' and \code{covInBlock} indicating whether a covariance is tested within a block of the complete covariance matrix}
##' \item{\code{dims}}{a list with the dimensions of the models (\code{nbFE1} and \code{nbFE0} the number of fixed effects in m1 and m0,
##' \code{nbRE1} and \code{nbRE0} the number of random effects in m1 and m0 and \code{dimSigma} the number of residual error parameters)}
##' \item{\code{structGamma}}{the structure of the covariance matrix of the random effects as a list of three logical elements: \code{diag},
##' \code{full} and \code{blockDiag}, equal to \code{TRUE} if the matrix is diagonal, full or block-diagonal respectively.}
NULL

##' @rdname funcStruct
modelStructnlme <- function(m1,m0){
  # name of the random grouping variable
  nameRE <- names(m1$groups)

  vc0 <- eval(parse(text=paste("m0$modelStruc$reStruct$",nameRE,sep="")))
  vc1 <- eval(parse(text=paste("m1$modelStruc$reStruct$",nameRE,sep="")))

  # Structure of the covariance matrix
  covStruct0 <- class(vc0)[1]
  covStruct1 <- class(vc1)[1]

  namesRE0 <- attr(vc0,"Dimnames")[[1]]
  namesRE1 <- attr(vc1,"Dimnames")[[1]]
  nameVarTested <- namesRE1[!(namesRE1%in%namesRE0)]

  # dimension of the parameters
  nbFixEff0 <- m0$dims$ncol[2]
  nbFixEff1 <- m1$dims$ncol[2]
  nbRanEff0 <- length(namesRE0)
  nbRanEff1 <- length(namesRE1)
  nbRanEffTest <- nbRanEff1-nbRanEff0

  # retrieve names of parameters to identify their order in the FIM and in the cone
  #nbTotParam <- ncol(invfim)
  int1 <- nlme::intervals(m1)
  int0 <- nlme::intervals(m0)
  nameParams1 <- c(rownames(int1$fixed),eval(parse(text=paste("rownames(int1$reStruc$",nameRE,")",sep=""))))
  nameParams0 <- c(rownames(int0$fixed),eval(parse(text=paste("rownames(int0$reStruc$",nameRE,")",sep=""))))
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
  if (nbFixEff1 != nbFixEff0) stop("Error: the current version of the package does not support simultaneously testing means and variances. Models should have the same fixed effects")
  if (!prod(namesRE0 %in% namesRE1)) stop("Error: the models should be nested, but it seems that some random effects are in m0 but not in m1")

  if(is.null(m1$modelStruct$corStruct)){ # get the dimension of the residual variance
    dimSigma=1
  }else{
    dimSigma=1+length(m1$modelStruct$corStruct) # to be refined
  }

  diag <- (covStruct1 == "pdDiag" || covStruct1 == "pdIdent" || !("co" %in% dd$type))
  blockDiag <- (covStruct1 == "pdBlocked")
  full <- !diag & !blockDiag

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
              dims=list(nbFE1=nbFixEff1,nbFE0=nbFixEff0,nbRE1=nbRanEff1,nbRE0=nbRanEff0,dimSigma=dimSigma),
              structGamma=list(diag=diag,full=full,blockDiag=blockDiag)))
}


##' @rdname funcStruct
modelStructlme4 <- function(m1,m0,linmodel){
  # name of the grouping factor
  nameRE <- names(m1@flist)
  if (length(nameRE)>1) stop("Error: the package does not currently support more than one level of random effects")

  # dimension of the parameters
  nbFixEff0 <- lme4::getME(m0,"p")
  nbFixEff1 <- lme4::getME(m1,"p")
  namesRE0 <- unlist(m0@cnms)
  namesRE1 <- unlist(m1@cnms)
  nameVarTested <- namesRE1[!(namesRE1%in%namesRE0)]
  nbRanEff0 <- length(namesRE0)
  nbRanEff1 <- length(namesRE1)
  nbRanEffTest <- nbRanEff1-nbRanEff0

  if (nbRanEff0==nbRanEff1) stop("Error: there are the same number of random effects in models m0 and m1. Please check the models' formulation.")

  # Structure of the covariance matrix (diagonal, blocked or full)
  nbCompVar1 <- lme4::getME(m1,"devcomp")$dims["nth"]
  nbCompVar0 <- lme4::getME(m0,"devcomp")$dims["nth"]
  diag <- (nbCompVar1==nbRanEff1)
  full <- (nbCompVar1==(nbRanEff1*(nbRanEff1+1)/2))
  blockDiag <- !diag & !full

  ## CHECK IF ML WAS USED AND NOT REML
  if (lme4::isREML(m1) || lme4::isREML(m0)) stop("Error: the models should be fitted using Maximum Likelihood (ML) instead of Restricted ML (REML)")

  # retrieve names of parameters to identify their order in the FIM and in the cone
  #nbTotParam <- ncol(invfim1)
  if (linmodel){
    invfim1 <- merDeriv::vcov.lmerMod(m1,full=T)
    invfim0 <- merDeriv::vcov.lmerMod(m0,full=T)
  }else{
    invfim1 <- merDeriv::vcov.glmerMod(m1,full=T)
    invfim0 <- merDeriv::vcov.glmerMod(m0,full=T)
  }
  nameParams1 <- dimnames(invfim1)[[2]]
  nameParams0 <- dimnames(invfim0)[[2]]
  paramTested <- !(nameParams1 %in% nameParams0)
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
              dims=list(nbFE1=nbFixEff1,nbFE0=nbFixEff0,nbRE1=nbRanEff1,nbRE0=nbRanEff0,dimSigma=1),
              structGamma=list(diag=diag,full=full,blockDiag=blockDiag)))
}


##' @rdname funcStruct
modelStructsaemix <- function(m1,m0){
  # dimension of the parameters
  nbFixEff0 <- sum(m0@model@fixed.estim>0)
  nbFixEff1 <- sum(m1@model@fixed.estim>0)
  namesRE0 <- m0@model@name.random
  namesRE1 <- m1@model@name.random
  nameVarTested <- namesRE1[!(namesRE1%in%namesRE0)]
  nameVarTested <- gsub("omega2.","",nameVarTested)
  nbRanEff0 <- length(namesRE0)
  nbRanEff1 <- length(namesRE1)
  nbRanEffTest <- nbRanEff1-nbRanEff0

  covStruct1 <- m1@model@covariance.model
  covStruct0 <- m0@model@covariance.model

  diag <- (sum(covStruct1)==sum(diag(covStruct1)))
  full <- (min(covStruct1)==1)
  blockDiag <- !diag & !full

  # dimension of residual error
  if (m1@model@error.model == "combined"){
    dimSigma <- 2
  }else{
    dimSigma <- 1
  }

  if(nbRanEff0==nbRanEff1) stop("Error: there are the same number of random effects in models m0 and m1. Please check the models' formulation.")
  if(min(covStruct1-covStruct0) < 0) stop("Error: the models should be nested, but it seems that some random effects are in m0 but not in m1")

  # retrieve names of parameters to identify their order in the FIM and in the cone
  #nbTotParam <- ncol(invfim1)
  nameFix1 <- c(m1@model@name.fixed[m1@model@fixed.estim>0])

  dd <- expand.grid(rownames(covStruct1),rownames(covStruct1))
  names(dd) <- c("var1","var2")
  dd$names <- paste(dd[,1],dd[,2],sep=".")
  dd$lowertri <- as.vector(lower.tri(covStruct1,diag = T))
  dd$zero1 <- as.vector(covStruct1==0)
  dd$diag <- as.vector(lower.tri(covStruct1,diag = T)) - as.vector(lower.tri(covStruct1,diag = F))
  dd$type <- ifelse(dd$diag==1,"sd","co")
  dd$tested <- as.vector(covStruct0==0)
  dd <- dd[dd$lowertri & !dd$zero1,]
  dd$zero1 <- dd$lowertri <- dd$diag <- NULL

  dd <- rbind(data.frame(names=nameFix1,var1=rep(NA,length(nameFix1)),var2=rep(NA,length(nameFix1)),type=rep("beta",length(nameFix1)),tested=rep(FALSE,length(nameFix1))),dd)

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

    ## CHECK BELOW
    if (nrow(ddco)>0){
      dd$covTested <- NA
      dd$covInBlock <- NA
      for (i in 1:nrow(ddco)){
        dd[dd$type=="co",][i,]$covTested <- dd[dd$type=="co",][i,]$tested*(!(dd$var1[dd$type=="co"][i] %in% nameVarTested) & !(dd$var2[dd$type=="co"][i] %in% nameVarTested))
        dd[dd$type=="co",][i,]$covInBlock <- dd[dd$type=="co",][i,]$tested*(dd$var1[dd$type=="co"][i] %in% nameVarTested) & (dd$var2[dd$type=="co"][i] %in% nameVarTested) # identify covariances that are tested as part of a block of the covariance matrix tested
      }
      dd$covInBlock[dd$type=="sd"] <- 1
    }
  }
  rownames(dd) <- 1:nrow(dd) # reinitialize rownames (used to be indices of elements in the covariance matrix)

  return(list(detailStruct=dd,
              nameVarTested=nameVarTested,
              dims=list(nbFE1=nbFixEff1,nbFE0=nbFixEff0,nbRE1=nbRanEff1,nbRE0=nbRanEff0,dimSigma=dimSigma),
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
