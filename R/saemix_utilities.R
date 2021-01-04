#' @name extractStruct.saemix
#' @rdname extractStruct.saemix
#'
#' @title Extract model structure
#'
#' @param m1 the fit under H1
#' @param m0 the fit under H0
#' @param randm0 a boolean indicating whether random effects are present in m0
extractStruct.saemix <- function(m1,m0,randm0){
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
  
  if (!prod(namesFE0 %in% namesFE1)) stop("Error: the models should be nested, but it seems that some fixed effects are in <m0> but not in <m1>")
  if (!prod(namesRE0 %in% namesRE1)) stop("Error: the models should be nested, but it seems that some random effects are in <m0> but not in <m1>")
  
  covStruct1 <- m1@model@covariance.model
  if (randm0) covStruct0 <- m0@model@covariance.model
  
  diag <- (sum(covStruct1)==sum(diag(covStruct1)))
  full <- (min(covStruct1)==1)
  blockDiag <- !diag & !full
  struct <- c("diag","full","blockDiag")[c(diag,full,blockDiag)]
  
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
              structGamma=struct))
}


#' @name bootinvFIM.saemix
#' @rdname bootinvFIM.saemix
#'
#' @title Compute the inverse of the Fisher Information Matrix using parametric bootstrap
#'
#' @param m the model under H1
#' @param B the bootstrap sample size
bootinvFIM.saemix <- function(m, B=1000){
  
  simul <- saemix::simul.saemix(m,nsim=B)
  m.data <- m@data
  options <- list(save=FALSE,save.graphs=FALSE,ll.is=FALSE,save=FALSE,ll.qg=FALSE,print=FALSE,warnings=FALSE,displayProgress=FALSE)
  
  diagSigma <- Matrix::isDiagonal(m@results@omega)
  
  fit.saemix <- lapply(1:B, FUN = function(i){
    data_i <- simul@sim.data@datasim[simul@sim.data@datasim$irep==i,]
    m.data@data[,m.data@name.response] <- data_i$ysim
    log <- utils::capture.output({
      res <- saemix::saemix(m@model,m.data,options)})
    
    if (diagSigma){
      theta <- c(res@results@betas,diag(res@results@omega),res@results@respar[res@results@respar>0])
    }else{
      theta <- c(res@results@betas,res@results@omega[lower.tri(res@results@omega,diag=T)],res@results@respar[res@results@respar>0])
    }
    return(theta)
  })  
  
  theta <- do.call(rbind,fit.saemix)
  
  invfim <- cov(theta)
  colnames(invfim) <- rownames(invfim) <- c(m@results@name.fixed,m@results@name.random,m@results@name.sigma[m@results@respar>0])
  
  return(invfim)
}
