#' Monte Carlo approximation of chi-bar-square weights
#'
#' Approximation of the chi-bar-square weigths via Monte Carlo approximation.
#'
#' The chi-bar-square distribution \eqn{\bar{\chi}^2(I,C)} is a mixture of chi-square distributions. The function provides
#' a method to approximate the weights of the mixture components, when the number of components is known as well as the
#' degrees of freedom of each chi-square distribution in the mixture, and given a vector of simulated values from the target
#' \eqn{\bar{\chi}^2(I,C)} distribution. Let us assume that there are \eqn{p} components in the mixture, with degrees of
#' freedom between \eqn{n_1} and \eqn{n_p}. By definition of a mixture distribution, we have :
#' \deqn{ P(\bar{\chi}^2(I,C) \leq c) = \sum_{i=n_1}^{n_p} w_i P(\chi^2_{i} \leq c)}
#' Choosing \eqn{p-2} values \eqn{c_1, \dots, c_{p-2}}, the function will generate a system of \eqn{p-2} equations
#' according to the above relationship, and add two additional relationships stating that the sum of all the weights is
#' equal to 1, and that the sum of odd weights and of even weights is equal to 1/2, so that we end up with a system a \eqn{p}
#' equations with \eqn{p} variables.
#'
#'
#' @name approxWeights
#'
#' @param x a vector of i.i.d. random realizations of the target chi-bar-square distribution
#' @param df a vector containing the degrees of freedom of the chi-squared components
#' @param q the empirical quantile of \code{x} used to choose the \eqn{p-2} values \eqn{c_1, \dots, c_{p-2}} (see Details)
#' @return A vector containing the estimated weights, as well as their covariance matrix.
#' @author Charlotte Baey <\email{charlotte.baey@univ-lille.fr}>
#'
#' @export approxWeights
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
