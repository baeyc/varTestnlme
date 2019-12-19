#' @name plot.varTestObject
#' @rdname plot.varTestObject
#'
#' @title Diagnostic plot for the approximation of the chi-bar-square distribution
#'
#' @description Plot the empirical cumulative distribution function (cdf) of the simulated chi-bar-square distributed variable, along with the exact cdf
#' of all the chi-square distributions involved in the mixture, and with the cdf based on the approximated weights. This function can only be used when
#' the weights were approximated by simulation.
#'
#' @param x a object of class \code{\link{varTestObject}} obtained from a call to function \code{\link{varTest}}
#'
#' @importFrom stats pchisq ecdf
#' @importFrom graphics plot points legend
NULL

plot.varTestObject <- function(x){
  if (length(x@chibarsquare)==0) stop("Only available when the weights were computed by simulation")

  if (requireNamespace("ggplot2", quietly = TRUE)){
    cbs <- y1 <- df <- NULL # solve "no visible binding for global variable" issue
    d <- data.frame(cbs=x@chibarsquare)
    x1 <- seq(0,max(x@chibarsquare),0.1)
    y <- as.vector(sapply(x@cbs@df,FUN = function(df){pchisq(x1,df)}))
    d2 <- data.frame(x1=rep(x1,length(x@cbs@df)),y1=y,df=rep(x@cbs@df,each=length(x1)))
    ymix <- 0*x1
    for (i in 1:length(x@cbs@df)){
      ymix <- ymix + x@weights[i]*pchisq(x1,x@cbs@df[i])
    }
    d3 <- data.frame(x1=x1,y=ymix)
    cols <- c("Components"="darkgrey","Mixture"="black")
    line_types <- c("Plain"=1,"Points"=3)
    ggplot2::ggplot(data=d,ggplot2::aes(x=cbs,colour="Mixture",linetype="Plain")) + ggplot2::stat_ecdf() + ggplot2::geom_line(data=d2,ggplot2::aes(x=x1,y=y1,group=df,colour="Components",linetype="Plain")) +
      ggplot2::geom_line(data=d3,ggplot2::aes(x=x1,y=y,linetype = "Points",colour="Mixture")) + ggplot2::scale_color_manual(values=cols,labels=c("CDF of each mixture component","Empirical CDF"),name="") +
      ggplot2::scale_linetype_manual(breaks=c("Points"),values=line_types,name="",labels=c("CDF of chi-bar-square using weights estimates")) +
      ggplot2::guides(color = ggplot2::guide_legend(order = 1), linetype = ggplot2::guide_legend(order = 2)) +  ggplot2::theme(legend.position = 'bottom',legend.spacing.x = ggplot2::unit(0.25, 'cm')) + ggplot2::ggtitle("Estimation of chi-bar-square distribution")
  }else{
    plot(ecdf(x@chibarsquare),cex=0.25,main="Estimation of chi-bar-square distribution")
    x1 <- seq(0,max(x@chibarsquare),0.1)
    ymix <- 0*x1
    for (i in 1:length(x@cbs@df)){
      points(x1,pchisq(x1,x@cbs@df[i]),type="l",col="darkgrey")
      ymix <- ymix + x@weights[i]*pchisq(x1,x@cbs@df[i])
    }
    points(x1,ymix,type="l",lty=2)
    legend("bottomright",c("Empirical cdf","CDF of chi-bar-square using weights estimates","CDF of each mixture component"),lty=c(1,2,1),col=c("black","black","darkgrey"),bty="n")
  }
}

