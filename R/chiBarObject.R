#' Class "chiBarSquareObject"
#'
#' An object of the \code{chiBarSquareObject} class, storing the parameters of the chi-bar-square distribution.
#'
#'
#' @name chiBarSquareObject-class
#' @docType class
#' @aliases chiBarSquareObject-class chiBarSquareObject
#' @slot V a positive-definite matrix
#' @slot dims the set of dimensions defining the cone
#' @slot orthant logical, equals \code{TRUE} is the cone is the nonnegative orthant of R^r
chiBarSquareObject <- setClass("chiBarSquareObject",
    slots = c(V = "matrix",
              invV = "matrix",
              dims = "list",
              df = "numeric",
              weights = "numeric",
              orthan = "logical"),

    # Set the default values for the slots
    prototype=list(
      V = diag(1),
      invV = diag(1),
      dims = list(dim0=0,dimR=0,dimSplus=numeric(),dim0=0),
      df = NULL,
      weights = NULL,
      orthan = T
    ),

    validity=function(object)
    {
      if(!is.positive.definite(object@invV)) {
        return("Matrix V is not positive definite")
      }
      return(TRUE)
    }
)


setGeneric("dfchisqbar", function(object){standardGeneric("dfchisqbar")})
setGeneric("pchisqbar", function(q,object,lower.tail){standardGeneric("pchisqbar")})

# Degrees of freedom

#' @name dfchisqbar
#' @aliases dfchisqbar,chiBarSquareObject-method
#' @title dfchisqbar
#' Computes the degrees of freedom of the chi-square involved in the mixture
#' @param object a chiBarSquareObject
#' @rdname dfchisqbar-methods
#' @exportMethod dfchisqbar


setMethod("dfchisqbar",
          c("chiBarSquareObject"),
          function(object){
            # dimLSincluded : dimension of the buggest linear space included in the cone
            # dimLScontaining : dimension of the smallest linear space containing the cone
            q <- sum(unlist(object@dims)) # total dimension
            dimLSincluded <- object@dims$dimBeta$dimR + object@dims$dimGamma$dimR # weights of chi-bar with df=0, ..., dimLSincluded-1 are null
            if (object@orthan) dimLScontaining <- object@dims$dimBeta$dimR + object@dims$dimGamma$dimR + sum(object@dims$dimGamma$dimSplus) # weights of chi-bar with df=dimLScontaining+1, ..., q are null
            if (!object@orthan) dimLScontaining <- object@dims$dimBeta$dimR + object@dims$dimGamma$dimR + sum(object@dims$dimGamma$dimSplus*(object@dims$dimGamma$dimSplus+1)/2)
            return(seq(dimLSincluded,dimLScontaining,1))
          })

# Cumulative distribution function

#' @name pchisqbar
#' @aliases pchisqbar,numeric,chiBarSquareObject,logical-method
#' @title pchisqbar
#' Cumulative distribution function of the chi-bar-square distribution
#' @param q the quantile
#' @param object a chiBarSquareObject
#' @param lower.tail logical, default to TRUE
#' @rdname pchisqbar-methods
#' @exportMethod pchisqbar
#'
setMethod("pchisqbar",
          c("numeric","chiBarSquareObject","logical"),
          function(q,object,lower.tail=T){
            paste(lower.tail)
            paste(q)
            sum(object@weights * stats::pchisq(q,df=object@df,lower.tail = lower.tail))
          }
          )
