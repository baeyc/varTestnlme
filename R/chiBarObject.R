#' Class "chiBarSquareObject"
#'
#' An object of the \code{chiBarSquareObject} class, storing the parameters of the chi-bar-square distribution.
#'
#'
#' @name chiBarSquareObject-class
#' @docType class
#' @aliases chiBarSquareObject-class chiBarSquareObject
#' @section Objects from the Class:
#' An object of the chiBarSquareObject contains the following slots:
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
