#' Class "chibarObject"
#'
#' An object of the \code{chibarObject} class, storing the
#'
#' Details of the algorithm can be found in the pdf file accompanying the package.
#'
#' @name chibarObject-class
#' @docType class
#' @aliases chibarObject-class chibarObject
#' @section Objects from the Class:
#' An object of the chibarObject contains the following slots:
#' @slot V a positive-definite matrix
#' @slot dims the set of dimensions defining the cone
#' @slot orthant logical, equals \code{TRUE} is the cone is the nonnegative orthant of R^r
chiBarObject <- setClass("chiBarObject",
    slots = c(V = "matrix",
              invV = "matrix",
              dims = "list",
              orthan = "logical"),

    # Set the default values for the slots
    prototype=list(
      V = diag(1),
      invV = diag(1),
      dims = list(dim0=0,dimR=0,dimSplus=numeric(),dim0=0),
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
