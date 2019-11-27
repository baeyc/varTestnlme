#' Class "varTestObject"
#'
#' An object of the \code{varTestObject} class, storing the results of the LRT
#'
#'
#' @name varTestObject-class
#' @docType class
#' @aliases varTestObject-class varTestObject
#' @section Objects from the Class:
#' An object of the varTestObject contains the following slots:
#' @slot lrt the likelihood ratio test statistics
#' @slot df the degrees of freedom of the chi-square distributions involved in the mixture
#' @slot weights the weights associated to the chi-square distributions involved in the mixture
#' @slot pvalue the p-value of the LRT
varTestObject <- setClass("varTestObject",
         slots = c(lrt = "numeric",
                   cbs = "chiBarSquareObject",
                   namesTestedParams = "list",
                   chibarsquare = "numeric",
                   weights = "numeric",
                   sdWeights = "numeric",
                   pvalueWeights = "numeric",
                   pvalueMC = "numeric"),

         validity=function(object)
         {
           return(TRUE)
         }
)

