#' @title Variance component testing
#'
#' @description Perform a likelihood ratio test to test whether a subset of the variances of the random effects
#' are equal to zero. The test is defined by two hypotheses, H0 and H1, and the model under H0 is
#' assumed to be nested within the model under H1. These functions can be used on objects of class
#' \code{lme}-, \code{nlme}-, \code{mer}-, \code{lmerMod}, \code{glmerMod}, \code{nlmerMord} or \code{SaemixObject}.
#'
#' It is possible to tests if any subset of the variances are equal to zero. However, the function does not
#' currently support nested random effects, and assumes that the random effects are Gaussian.
#' 
#' @details The asymptotic distribution of the likelihood ratio test is a chi-bar-square, with weights that need to be
#' approximated by Monte Carlo methods, apart from some specific cases where they are available explicitly. 
#' Therefore, the p-value of the test is not exact but approximated. This computation can be time-consuming, so
#' the default behaviour of the function is to provide bounds on the exact p-value, which can be enough in practice
#' to decide whether to reject or not the null hypothesis. This is triggered by the option \code{pval.comp="bounds"}.
#' To compute an approximation of the exact p-value, one should use the option \code{pval.comp="approx"} or \code{pval.comp="both"}.
#' 
#' When \code{pval.comp="approx"} or \code{pval.comp="both"}, the weights of the chi-bar-square distribution are computed using Monte Carlo,
#' which might involve a larger computing time.
#' 
#' The \code{control} argument controls the options for chi-bar-square weights computation. It is a list with the
#' following elements: \code{M} the size of the Monte Carlo simulation, i.e. the number of samples generated, \code{parallel} a
#' boolean to specify if parallel computing should be used, and \code{nbcores} the number of cores to be used in case of 
#' parallel computing. Default is \code{M=5000}, \code{parallel=FALSE} and \code{nb_cores=1}.
#' If \code{parallel=TRUE} but the value of \code{nb_cores} is not given, then it is set to the number of detected
#' cores minus 1
#'
#' @name varCompTest
#' @aliases varCompTest.lme varCompTest.lme4 varCompTest.saemix
#'
#' @param  m1 a fit of the model under H1, obtained from \code{nlme}, \code{lme4}
#' or \code{saemix}
#' @param m0 a fit of the model under H0, obtained from the same package as \code{m0}
#' @param control (optional) a list of control options for the computation of the chi-bar-weights (see Details section)
#' @param pval.comp (optional) the method to be used to compute the p-value, one of: \code{"bounds"} (the default),
#' \code{"approx"} or \code{"both"} (see Details section)
#' @param fim (optional) the method to compute the Fisher Information Matrix. Options are: \code{fim="extract"} to extract the
#' FIM computed by the package which was used to fit the models, \code{fim="compute"} to evaluate the FIM using parametric
#' bootstrap, and \code{fim=I} with \code{I} a positive semidefinite matrix, for a FIM provided by the user.
#' @param  output a boolean specifying if any output should be printed in the console (default to TRUE) 
#' 
#' @return An object of class \code{htest} with the following components:
#' \itemize{
#' \item \code{statistic} the likelihood ratio test statistics
#' \item \code{null.value} 
#' \item \code{alternative} 
#' \item \code{parameters} the parameters of the limiting chi-bar-square distribution: the degrees of freedom and
#' the weights of the chi-bar-square components and the Fisher Information Matrix
#' \item \code{method} a character string indicating the name of the test
#' \item \code{pvalue} a named vector containing the different p-values computed by the function: using the 
#' (estimated) weights, using the random sample from the chi-bar-square distribution, and the two bounds on
#' the p-value.
#' }
#'
#' @author Charlotte Baey <\email{charlotte.baey@univ-lille.fr}>
#'
#' @examples 
#' # load lme4 package and example dataset
#' library(lme4)
#' data(Orthodont, package = "nlme")
#' 
#' # fit the two models under H1 and H0
#' m1 <- lmer(distance ~ 1 + Sex + age + age*Sex + 
#' (0 + age | Subject), data = Orthodont, REML = FALSE)
#' m0 <- lm(distance ~ 1 + Sex + age + age*Sex, data = Orthodont)
#' 
#' # compare them (order is important: m1 comes first)
#' varCompTest(m1,m0,pval.comp="bounds")
#' 
#' # using nlme
#' library(nlme)
#' m1 <- lme(distance ~ 1 + Sex + age + age*Sex, 
#' random = pdSymm(Subject ~ 1 + age), data = Orthodont, method = "ML")
#' m0 <- lme(distance ~ 1 + Sex, random = ~ 1 | Subject, data = Orthodont, method = "ML")
#' 
#' varCompTest(m1,m0)
#' 
#' @references Baey C, Cournède P-H, Kuhn E, 2019. Asymptotic distribution of likelihood ratio test
#' statistics for variance components in nonlinear mixed effects models. \emph{Computational
#' Statistics and Data Analysis} 135:107-122.
#'
#' Silvapulle  MJ, Sen PK, 2011. Constrained statistical inference: order, inequality and shape constraints.
#' @export varCompTest
#' @importFrom stats formula pchisq
varCompTest <- function(m1, m0, control = list(M=5000,parallel=T,nb_cores=1,B=1000),pval.comp = "bounds",fim = "extract", output=TRUE) {

  pkg1 <- pckName(m1)
  pkg0 <- pckName(m0)
  
  samePkg <- max(pkg1 == pkg0)
  
  if(!(pkg1 %in% c("nlme", "lme4", "saemix")))
    stop("'m1' must be fitted using 'nlme', 'lme4', or 'saemix' package")
  if(!samePkg & !max(pkg0 %in% c("lm","nls","glm")))
    stop("'m1' and 'm0' must be fitted with the same package")
  
  UseMethod("varCompTest", m1)
}

varTest <- function(m1, m0, control = list(M=5000,parallel=T,nb_cores=1,B=1000),pval.comp = "bounds",fim = "extract") {
  .Deprecated("varCompTest")
  
  varCompTest(m1, m0, control = list(M=5000,parallel=T,nb_cores=1,B=1000),pval.comp = "bounds",fim = "extract", output=TRUE)
}


#' @title Extracting models' structures
#'
#' @description Functions extracting the structure of the models under both hypothesis: the number of fixed and random effects, 
#' the number of tested fixed and random effects, and the residual dimension, as well as the random effects covariance structure
#'
#' @param m1 the model under H1
#' @param m0 the model under H0
#' @param randm0 a boolean stating whether the model under H0 contains any random effect

#' @return A list with the following components:
#' \item{\code{detailStruct}}{a data frame containing the list of the parameters and whether they are tested or not}
#' \item{\code{nameVarTested}}{the name of the variance components being tested}
#' \item{\code{nameFixedTested}}{the name of the fixed effects being tested}
#' \item{\code{dims}}{a list with the dimensions of fixed and random effects, tested or not tested}
#' \item{\code{structGamma}}{the structure of the covariance matrix of the random effects \code{diag}, \code{full} or 
#' \code{blockDiag}}
#' 
#' @name extractStruct
#' @export
extractStruct <- function(m1,m0,randm0){
  UseMethod("extractStruct",m1)
}


#' @title Extract covariance matrix 
#'
#' @description Extract covariance matrix of the random effects for a model fitted with lme4.
#'
#' @param m a fit from lme4 package (either linear or nonlinear)
#' 
#' @export extractVarCov
extractVarCov <- function(m){
  UseMethod("extractVarCov",m)
}


#' @title Approximation of the inverse of the Fisher Information Matrix via parametric bootstrap
#'
#' @description When the FIM is not available, this function provides an approximation of the FIM based on an estimate
#' of the covariance matrix of the model's parameters obtained via parametric bootstrap.
#'
#' @name bootinvFIM
#' @aliases invFIM bootstrap
#'
#' @param m a fitted model that will be used as the basis of the parametric bootstrap (providing the initial maximum
#' likelihood estimate of the parameters and the modelling framework)
#' @param B the size of the bootstrap sample
#' @return the empirical covariance matrix of the parameter estimates obtained on the bootstrap sample
#' @author Charlotte Baey <\email{charlotte.baey@univ-lille.fr}>
#'
#' @importFrom foreach %dopar%
#' @export bootinvFIM
bootinvFIM <- function(m,B=1000){
  UseMethod("bootinvFIM",m)
}
