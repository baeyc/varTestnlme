#' Variance components testing in linear and nonlinear mixed effects models
#'
#' \code{varTestnlme} implements likelihood ratio tests for variance components in general (linear and nonlinear) mixed effects models,
#' assuming a multivariate Gaussian distribution for the random effects.
# More preciely, we consider models of the form:
# \deqn{y_{i} = g(\varphi_i, x_{i}) + \varepsilon_{i}
# \varphi_i =U_i \beta + V_i b_i \ \ ,  \ 1 \leq i \leq n,}
# with \eqn{y_{i}} denotes the vector of \eqn{n_i} observations of individual \eqn{i} , \eqn{1 \leq i \leq n}, \eqn{\missing_i} the
# vector of individual parameters of individual \eqn{i}, \eqn{x_{i}} a vector of covariates, and \eqn{\varepsilon_i} an error term.
#'
#' @name varTestnlme-package
#' @aliases varTestnlme-package varTestnlme
#' @docType package
#' @author Charlotte Baey (\email{charlotte.baey@univ-lille.fr})
#' @references Charlotte Baey, Paul-Henry Cournède, Estelle Kuhn (2019). Asymptotic distribution of likelihood ratio test statistics for
#' variance components in nonlinear mixed effects models. \emph{Computational Statistics and Data Analysis.}
#' @keywords package
#' @import methods
#' @noRd
NULL
