
<!-- README.md is generated from README.Rmd. Please edit that file -->

# varTestnlme

[varTestnlme](https://baeyc.github.io/varTestnlme/index.html) implements
the likelihood ratio test (LRT) for testing the presence of random
effects in linear, generalized linear and nonlinear mixed-effects model.
The test can be used to answer questions of the type:

  - should a certain subset of random effects be in fact considered as
    fixed effects?
  - is there any random effects in the model?
  - are there any correlation between two subsets of random effects?

It is possible to compare two models with different random effects,
provided that the random structures of the two models are **nested**.

The package works on models that were fitted using
[nlme](https://CRAN.R-project.org/package=nlme),
[lme4](https://CRAN.R-project.org/package=lme4), or
[saemix](https://CRAN.R-project.org/package=saemix) packages.

## Reference

Baey C, Cournède P-H, Kuhn E, 2019. Asymptotic distribution of
likelihood ratio test statistics for variance components in nonlinear
mixed effects models.  135:107–122 (2019),
<https://doi.org/10.1016/j.csda.2019.01.014>

## Installation

Install the development version from Github:

``` r
# install.packages("devtools")
devtools::install_github("baeyc/varTestnlme")
```

## Example

``` r
# First, fit models from your preferred mixed-effects models package
library(nlme)
fm1Indom.lis <- nlsList(conc ~ SSbiexp(time, A1, lrc1, A2, lrc2),data = Indometh)
m1 <- nlme(fm1Indom.lis,random = pdDiag(A1 + lrc1 + A2 + lrc2 ~ 1))
m0 <- nlme(fm1Indom.lis,random = pdDiag(A1 + lrc1 ~ 1))

# Then, simply use the varTest function to compare the two models
# It will automatically detect the covariance structures
varTest(m1,m0)
#> Variance components testing in mixed-effects models
#> (models fitted using the nlme package)
#> 
#> Testing that the variances of A2 and lrc2 are null
#> model under H1: list(A1 ~ 1, lrc1 ~ 1, A2 ~ 1, lrc2 ~ 1) (fixed effects) , structure(list(Subject = structure(c(3.43637993576256, 2.44739674232125, 2.05693983560385, 1.96340665027553), formula = structure(list(A1 ~ 1, lrc1 ~ 1, A2 ~ 1, lrc2 ~ 1), class = "listForm"), Dimnames = list(c("A1", "lrc1", "A2", "lrc2"), c("A1", "lrc1", "A2", "lrc2")), class = c("pdDiag", "pdMat"))), settings = c(0, 1, 0, 1), class = "reStruct") (random effects)
#> model under H0: list(A1 ~ 1, lrc1 ~ 1, A2 ~ 1, lrc2 ~ 1) (fixed effects) , structure(list(Subject = structure(c(3.43637993576256, 2.44739674232125), formula = structure(list(A1 ~ 1, lrc1 ~ 1), class = "listForm"), Dimnames = list(c("A1", "lrc1"), c("A1", "lrc1")), class = c("pdDiag", "pdMat"))), settings = c(0, 1, 0, 1), class = "reStruct") (random effects)
#> 
#> Likelihood ratio test statistics: 
#>  LRT =  4.5608 
#> 
#> Limiting distribution: 
#> mixture of 3 chi-bar-square distributions with degrees of freedom 0, 1, 2 
#> 
#> lower-bound for p-value: 0.016356  upper bound for p-value: 0.067476
```
