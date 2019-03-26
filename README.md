
<!-- README.md is generated from README.Rmd. Please edit that file -->
varTestnlme
===========

Implements the likelihood ratio test (LRT) for testing the presence of random effects in linear, generalized linear and nonlinear mixed-effects model. The test can be used to answer questions of the type:

-   should a certain subset of random effects be in fact considered as fixed effects?
-   is there any random effects in the model?
-   are there any correlation between two subsets of random effects?

It is possible to compare two models with different random effects, provided that the random structures of the two models are **nested**.

The package works on models that were fitted using `nlme`, `lme4` or `saemix` packages.

Installation
------------

Install the development version from Github:

``` r
# install.packages("devtools")
devtools::install_github("baeyc/varTestnlme")
#> Skipping install of 'varTestnlme' from a github remote, the SHA1 (21414cf5) has not changed since last install.
#>   Use `force = TRUE` to force installation
```

Example
-------

``` r
# First, fit models from your preferred mixed-effects models package
require(nlme)
#> Loading required package: nlme
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
#> model under H1: conc ~ SSbiexp(time, A1, lrc1, A2, lrc2) (model), (A1 ~ 1, lrc1 ~ 1, A2 ~ 1, lrc2 ~ 1) (random effects)
#> model under H0: conc ~ SSbiexp(time, A1, lrc1, A2, lrc2) (model), (A1 ~ 1, lrc1 ~ 1) (random effects)
#> 
#> Likelihood ratio test statistics: 
#>  LRT =  4.560838 
#> 
#> Limiting distribution: 
#>   mixture of 3 chi-bar-square distributions with degrees of freedom 0, 1, 2 
#>   associated weights 0.249, 0.5, 0.251 
#> 
#> p-value: 0.04206452
```
