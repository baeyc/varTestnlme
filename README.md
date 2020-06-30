
<!-- README.md is generated from README.Rmd. Please edit that file -->

[![](https://www.r-pkg.org/badges/version/varTestnlme)](https://CRAN.R-project.org/package=varTestnlme)

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

Install from CRAN:

``` r
install.packages("varTestnlme")
```

Install the development version from Github:

``` r
# install.packages("devtools")
#devtools::install_github("baeyc/varTestnlme")

devtools::load_all()
```

## Example

### With default options

The first example is based on the `nlme` package.

``` r
# First, fit models from your preferred mixed-effects models package
library(nlme)
fm1.lis <- nlsList(rate ~ SSasympOff(pressure, Asym, lrc, c0), data=Dialyzer)
nm1e <- nlme(fm1.lis, random = pdSymm(Asym + lrc ~ 1))
nm0e <- nlme(fm1.lis, random = pdDiag(lrc ~ 1))

# Then, simply use the varTest function to compare the two models
# It will automatically detect the covariance structures
vte <- varTest(nm1e, nm0e, fim = "compute", pval.comp = "both", control=list(B=10))
#> Variance components testing in mixed-effects models
#> (models fitted using the nlme package)
#> Testing that the variance of Asym is null
#> model under H1: rate ~ SSasympOff(pressure, Asym, lrc, c0)  (non linear model),  structure(list(Subject = structure(c(2.04719313529589, 0, -1.19146681079193  (random effects)model under H1: rate ~ SSasympOff(pressure, Asym, lrc, c0)  (non linear model),  ), formula = structure(list(Asym ~ 1, lrc ~ 1), class = "listForm"), Dimnames = list(  (random effects)model under H1: rate ~ SSasympOff(pressure, Asym, lrc, c0)  (non linear model),      c("Asym", "lrc"), c("Asym", "lrc")), class = c("pdSymm",   (random effects)model under H1: rate ~ SSasympOff(pressure, Asym, lrc, c0)  (non linear model),  "pdMat"))), settings = c(0, 1, 0, 0), class = "reStruct")  (random effects)
#> model under H1: rate ~ SSasympOff(pressure, Asym, lrc, c0)  (non linear model),  structure(list(Subject = structure(-1.19146681079193, formula = structure(list(  (random effects)model under H1: rate ~ SSasympOff(pressure, Asym, lrc, c0)  (non linear model),      lrc ~ 1), class = "listForm"), Dimnames = list("lrc", "lrc"), class = c("pdDiag",   (random effects)model under H1: rate ~ SSasympOff(pressure, Asym, lrc, c0)  (non linear model),  "pdMat"))), settings = c(0, 1, 0, 1), class = "reStruct")  (random effects)
#> 
#> Likelihood ratio test statistics: 
#>  LRT =  151.85 
#> 
#> Limiting distribution:
#> mixture of 2 chi-bar-square distributions with degrees of freedom 1, 2
#>  associated weights and sd:  0.5 (0), 0.5 (0) 
#> 
#> p-value (from estimated weights): 5.6523e-34
#> lower-bound for p-value: 5.6523e-34  upper bound for p-value: 5.6523e-34
```

It works similarly with lme4 package or saemix.

### Estimation of Fisher Information Matrix

The imiting distribution of the test statistic is a chi-bar-square
distribution whose weights depend on the Fisher Information Matrix (FIM)
of the model. Depending on the package used to fit the model and on your
preferences, it is possible to specify different options with the
argument `fim`: \* `fim="extract"`: (the default option) to extract the
FIM computed by the package used to fit the models. Please note that it
might not be available in some cases. For example, for nonlinear models
fitted with `lme4`. \* `fim="compute"`: to estimate the FIM using
parametric bootstrap. Note that it can take some time to run. In this
case, one must provide the bootstrap sample size with the `control`
argument \* `fim=I`: with `I` a symmetric positive definite matrix, to
provide your own FIM

``` r
# Bootstrap estimation of the FIM
varTest(nm1e, nm0e, fim="compute", pval.comp="both", control=list(B=100))
#> Variance components testing in mixed-effects models
#> (models fitted using the nlme package)
#> Testing that the variance of Asym is null
#> model under H1: rate ~ SSasympOff(pressure, Asym, lrc, c0)  (non linear model),  structure(list(Subject = structure(c(2.04719313529589, 0, -1.19146681079193  (random effects)model under H1: rate ~ SSasympOff(pressure, Asym, lrc, c0)  (non linear model),  ), formula = structure(list(Asym ~ 1, lrc ~ 1), class = "listForm"), Dimnames = list(  (random effects)model under H1: rate ~ SSasympOff(pressure, Asym, lrc, c0)  (non linear model),      c("Asym", "lrc"), c("Asym", "lrc")), class = c("pdSymm",   (random effects)model under H1: rate ~ SSasympOff(pressure, Asym, lrc, c0)  (non linear model),  "pdMat"))), settings = c(0, 1, 0, 0), class = "reStruct")  (random effects)
#> model under H1: rate ~ SSasympOff(pressure, Asym, lrc, c0)  (non linear model),  structure(list(Subject = structure(-1.19146681079193, formula = structure(list(  (random effects)model under H1: rate ~ SSasympOff(pressure, Asym, lrc, c0)  (non linear model),      lrc ~ 1), class = "listForm"), Dimnames = list("lrc", "lrc"), class = c("pdDiag",   (random effects)model under H1: rate ~ SSasympOff(pressure, Asym, lrc, c0)  (non linear model),  "pdMat"))), settings = c(0, 1, 0, 1), class = "reStruct")  (random effects)
#> 
#> Likelihood ratio test statistics: 
#>  LRT =  151.85 
#> 
#> Limiting distribution:
#> mixture of 2 chi-bar-square distributions with degrees of freedom 1, 2
#>  associated weights and sd:  0.5 (0), 0.5 (0) 
#> 
#> p-value (from estimated weights): 5.6523e-34
#> lower-bound for p-value: 5.6523e-34  upper bound for p-value: 5.6523e-34
```

Another example using the `lme4` package and the Orange trees data.

``` r
library(lme4)
#> Le chargement a nécessité le package : Matrix
#> 
#> Attachement du package : 'lme4'
#> The following object is masked from 'package:nlme':
#> 
#>     lmList
dial.start <- unlist(nm1e$call$start)
names(dial.start) <- c("Asym","lrc","c0")
dial.start  <- c(Asym = 50, lrc = 0.5, c0 = 0.2)
nm1  <- nlmer(rate ~ SSasympOff(pressure, Asym, lrc, c0) ~ (0 + Asym + lrc || Subject), Dialyzer, start = dial.start)
nm0  <- nlmer(rate ~ SSasympOff(pressure, Asym, lrc, c0) ~ (0 + Asym | Subject), Dialyzer, start = dial.start)
vt <- varTest(nm1,nm0,fim="compute",pval.comp = "both",control=list(M=10,B=20))
#> Variance components testing in mixed-effects models
#> (models fitted using the lme4 package)
#> Testing that the variance of lrc is null
#> model under H1: rate ~ SSasympOff(pressure, Asym, lrc, c0) ~ (0 + Asym + lrc || Subject)
#> model under H0: rate ~ SSasympOff(pressure, Asym, lrc, c0) ~ (0 + Asym | Subject)
#> 
#> Likelihood ratio test statistics: 
#>  LRT =  28.962 
#> 
#> Limiting distribution:
#> mixture of 2 chi-bar-square distributions with degrees of freedom 0, 1
#>  associated weights and sd:  0.5 (0), 0.5 (0) 
#> 
#> p-value (from estimated weights): 3.6907e-08
#> lower-bound for p-value: 3.6907e-08  upper bound for p-value: 3.6907e-08


library(nlme)
nm1e <- nlme(rate ~ SSasympOff(pressure, Asym, lrc, c0),
             random = pdDiag(Asym + lrc + c0 ~ 1),
             fixed = Asym + lrc + c0 ~ 1,
             data = Dialyzer,
             start = dial.start)
nm0e <- nlme(rate ~ SSasympOff(pressure, Asym, lrc, c0),
             random = pdDiag(Asym + c0 ~ 1),
             fixed = Asym + lrc + c0 ~ 1,
             data = Dialyzer,
             start = dial.start)
varTest(nm1e,nm0e)
#> Variance components testing in mixed-effects models
#> (models fitted using the nlme package)
#> Testing that the variance of lrc is null
#> model under H1: rate ~ SSasympOff(pressure, Asym, lrc, c0)  (non linear model),  pdDiag(Asym + lrc + c0 ~ 1)  (random effects)
#> model under H1: rate ~ SSasympOff(pressure, Asym, lrc, c0)  (non linear model),  pdDiag(Asym + c0 ~ 1)  (random effects)
#> 
#> Likelihood ratio test statistics: 
#>  LRT =  28.604 
#> 
#> Limiting distribution:
#> mixture of 2 chi-bar-square distributions with degrees of freedom 0, 1
#> lower-bound for p-value: 4.4405e-08  upper bound for p-value: 4.4405e-08
```
