
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

An example using the `nlme` package.

``` r
library(nlme)
data("Orthodont")

# using nlme, with correlated slope and intercept
m1 <- lme(distance ~ 1 + Sex + age + age*Sex, random = pdSymm(Subject ~ 1 + age), data = Orthodont, method = "ML")
m0 <- lme(distance ~ 1 + Sex + age + age*Sex, random = ~ 1 | Subject, data = Orthodont, method = "ML")

vt <- varTest(m1,m0)
#> Variance components testing in mixed-effects models
#> (models fitted using the nlme package)
#> Testing that the variance of age is null
#> model under H1: distance ~ 1 + Sex + age + age * Sex  (fixed effects),  pdSymm(Subject ~ 1 + age)  (random effects)
#> model under H0: distance ~ 1 + Sex + age + age * Sex  (fixed effects),  ~1 | Subject  (random effects)
#> 
#> Likelihood ratio test statistics: 
#>  LRT =  0.83311 
#> 
#> Limiting distribution:
#> mixture of 2 chi-bar-square distributions with degrees of freedom 1, 2
#> lower-bound for p-value: 0.51035  upper bound for p-value: 0.51035
```

It works similarly with lme4 package or saemix.
