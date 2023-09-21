
<!-- README.md is generated from README.Rmd. Please edit that file -->

[![](https://www.r-pkg.org/badges/version-last-release/varTestnlme)](https://CRAN.R-project.org/package=varTestnlme)
[![](https://cranlogs.r-pkg.org/badges/grand-total/varTestnlme)](https://cranlogs.r-pkg.org/badges/grand-total/varTesnlme)

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

Baey C, Kuhn E, 2023. varTestnlme: An R Package for Variance Components
Testing in Linear and Nonlinear Mixed-Effects Models. **Journal of
Statistical Software**. <https://doi.org/10.18637/jss.v107.i06>

Baey C, Cournède P-H, Kuhn E, 2019. Asymptotic distribution of
likelihood ratio test statistics for variance components in nonlinear
mixed effects models. **Computational Statistic and Data Analysis**.
135:107–122 (2019), <https://doi.org/10.1016/j.csda.2019.01.014>

## Installation

Install from CRAN:

``` r
install.packages("varTestnlme")
```

Or install the development version from Github:

``` r
install.packages("devtools")
devtools::install_github("baeyc/varTestnlme")
```

## Example

An example using the `nlme` package.

**Since version 1.0.0, the name of the main function has been changed
from `varTest` to `varCompTest` due to a conflict with an existing
function from package `EnvStats`.**

``` r
library(nlme)
data("Orthodont")

# using nlme, with correlated slope and intercept
m1 <- lme(distance ~ 1 + Sex + age + age*Sex, random = pdSymm(Subject ~ 1 + age), data = Orthodont, method = "ML")
m0 <- lme(distance ~ 1 + Sex + age + age*Sex, random = ~ 1 | Subject, data = Orthodont, method = "ML")
vt <- varCompTest(m1,m0)
#> Variance components testing in mixed effects models
#> Testing that the variance of the random effect associated to age is equal to 0
#> Likelihood ratio test statistic:
#>  LRT = 0.8331072
#> 
#> p-value from exact weights: 0.5103454
#> 
```

It works similarly with lme4 package or saemix.
