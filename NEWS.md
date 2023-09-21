# varTestnlme 1.3.5

Add reference (doi number) for the associated JSS paper.

# varTestnlme 1.3.4

Add an export for function extractFIM.lme 

# varTestnlme 1.3.2

Correct a bug in the comparison of residual variance structure, and coming back to lmersampler for parametric bootstrap in linear mixed model after the package being un-archived.

# varTesnlme 1.3.1

Change the computation of parametric bootstrap for linear mixed models fitted with nlme after lmersampler package was archived

# varTestnlme 1.3.0

Correct a bug for Windows platforms, where new control arguments in varCompTest function could not e accounted for.

# varTestnlme 1.2.0

Correct a typo in the definition of generic functions for SaemixObject.

# varTestnlme 1.1.1

Due to the updated version of package lmeresampler, some bugs appeared in the previous version and were corrected. More precisely, the mySumm function used for the parametric bootstrap now outputs named vectors to avoid bugs when using dplyr::bind_rows function inside lmeresampler package.

# varTestnlme 1.1.0

Correct some typos in documentation, and export print and summary functions. Also export some functions needed to compute the Fisher Information Matrix. Small improvements in the outputs of print and summary functions. Also add some warning when the residual variance structure is not the same between the two models, when using nlme.
Correct a bug in chi-bar-square weights estimation.

# varTestnlme 1.0.0

Change the structure of the package for a more modular one. Now the main function returns an object of class htest, which can be printed for example via the print() function of package EnvStats. To avoid conflicts with the varTest function of EnvStats package, the main function of our package was renamed from varTest to varCompTest.
Parametric bootstrap for the estimation of the FIM is now available for the saemix package also.

# varTestnlme 0.2.0

Add the computation of the FIM via parametric bootstrap for models fitted with nlme and lme4 packages

# varTestnlme 0.1.0

This is the first release of varTestnlme
