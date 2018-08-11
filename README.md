**spef**
--------

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active) [![minimal R version](https://img.shields.io/badge/R%3E%3D-3.4.4-6666ff.svg)](https://cran.r-project.org/) [![Travis-CI Build Status](https://travis-ci.org/stc04003/spef.svg?branch=master)](https://travis-ci.org/stc04003/spef) [![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/stc04003/spef?branch=master&svg=true)](https://ci.appveyor.com/project/stc04003/spef) [![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/spef)](https://cran.r-project.org/package=spef) [![packageversion](https://img.shields.io/badge/Package%20version-1.0.8-orange.svg?style=flat-square)](commits/master) [![Last-changedate](https://img.shields.io/badge/last%20change-2018--08--11-yellowgreen.svg)](/commits/master) ----

Panel count data arise in many applications when the event history of a recurrent event process is only examined at a sequence of discrete time points. Methods in panel count data implemented in this package are grouped into two categories depending on whether the examination times are associated with the recurrent event process after conditioning on covariates. Estimating procedures include: Wang-Yan’s augmented estimating equations (`AEE`, `AEEX`), HuangWang-Zhang’s method (`HWZ`), Zhang’s maximum pseudolikelihood (`MPL`), Maximum pseudolikelihood with I-Splines (`MPLs`), Maximum likelihood with I-Splines (`MLs`), Sun-Wei’s method (`EE.SWa`, `EE.SWb`, `EE.SWc`), Hu-Sun-Wei's method (`EE.HSWc`, `EE.HSWm`), and accelerated mean model (`AMM`).

A comprehensive review can be found in <https://doi.org/10.1111/insr.12271>.

### Installation

You can install spef from CRAN using

``` r
install.packages("spef")
library(spef)
```

You can install spef from GitHub with:

``` r
## install.packages("devtools")
devtools::install_github("stc04003/spef")
```

### Reference

Sun, J. and Wei, L. J. (2000). Regression analysis of panel count data with covariates-dependent observation and censoring times. *Journal of the Royal Statistical Society, Series B: Statistical Methodology*, **62**(2), 293--302.

Zhang, Y. (2002). A Semiparametric pseudolikelihood estimation method for panel count data. *Biometrika*, **89**(1), 39--48.

Huang, C., Wang, M., and Zhang, Y. (2006). Analyzing panel count data with informative observation times. *Biometrika*, **93**(4), 763--776.

Lu, M., Zhang, Y., and Huang, J. (2007). Estimation of the mean function with panel count data using monotone polynomial splines. *Biometrika* **94**(3), 705--718.

Wang, X. and Yan, J. (2011). Fitting semiparametric regressions for panel count survival data with an R package spef. *Computer methods and programs in biomedicine* **104**(2), 278--285.

Wang, X. and Yan, J. (2013). Augmented estimating equations for semiparametric panel count regression with informative observation times and censoring time. *Statistica Sinica*, **23**(1), 359--381.

Chiou, S., Huang, C.-Y., Xu, G., and Yan, J. (2018). Semiparametric regression analysis of panel count data: A practical review. *International Statistical Review*, <https://doi.org/10.1111/insr.12271>.
