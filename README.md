
<!-- README.md is generated from README.Rmd. Please edit that file -->
mnt
===

<!-- badges: start -->
<!-- badges: end -->
The package mnt is designed to give users access to state of the art tests of multivariate normality. It accompanies the survey paper on goodness of fit tests of multivariate normality by Ebner, B. and Henze, N. (2020) Tests for multivariate normality -- a critical review with emphasis on weighted *L*<sup>2</sup>-statistics, that will appear in TEST. All of the described tests can be performed by functions provided in mnt.

Installation
------------

You can install the released version of mnt from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("mnt")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("LBPy/mnt")
```

Example
-------

This is a basic example on how to use the mnt package: We generate a multivariate data set X.data and perform the BHEP test of normality for the generated X.data and using the tuning parameter a=3. The significance level is alpha. Note that the critical values are simulated by a Monte Carlo method.

``` r
library(mnt)
X.data = MASS::mvrnorm(50,c(3,4,5),diag(3,3)) 
X.BHEP = test.BHEP(X.data,a=3,alpha=0.05) 
```

``` r
X.BHEP 
#> 
#> ------------------------------------------------------------------------- 
#> 
#>          Test for multivariate normality with the BHEP  teststatistic.
#> 
#> tuning parameter = 3  
#> BHEP  =  0.9514364  
#> critical value =   1.09841  (via monte carlo) 
#> 
#> 
#> -------------------------------------------------------------------------
```

The value of the test statistic can directly be computed by

``` r
BHEP(X.data,a=3)                       
#> [1] 0.9514364
```

This also works in the univariate case:

``` r
X.data = stats::rnorm(25,3,5)
X.BHEP = test.BHEP(X.data,a=2,alpha=0.05) 
BHEP(X.data,a=2)     
```

``` r
X.BHEP 
#> 
#> ------------------------------------------------------------------------- 
#> 
#>          Test for multivariate normality with the BHEP  teststatistic.
#> 
#> tuning parameter = 2  
#> BHEP  =  0.6427705  
#> critical value =   0.9922879  (via monte carlo) 
#> 
#> 
#> -------------------------------------------------------------------------
```

And for other test statistics too:

``` r
X.data = stats::rnorm(25,3,5)
X.DEHT = test.DEHT(X.data,a=2,alpha=0.05) 
DEHT(X.data,a=2)     
```

``` r
X.DEHT 
#> 
#> ------------------------------------------------------------------------- 
#> 
#>          Test for multivariate normality with the DEH based on harmonic oscillator  teststatistic.
#> 
#> tuning parameter = 2  
#> DEH based on harmonic oscillator  =  1.303027  
#> critical value =   1.464164  (via monte carlo) 
#> 
#> 
#> -------------------------------------------------------------------------
```
