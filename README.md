
# rsvddpd

<!-- badges: start -->

[![Build
Status](https://www.travis-ci.com/subroy13/rsvddpd.svg?branch=master)](https://www.travis-ci.com/subroy13/rsvddpd)
<!-- badges: end -->

The `R` package `rsvddpd` is an acronym for Robust Singular Value
Decomposition using Density Power Divergence. As the name suggests, the
package mainly concerns with a special function for performing SVD in a
robust way in presence of outliers.

There are 2 primary functions in the package.

-   rSVDdpd - Performs the robust SVD of a numeric matrix *X*.
-   simSVD - Simulates different scenarios to compare performances of
    SVD algorithms under outlier contamination.

## Installation

You can install the development version from
[GitHub](https://github.com/subroy13/rsvddpd) with:

``` r
# install.packages("devtools")
devtools::install_github("subroy13/rsvddpd")
```

Use the following to install the development version with manuals and
vignettes, which provides useful information about the structure of the
function.

``` r
devtools::install_github("subroy13/rsvddpd", build_opts = c("--no-resave-data"), build_manual = TRUE, build_vignettes = TRUE)
```

## Examples

This is a basic example usages which shows the need for the package.

``` r
library(rsvddpd)

X <- matrix(1:20, nrow = 4, ncol = 5)
svd(X)
#> $d
#> [1] 5.352022e+01 2.363426e+00 4.870683e-15 7.906968e-16
#> 
#> $u
#>            [,1]       [,2]        [,3]       [,4]
#> [1,] -0.4430188 -0.7097424 -0.52426094  0.1585890
#> [2,] -0.4798725 -0.2640499  0.81721984  0.1793091
#> [3,] -0.5167262  0.1816426 -0.06165685 -0.8343851
#> [4,] -0.5535799  0.6273351 -0.23130204  0.4964870
#> 
#> $v
#>             [,1]        [,2]       [,3]       [,4]
#> [1,] -0.09654784  0.76855612 -0.6000256  0.1704800
#> [2,] -0.24551564  0.48961420  0.5577664 -0.5560862
#> [3,] -0.39448345  0.21067228  0.2312115  0.1606664
#> [4,] -0.54345125 -0.06826963  0.2643802  0.6650059
#> [5,] -0.69241905 -0.34721155 -0.4533325 -0.4400661
```

As you can see, the first two singular values are 53.5 and 2.36, and the
third and fourth singular values are very small positive reals.

Let us see what happens when you contaminate just one entry of the
matrix by a large value say 100.

``` r
X[2, 3] <- 100
svd(X)
#> $d
#> [1] 1.070340e+02 3.617861e+01 2.200002e+00 1.851858e-15
#> 
#> $u
#>            [,1]       [,2]         [,3]          [,4]
#> [1,] -0.1472125 -0.4893994  0.816938282  2.672612e-01
#> [2,] -0.9548191  0.2971284  0.005940614 -1.110223e-16
#> [3,] -0.1753739 -0.5614500 -0.105644232 -8.017837e-01
#> [4,] -0.1894546 -0.5974754 -0.566935489  5.345225e-01
#> 
#> $v
#>             [,1]       [,2]         [,3]          [,4]
#> [1,] -0.03121244 -0.1097166 -0.798115312  3.357170e-01
#> [2,] -0.08603093 -0.2591083 -0.524844308 -6.516983e-01
#> [3,] -0.94371346  0.3306537 -0.008548386  1.387779e-16
#> [4,] -0.19566792 -0.5578918  0.021697699  6.122270e-01
#> [5,] -0.25048641 -0.7072836  0.294968703 -2.962457e-01
```

Note that, the first singular value changes drastically, being 107,
while second and third singular values 36.1 and 2.2 respectively.
However, such error is very common in practice, and can pose serious
problem in many statistical estimation techniques. `rSVDdpd` solves the
problem as shown in the following code.

``` r
rSVDdpd(X, alpha = 0.3)
#> $d
#> [1] 5.355988e+01 2.358919e+00 1.491240e-01 6.679025e-11
#> 
#> $u
#>           [,1]       [,2]       [,3]          [,4]
#> [1,] 0.4426835 -0.7124329  0.4743859  2.672615e-01
#> [2,] 0.4810552 -0.2588262 -0.8376126  2.694772e-07
#> [3,] 0.5163460  0.1804024  0.2408009 -8.017834e-01
#> [4,] 0.5531763  0.6268200  0.1240084  5.345228e-01
#> 
#> $v
#>            [,1]        [,2]       [,3]          [,4]
#> [1,] 0.09646644  0.77032333  0.2174772 -5.826874e-01
#> [2,] 0.24532587  0.49133606  0.2200082  7.029334e-01
#> [3,] 0.39606258  0.20320671 -0.8954560  5.541357e-07
#> [4,] 0.54304952 -0.06663851  0.2250699  2.219514e-01
#> [5,] 0.69191119 -0.34562580  0.2276008 -3.421953e-01
```

# Author and Maintainer

-   Subhrajyoty Roy, Indian Statistical Institute, Kolkata

# Getting help

If you encounter a clear bug, please file an issue with a minimal
reproducible example on
[GitHub](https://github.com/subroy13/rsvddpd/issues).

------------------------------------------------------------------------

This package is distributed under
[MIT](https://github.com/subroy13/rsvddpd/blob/master/LICENSE.md)
license.
