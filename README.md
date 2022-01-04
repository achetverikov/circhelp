
<!-- README.md is generated from README.Rmd. Please edit that file -->

# circhelp

<!-- badges: start -->

A small helper package for circular data analyses in R, particularly
useful for cognitive studies on orientation, motion direction, and other
circular features. Contains functions for descriptive statistics for
circular data (computing means, SD, and skewness), angular differences,
and correlation. Also includes a function to correct for cardinal biases
in the human estimates of circular features (e.g., orientation).
<!-- badges: end -->

The goal of circhelp is to …

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("achetverikov/circhelp")
```

## Usage

Most of the functions are self-explanatory.

``` r
library(circhelp)
#> Loading required package: gamlss
#> Loading required package: splines
#> Loading required package: gamlss.data
#> 
#> Attaching package: 'gamlss.data'
#> The following object is masked from 'package:datasets':
#> 
#>     sleep
#> Loading required package: gamlss.dist
#> Loading required package: MASS
#> Loading required package: nlme
#> Loading required package: parallel
#>  **********   GAMLSS Version 5.2-0  **********
#> For more on GAMLSS look at https://www.gamlss.com/
#> Type gamlssNews() to see new features/changes/bug fixes.
#> Loading required package: data.table
#> Loading required package: mgcv
#> This is mgcv 1.8-33. For overview type 'help("mgcv-package")'.
library(mgcv)
# compute a set of descriptive statistics
x <- rnorm(500)
circ_descr(x)
#> $mu
#> [1] 5.045977e-05
#> 
#> $sigma
#> [1] 0.9727528
#> 
#> $skew_pewsey
#> [1] 0.01379079
#> 
#> $skew_fischer
#> [1] 0.05952331
#> 
#> $rho
#> [1] 0.6230528
#> 
#> $skew_rel_to_zero
#> [1] 0.01380606

# compute difference in orientations
a <- 5 
b <- 170
angle_diff_180(a, b)
#> [1] 15

# compute difference in 360° space (e.g., motion directions)
angle_diff_360(a, b)
#> [1] -165

# compute correlation between angles
data <- rmvn(10000, c(0,0), V = matrix(c(1,0.5,0.5,1), ncol = 2))
circ_corr(data[,1], data[,2])
#> [1] 0.4393806
```

The only (somewhat) complicated function is `remove_cardinal_biases`,
see the [help
files](https://achetverikov.github.io/circhelp/reference/remove_cardinal_biases.html)
and [the
vignette](https://achetverikov.github.io/circhelp/articles/cardinal_biases.html)
for an example use case.
<img src="https://achetverikov.github.io/circhelp/articles/cardinal_biases_files/figure-html/correct_biases_ex-1.png" title="Example of cardinal biases processing" width="100%"/>
