# circhelp

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/circhelp)](https://CRAN.R-project.org/package=circhelp)
<!-- badges: end -->

`circhelp` is a small helper package for circular data analyses in R,
particularly useful for cognitive studies on orientation, motion
direction, and other circular features. The package contains functions
for descriptive statistics for circular data (computing means, SD, and
skewness), angular differences, and correlation. It also includes a
function to correct for cardinal biases in the human estimates of
circular features (e.g., orientation). The documentation is available
from the package help files and at
<https://achetverikov.github.io/circhelp/>.

## Installation

You can install the latest released version from CRAN with:

``` r
install.packages("circhelp")
```

The current developmental version from [GitHub](https://github.com/)
could be installe with:

``` r
# install.packages("devtools")
devtools::install_github("achetverikov/circhelp")
```

## Usage

Most functions are self-explanatory.

``` r
library(circhelp)
#> Loading required package: data.table
#> Loading required package: ggplot2
```

``` r
library(mgcv)
#> Loading required package: nlme
#> This is mgcv 1.9-1. For overview type 'help("mgcv-package")'.
```

``` r
# compute a set of descriptive statistics
x <- rnorm(500)
circ_descr(x)
#> $mu
#> [1] -0.00686282
#> 
#> $sigma
#> [1] 0.9944722
#> 
#> $skew_pewsey
#> [1] -0.01750534
#> 
#> $skew_fischer
#> [1] -0.06237638
#> 
#> $rho
#> [1] 0.6098834
#> 
#> $skew_rel_to_zero
#> [1] -0.01980849
```

``` r

# compute difference in orientations
a <- 5
b <- 170
angle_diff_180(a, b)
#> [1] 15
```

``` r

# compute difference in 360° space (e.g., motion directions)
angle_diff_360(a, b)
#> [1] -165
```

``` r

# compute correlation between angles
data <- rmvn(10000, c(0, 0), V = matrix(c(1, 0.5, 0.5, 1), ncol = 2))
circ_corr(data[, 1], data[, 2])
#> [1] 0.4391056
```

The only (somewhat) complicated function is `remove_cardinal_biases`,
see the [help
files](https://achetverikov.github.io/circhelp/reference/remove_cardinal_biases.html)
and [the
vignette](https://achetverikov.github.io/circhelp/articles/cardinal_biases.html)
for an example use case.
<img src="https://achetverikov.github.io/circhelp/articles/cardinal_biases_files/figure-html/correct_biases_ex-1.png" title="Example of cardinal biases processing" width="100%"/>
