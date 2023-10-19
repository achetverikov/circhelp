
<!-- README.md is generated from README.Rmd. Please edit that file -->

# circhelp

<!-- badges: start -->

`circhelp` is a small helper package for circular data analyses in R,
particularly useful for cognitive studies on orientation, motion
direction, and other circular features. The package contains functions
for descriptive statistics for circular data (computing means, SD, and
skewness), angular differences, and correlation. It also includes a
function to correct for cardinal biases in the human estimates of
circular features (e.g., orientation). <!-- badges: end -->

## Installation

You can install the current version from [GitHub](https://github.com/)
with:

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
library(mgcv)
#> Loading required package: nlme
#> This is mgcv 1.9-0. For overview type 'help("mgcv-package")'.
# compute a set of descriptive statistics
x <- rnorm(500)
circ_descr(x)
#> $mu
#> [1] -0.01653675
#> 
#> $sigma
#> [1] 1.004497
#> 
#> $skew_pewsey
#> [1] 0.001221028
#> 
#> $skew_fischer
#> [1] 0.02368438
#> 
#> $rho
#> [1] 0.6038032
#> 
#> $skew_rel_to_zero
#> [1] -0.003465756

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
#> [1] 0.4330472
```

The only (somewhat) complicated function is `remove_cardinal_biases`,
see the [help
files](https://achetverikov.github.io/circhelp/reference/remove_cardinal_biases.html)
and [the
vignette](https://achetverikov.github.io/circhelp/articles/cardinal_biases.html)
for an example use case.
<img src="https://achetverikov.github.io/circhelp/articles/cardinal_biases_files/figure-html/correct_biases_ex-1.png" title="Example of cardinal biases processing" width="100%"/>
