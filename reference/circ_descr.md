# A set of descriptive statistics for circular data

A set of descriptive statistics for circular data

## Usage

``` r
circ_descr(x, w = NULL, d = NULL, na.rm = FALSE)
```

## Arguments

- x:

  vector of angles

- w:

  weights for the values in the vector

- d:

  correction for the bias for data with known spacing

- na.rm:

  a logical value indicating whether NA values should be removed before
  the computation proceeds

## Value

a list with descriptive statistics

- mu - mean

- sigma - standard deviation

- skew_pewsey - skewness as defined by Pewsey

- skew_fischer - skewness as defined by Fischer

- rho - mean resultant length

- skew_rel_to_zero - skewness relative to zero

## Examples

``` r
x <- c(rnorm(50, 0, 0.5), rnorm(20, 1, 0.5))
circ_descr(x)
#> $mu
#> [1] 0.2206685
#> 
#> $sigma
#> [1] 0.6816897
#> 
#> $skew_pewsey
#> [1] -0.01884082
#> 
#> $skew_fischer
#> [1] -1.791742
#> 
#> $rho
#> [1] 0.7926683
#> 
#> $skew_rel_to_zero
#> [1] 0.1350797
#> 
```
