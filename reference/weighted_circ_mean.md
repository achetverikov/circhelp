# Weighted circular parameters

Weighted circular parameters

## Usage

``` r
weighted_circ_mean(x, w, na.rm = FALSE)

weighted_circ_mean2(x, w, na.rm = FALSE)

weighted_circ_sd(x, w, na.rm = FALSE)

weighted_circ_rho(x, w, na.rm = FALSE)
```

## Arguments

- x:

  vector of values (in radians)

- w:

  vector of weights

- na.rm:

  a logical value indicating whether NA values should be removed before
  the computation proceeds

## Value

weighted mean of values in the vector

## Functions

- `weighted_circ_mean()`: weighted circular mean

- `weighted_circ_mean2()`: an alternative way to compute weighted
  circular mean (the results are the same)

- `weighted_circ_sd()`: weighted circular SD

- `weighted_circ_rho()`: weighted mean resultant length

## Examples

``` r
x <- rnorm(1000, 0, 0.5)
w <- runif(1000, 0, 1)
weighted.mean(x, w)
#> [1] -0.0002040864
weighted_circ_mean(x, w)
#> [1] 0.001082165
```
