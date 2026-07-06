# Weighted circular parameters

Weighted circular parameters

## Usage

``` r
weighted_circ_mean(x, w, period = NULL, na.rm = FALSE)

weighted_circ_mean2(x, w, na.rm = FALSE)

weighted_circ_sd(x, w, na.rm = FALSE)

weighted_circ_rho(x, w, na.rm = FALSE)
```

## Arguments

- x:

  vector of values (in radians, unless `period` is set)

- w:

  vector of weights

- period:

  if not `NULL`, the period of the circular space (e.g., 180 or 360);
  `x` is then converted to radians as `x / period * 2 * pi` and the
  result is converted back to the same period

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
#> [1] -0.0007545474
weighted_circ_mean(x, w)
#> [1] 6.034346e-05
```
