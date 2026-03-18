# Circular correlation coefficient

Computes a circular correlation coefficient as defined in Jammalamadaka
& SenGupta (2001).

## Usage

``` r
circ_corr(a, b, ill_defined = FALSE, mu = NULL, na.rm = FALSE)
```

## Arguments

- a:

  first variable

- b:

  second variable

- ill_defined:

  is one of the variables mean is not well-defined (e.g., it is
  uniformly distributed)?

- mu:

  fix the mean parameter of both vectors to a certain value

- na.rm:

  a logical value indicating whether NA values should be removed before
  the computation proceeds

## Value

correlation coefficient

## References

Jammalamadaka, S. R., & SenGupta, A. (2001). Topics in Circular
Statistics. WORLD SCIENTIFIC.
[doi:10.1142/4031](https://doi.org/10.1142/4031)

## Examples

``` r
requireNamespace("mgcv")
data <- mgcv::rmvn(10000, c(0, 0), V = matrix(c(1, 0.5, 0.5, 1), ncol = 2))
circ_corr(data[, 1], data[, 2])
#> [1] 0.4472373
```
