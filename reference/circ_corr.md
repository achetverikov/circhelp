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
set.seed(1)
x <- rnorm(1000)
y <- 0.5 * x + sqrt(0.75) * rnorm(1000)
circ_corr(x, y)
#> [1] 0.4095552
```
