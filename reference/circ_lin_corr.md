# Circular-linear correlation

Implementation of the circular-linear correlation measure introduced by
Mardia (1976) and Johnson and Wehrly (1977) as cited in Jammalamadaka &
Sengupta (2001).

## Usage

``` r
circ_lin_corr(circ_x, lin_x, na.rm = FALSE)
```

## Arguments

- circ_x:

  circular variable

- lin_x:

  linear variable

- na.rm:

  a logical value indicating whether NA values should be removed before
  the computation proceeds

## Value

circular-linear correlation measure

## Details

This measure is computed as \\r^2 = (r\_{xc}^2+r\_{xs}^2-2 r\_{xc}
r\_{xs}r\_{cs})/(1-r\_{cs}^2)\\ where \\r\_{xc} = corr(x,
cos(\alpha))\\, \\r\_{xs} = corr(x, sin(\alpha))\\, \\r\_{cs} =
corr(cos(\alpha), sin(\alpha))\\, and \\\alpha\\ and \\x\\ are the
circular and linear variables, respectively.

## References

Jammalamadaka, S. R., & SenGupta, A. (2001). Topics in Circular
Statistics. WORLD SCIENTIFIC.
[doi:10.1142/4031](https://doi.org/10.1142/4031)

## Examples

``` r

set.seed(1)
x <- rnorm(50)
a <- rnorm(50, sd = 0.4)
circ_lin_corr(x + a, x)
#> [1] 0.864586
```
