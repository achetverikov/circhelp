# Pad circular data on both ends

Pad circular data on both ends

## Usage

``` r
pad_circ(
  data,
  circ_var,
  circ_borders = c(-90, 90),
  circ_part = 1/6,
  verbose = FALSE
)
```

## Arguments

- data:

  data.table to pad

- circ_var:

  circular variable

- circ_borders:

  range of the circular variable

- circ_part:

  padding proportion

- verbose:

  print extra info

## Value

a padded data.table

## Details

Pads the data by adding a part of the data (default: 1/6th) from one end
to another end. Useful to roughly account for circularity when using
non-circular methods.

## Examples

``` r

dt <- data.table::data.table(x = runif(1000, -90, 90), y = rnorm(1000))
pad_circ(dt, "x", verbose = TRUE)
#> Rows in original DT: 1000, padded on the left: 190, padded on the right: 170
#>                 x          y
#>             <num>      <num>
#>    1:  -46.870792 -1.5570357
#>    2:   26.597677  1.9231637
#>    3:   85.620737 -1.8568296
#>    4:  -21.960215 -2.1061184
#>    5:   -6.454057  0.6976485
#>   ---                       
#> 1356:  -93.428947 -0.9763090
#> 1357:  -96.305980  1.1929717
#> 1358: -113.204793 -0.3776643
#> 1359: -109.728215  0.7502717
#> 1360: -110.079715 -1.7150032
```
