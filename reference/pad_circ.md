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
#> Rows in original DT: 1000, padded on the left: 172, padded on the right: 180
#>                x           y
#>            <num>       <num>
#>    1:   87.28217  0.36960518
#>    2:   12.93592  0.43279388
#>    3:  -22.05991  0.05649643
#>    4:   33.47731  0.21501491
#>    5:   17.92069  1.39877788
#>   ---                       
#> 1348:  -95.52749 -0.37747204
#> 1349:  -91.50992  0.18173318
#> 1350: -112.62043  0.59092248
#> 1351:  -95.47732  0.19596330
#> 1352: -111.12037  1.37574918
```
