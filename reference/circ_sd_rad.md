# Circular standard deviation

Circular standard deviation

## Usage

``` r
circ_sd_rad(x, na.rm = FALSE)

circ_sd_360(x, na.rm = FALSE)

circ_sd_180(x, na.rm = FALSE)
```

## Arguments

- x:

  vector of angles

- na.rm:

  a logical value indicating whether NA values should be removed before
  the computation proceeds

## Value

standard deviation of values in the vector

## Functions

- `circ_sd_rad()`: SD of angles in radians

- `circ_sd_360()`: SD of angles in 360 degree space

- `circ_sd_180()`: SD of angles in 180 degree space

## Examples

``` r
circ_sd_rad(rnorm(50))
#> [1] 1.004706
circ_sd_180(rnorm(50))
#> [1] 0.8488329
```
