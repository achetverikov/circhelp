# Circular mean

Circular mean

## Usage

``` r
circ_mean_rad(x, na.rm = FALSE)

circ_mean_180(x, na.rm = FALSE)

circ_mean_360(x, na.rm = FALSE)

circ_mean(x, period = 360, na.rm = FALSE)
```

## Arguments

- x:

  vector of values

- na.rm:

  a logical value indicating whether NA values should be removed before
  the computation proceeds

- period:

  the period of the circular space (e.g., 180 for line orientations, 360
  for compass directions)

## Value

mean of values in the vector

## Functions

- `circ_mean_rad()`: circular mean in 2pi space

- `circ_mean_180()`: circular mean in 180° space (e.g., line
  orientation)

- `circ_mean_360()`: circular mean in 360° space

- `circ_mean()`: circular mean in a space with an arbitrary period
  (e.g., 180 or 360)

## Examples

``` r
x <- runif(1000, -pi, pi)
mean(x)
#> [1] -0.04127712
circ_mean_rad(x)
#> [1] -0.9695238
```
