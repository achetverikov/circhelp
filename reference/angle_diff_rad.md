# Differences between angles in different circular spaces

Differences between angles in different circular spaces

## Usage

``` r
angle_diff_rad(a, b)

angle_diff_360(a, b)

angle_diff_180(a, b)

angle_diff_90(a, b)

angle_diff_180_45(a, b)

angle_diff_360_90(a, b)
```

## Arguments

- a:

  first angle

- b:

  second angle

## Value

difference between a and b

## Details

By default, all functions return values in ± half-range space (e.g., -pi
to pi for 2pi radian space used by `angle_diff_rad()`) but
`angle_diff_180_45()` and `angle_diff_360_90()` return values in \[-1/4
range, 3/4 range\] space

## Functions

- `angle_diff_rad()`: angle difference in radians

- `angle_diff_360()`: angle difference in 360 degree space

- `angle_diff_180()`: angle difference in 180 degree space (e.g., line
  orientation)

- `angle_diff_90()`: angle difference in 90 degree space

- `angle_diff_180_45()`: angle difference in 180 degree space from -45
  to 135

- `angle_diff_360_90()`: angle difference in 360 degree space from -90
  to 270

## Examples

``` r
angle_diff_180(5, 175)
#> [1] 10
angle_diff_360(5, 175)
#> [1] -170
angle_diff_90(5, 175)
#> [1] 10
angle_diff_rad(5, 175)
#> [1] -0.3539967

angle_diff_360(300, 0)
#> [1] -60
angle_diff_360_90(300, 0)
#> [1] -60
```
