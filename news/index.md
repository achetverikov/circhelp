# Changelog

## circhelp 1.3.1

- Added
  [`circ_mean()`](https://achetverikov.github.io/circhelp/index.html/reference/circ_mean_rad.md),
  [`circ_sd()`](https://achetverikov.github.io/circhelp/index.html/reference/circ_sd_rad.md),
  and
  [`angle_diff()`](https://achetverikov.github.io/circhelp/index.html/reference/angle_diff_rad.md),
  generalized versions of the `_rad`/`_180`/`_360` family that accept an
  arbitrary `period` argument (e.g., 180, 360, or anything else).
- Added a `period` argument to
  [`weighted_circ_mean()`](https://achetverikov.github.io/circhelp/index.html/reference/weighted_circ_mean.md)
  to compute the weighted circular mean directly on data in an arbitrary
  circular space, without manually converting to radians first.

## circhelp 1.3

- [`remove_cardinal_biases()`](https://achetverikov.github.io/circhelp/index.html/reference/remove_cardinal_biases.md):
  replaced `gamlss` dependency with a native P-spline + REML/EM
  implementation. Results are numerically equivalent; runtime is
  approximately 2× faster.
- Removed `gamlss` from package imports.

## circhelp 1.2

- Added
  [`density_asymmetry_discrete()`](https://achetverikov.github.io/circhelp/index.html/reference/density_asymmetry_discrete.md).
- Added
  [`smoothed_circ_sd()`](https://achetverikov.github.io/circhelp/index.html/reference/smoothed_circ_sd.md).
- Improved
  [`density_asymmetry()`](https://achetverikov.github.io/circhelp/index.html/reference/density_asymmetry.md)
  speed and options.

## circhelp 1.1.3

- Added
  [`smoothed_circ_sd()`](https://achetverikov.github.io/circhelp/index.html/reference/smoothed_circ_sd.md).
- Improved speed of
  [`density_asymmetry()`](https://achetverikov.github.io/circhelp/index.html/reference/density_asymmetry.md)
  by using matrix multiplication.
- Added `"average"` as a bandwidth option in
  [`density_asymmetry()`](https://achetverikov.github.io/circhelp/index.html/reference/density_asymmetry.md).
- Updated imports (`ggplot2`, `stats`) and documentation.

## circhelp 1.1.2

- Added
  [`density_asymmetry_discrete()`](https://achetverikov.github.io/circhelp/index.html/reference/density_asymmetry_discrete.md).
- Added `return_full_density` option to
  [`density_asymmetry()`](https://achetverikov.github.io/circhelp/index.html/reference/density_asymmetry.md).
- Updated vignettes and package documentation/site build.

## circhelp 1.1.1

- Added
  [`density_asymmetry()`](https://achetverikov.github.io/circhelp/index.html/reference/density_asymmetry.md)
  to analyze the asymmetry in response probabilities.

## circhelp 1.1

CRAN release: 2024-07-04

- Initial CRAN submission.
