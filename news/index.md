# Changelog

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
