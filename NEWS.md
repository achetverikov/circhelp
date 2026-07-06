# circhelp 1.3.1

* Added `circ_mean()`, `circ_sd()`, and `angle_diff()`, generalized versions
  of the `_rad`/`_180`/`_360` family that accept an arbitrary `period`
  argument (e.g., 180, 360, or anything else).
* Added a `period` argument to `weighted_circ_mean()` to compute the weighted
  circular mean directly on data in an arbitrary circular space, without
  manually converting to radians first.

# circhelp 1.3

* `remove_cardinal_biases()`: replaced `gamlss` dependency with a native
  P-spline + REML/EM implementation. Results are numerically equivalent;
  runtime is approximately 2× faster.
* Removed `gamlss` from package imports.

# circhelp 1.2

* Added `density_asymmetry_discrete()`.
* Added `smoothed_circ_sd()`.
* Improved `density_asymmetry()` speed and options.

# circhelp 1.1.3

* Added `smoothed_circ_sd()`.
* Improved speed of `density_asymmetry()` by using matrix multiplication.
* Added `"average"` as a bandwidth option in `density_asymmetry()`.
* Updated imports (`ggplot2`, `stats`) and documentation.

# circhelp 1.1.2

* Added `density_asymmetry_discrete()`.
* Added `return_full_density` option to `density_asymmetry()`.
* Updated vignettes and package documentation/site build.

# circhelp 1.1.1

* Added `density_asymmetry()` to analyze the asymmetry in response probabilities.

# circhelp 1.1

* Initial CRAN submission.
