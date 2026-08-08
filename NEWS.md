# circhelp 1.4.0

* `density_asymmetry()`: the kernel density estimate is now wrapped around the
  circle (`wrap = TRUE`), so mass near +/- `circ_space`/2 is no longer truncated.
  This is a no-op to ~1e-13 for orientation or colour errors, but it matters for
  data with a 180-degrees-off reversal mode, such as motion direction.
* `density_asymmetry()`: the +/- `circ_space`/2 point is now excluded from the
  positive and negative sums (`exclude_antipode = TRUE`), like zero. It is the
  same angle reached from two sides, and with an odd `n` it was counted twice.
* `density_asymmetry()`: a bandwidth narrower than the grid spacing now rescales
  `yvar` and the bandwidth together (`rescale_narrow = TRUE`) so the kernel spans
  at least one cell. Without it, the rectangle sum aliases: on a narrow error
  distribution with `bw.SJ()`-selected bandwidth the returned asymmetry was off by
  7.5% of its value at the default `n = 181`, now 0.4%.
* `density_asymmetry()`: the two signs are now paired on a rounded `abs(x)`, since
  `seq()` is not bitwise symmetric for every `n` (on `n = 1801`, 604 of 900 pairs
  failed to match and dropped out of the sums entirely).
* `density_asymmetry_discrete()`: the same four fixes, with the same argument
  names and defaults. Wrapping there is done by evaluating `stats::density()` over
  the neighbouring periods and folding the images back, since augmenting the sample
  with shifted copies would run into `density()`'s binning range.
* The fixes above default to the corrected behaviour; pass `wrap = FALSE`,
  `exclude_antipode = FALSE`, `rescale_narrow = FALSE` for pre-1.4.0 results.

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
