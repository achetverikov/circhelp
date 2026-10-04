# Submission

This is a new release, version 1.4.0 (previous CRAN version: 1.1).

## Main changes since 1.1

* Added `density_asymmetry()` and `density_asymmetry_discrete()` for analyses
  of asymmetry in response probabilities.
* Added `smoothed_circ_sd()`.
* Added general `circ_mean()`, `circ_sd()`, and `angle_diff()` helpers for
  arbitrary circular periods.
* Replaced the `gamlss` dependency in `remove_cardinal_biases()` with a native
  P-spline + REML/EM implementation.
* Fixed circular wrapping, antipode handling, grid pairing, and narrow-bandwidth
  aliasing in the density-asymmetry estimators.
* Expanded tests and vignettes and reduced check-time test workloads.

## R CMD check results

Final CRAN-style check results are recorded by the R-CMD-check GitHub Actions
workflow for this release.
