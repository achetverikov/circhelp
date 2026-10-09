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

Checked with `R CMD check --as-cran` on the built tarball:

* Local Windows 11, R 4.6.1: 0 errors | 0 warnings | 1 note (a `lastMiKTeXException`
  file left in the temp directory by the local MiKTeX installation, not by the package)
* win-builder, R 4.6.1: Status OK

The R-CMD-check GitHub Actions workflow additionally covers the other platforms.
This package has no reverse dependencies.
