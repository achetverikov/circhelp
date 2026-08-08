make_dt <- function(n = 600, err_sd = 12, boundary_frac = 0, seed = 1) {
  set.seed(seed)
  fd <- runif(n, 0, 90)
  err <- 6 * sin(fd * 4 * pi / 180) + rnorm(n, 0, err_sd)
  if (boundary_frac > 0) {
    flip <- runif(n) < boundary_frac
    err[flip] <- 90 + rnorm(sum(flip), 0, 5)
  }
  data.table(abs_td_dist = fd, bias_to_distr_corr = (err + 90) %% 180 - 90)
}

test_that("the asymmetry is a mass and stops depending on the grid", {
  # Narrow errors drive bw.SJ() below the grid spacing, where the rectangle sum
  # aliases unless the data are rescaled to fit the grid.
  dt <- make_dt(err_sd = 2)
  bw <- 0.05 # a twentieth of the default grid spacing; bw.SJ() reaches this on real data
  mass <- function(n, rescale) {
    r <- density_asymmetry(dt, kernel_bw = bw, n = n, normalize = FALSE, rescale_narrow = rescale)
    setorder(r, dist)$delta * (180 / (n - 1))
  }
  converged <- mass(180001, FALSE) # dx << bw: rescaling is inactive either way
  expect_equal(mass(180001, TRUE), converged)
  err <- function(rescale) max(abs(mass(181, rescale) - converged)) / diff(range(converged))
  # measured: 1.1% of range with rescaling, 12% without (the safety cap binds here,
  # so the bandwidth floor still does part of the work)
  expect_lt(err(TRUE), 0.02)
  expect_gt(err(FALSE), 0.1)
  expect_lt(err(TRUE), err(FALSE) / 5)

  # and at a bandwidth bw.SJ() actually picks, the rescaled estimate is still the
  # closer one by a wide margin
  dt2 <- make_dt(err_sd = 0.1)
  sj <- bw.SJ(dt2$bias_to_distr_corr)
  mass2 <- function(n, rescale) {
    r <- density_asymmetry(dt2, kernel_bw = sj, n = n, normalize = FALSE, rescale_narrow = rescale)
    setorder(r, dist)$delta * (180 / (n - 1))
  }
  converged2 <- mass2(18001, FALSE)
  err2 <- function(rescale) max(abs(mass2(181, rescale) - converged2)) / diff(range(converged2))
  expect_lt(err2(TRUE), err2(FALSE) / 3)
})

test_that("rescaling does not fire when the bandwidth already spans a cell", {
  dt <- make_dt()
  r <- density_asymmetry(dt, n = 181, normalize = FALSE)
  expect_equal(attr(r, "bias_scale"), 1)
  expect_equal(r$delta, density_asymmetry(dt, n = 181, normalize = FALSE, rescale_narrow = FALSE)$delta)
})

test_that("the two signs are paired on every grid, not just the lucky ones", {
  # seq() is not bitwise symmetric for every n; unpaired points used to drop out.
  dt <- make_dt()
  for (n in c(181, 361, 1801)) {
    r <- density_asymmetry(dt, kernel_bw = 5, n = n, normalize = FALSE)
    expect_equal(nrow(r), 90)
    expect_false(any(is.na(r$delta)))
    # sum of densities scales with 1/dx, so the mass must not
    expect_equal(r$delta[1] * (180 / (n - 1)),
      density_asymmetry(dt, kernel_bw = 5, n = 18001, normalize = FALSE)$delta[1] * (180 / 18000),
      tolerance = 1e-3
    )
  }
})

test_that("wrapping keeps mass that sits on the boundary", {
  dt <- make_dt(boundary_frac = 0.15)
  bw <- bw.SJ(dt$bias_to_distr_corr)
  wrapped <- density_asymmetry(dt, kernel_bw = bw, normalize = FALSE)
  truncated <- density_asymmetry(dt, kernel_bw = bw, normalize = FALSE, wrap = FALSE)
  expect_gt(max(abs(wrapped$delta - truncated$delta)), 1e-3)
  # errors away from the boundary: wrapping changes nothing that matters
  quiet <- make_dt()
  expect_equal(
    density_asymmetry(quiet, kernel_bw = 5, normalize = FALSE)$delta,
    density_asymmetry(quiet, kernel_bw = 5, normalize = FALSE, wrap = FALSE)$delta,
    tolerance = 1e-10
  )
})

test_that("the antipode is signless under wrapping but still inflates the total", {
  dt <- make_dt(boundary_frac = 0.15)
  bw <- bw.SJ(dt$bias_to_distr_corr)
  raw <- function(...) density_asymmetry(dt, kernel_bw = bw, normalize = FALSE, ...)$delta
  # its two images cancel in the difference ...
  expect_equal(raw(), raw(exclude_antipode = FALSE), tolerance = 1e-12)
  # ... but not in the denominator normalize = TRUE divides by
  norm <- function(...) density_asymmetry(dt, kernel_bw = bw, normalize = TRUE, ...)$delta
  expect_gt(max(abs(norm() - norm(exclude_antipode = FALSE))), 1e-4)
  # with a truncated kernel it does not even cancel in the difference
  expect_gt(max(abs(raw(wrap = FALSE) - raw(wrap = FALSE, exclude_antipode = FALSE))), 1e-4)
})

test_that("the estimate is equivariant to the circular space it is expressed in", {
  dt <- make_dt()
  d180 <- density_asymmetry(dt, circ_space = 180, weights_sd = 10)
  d360 <- density_asymmetry(dt[, .(abs_td_dist = abs_td_dist * 2, bias_to_distr_corr = bias_to_distr_corr * 2)],
    circ_space = 360, weights_sd = 20
  )
  setorder(d180, dist)
  setorder(d360, dist)
  expect_equal(d180$delta, d360[dist %% 2 == 0]$delta, tolerance = 1e-12)
})

test_that("the wrapped discrete density conserves mass wherever the data sit", {
  dt <- make_dt(boundary_frac = 0.15)
  for (rot in c(0, 30, 70, 120)) {
    rotated <- dt[, .(bias_to_distr_corr = (bias_to_distr_corr + rot + 90) %% 180 - 90)]
    d <- density_asymmetry_discrete(rotated, kernel_bw = 5, return_full_density = TRUE)
    # sum over the half-open circle: the +90 cell is the -90 cell, and counting both
    # is exactly the error the naive 181-point sum makes
    expect_equal(sum(d[abs(x) < 90]$y) + d[x == -90]$y, 1, tolerance = 1e-5)
  }
  # truncating instead loses mass once the data reach the boundary
  rotated <- dt[, .(bias_to_distr_corr = (bias_to_distr_corr + 70 + 90) %% 180 - 90)]
  expect_lt(sum(density(rotated$bias_to_distr_corr, from = -90, to = 90, n = 181, bw = 5)$y), 0.97)
})

test_that("the discrete asymmetry is a mass and stops depending on the grid", {
  dt <- make_dt(err_sd = 2)
  mass <- function(n, rescale) {
    density_asymmetry_discrete(dt,
      kernel_bw = 0.05, n = n, normalize = FALSE, rescale_narrow = rescale
    )$delta * (180 / (n - 1))
  }
  converged <- mass(18001, FALSE)
  # measured: 2.9% off with rescaling, 131% off without
  expect_equal(mass(181, TRUE), converged, tolerance = 0.05)
  expect_gt(abs(mass(181, FALSE) - converged), abs(converged))
})

test_that("discrete: wrapping and antipode exclusion behave as in density_asymmetry", {
  dt <- make_dt(boundary_frac = 0.15)
  bw <- bw.SJ(dt$bias_to_distr_corr)
  raw <- function(...) density_asymmetry_discrete(dt, kernel_bw = bw, normalize = FALSE, ...)$delta
  expect_gt(abs(raw() - raw(wrap = FALSE)), 1e-3)
  expect_equal(raw(), raw(exclude_antipode = FALSE), tolerance = 1e-12)
  # errors away from the boundary: the wrapped images contribute nothing. The residual
  # is not zero because the two calls evaluate `density()` over different ranges and so
  # bin the data slightly differently; the wrapped mass itself is ~1e-51 here.
  quiet <- make_dt()
  expect_lt(abs(
    density_asymmetry_discrete(quiet, kernel_bw = 5, normalize = FALSE)$delta -
      density_asymmetry_discrete(quiet, kernel_bw = 5, normalize = FALSE, wrap = FALSE)$delta
  ), 1e-4)
})
