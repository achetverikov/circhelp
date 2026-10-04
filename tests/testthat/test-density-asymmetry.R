make_dt <- function(n = 300, err_sd = 12, boundary_frac = 0, seed = 1) {
  set.seed(seed)
  fd <- runif(n, 0, 90)
  err <- 6 * sin(fd * 4 * pi / 180) + rnorm(n, 0, err_sd)
  if (boundary_frac > 0) {
    flip <- runif(n) < boundary_frac
    err[flip] <- 90 + rnorm(sum(flip), 0, 5)
  }
  data.table(abs_td_dist = fd, bias_to_distr_corr = (err + 90) %% 180 - 90)
}

signed_normal_mass <- function(mu, bw) {
  2 * pnorm(mu / bw) - 1
}

weighted_signed_mass <- function(dt, bw, weights_sd = 10) {
  dist_grid <- 1:90
  weights <- exp(-0.5 * outer(
    dt$abs_td_dist, dist_grid,
    FUN = function(x, d) ((x - d) / weights_sd)^2
  ))
  weights <- sweep(weights, 2, colSums(weights), FUN = "/")
  as.vector(crossprod(signed_normal_mass(dt$bias_to_distr_corr, bw), weights))
}

test_that("the asymmetry is a mass and narrow kernels do not alias", {
  dt <- make_dt(err_sd = 2)
  bw <- 0.05
  truth <- weighted_signed_mass(dt, bw)

  mass <- function(rescale) {
    r <- density_asymmetry(
      dt, kernel_bw = bw, n = 181,
      normalize = FALSE, rescale_narrow = rescale
    )
    setorder(r, dist)
    r$delta
  }

  scaled_error <- max(abs(mass(TRUE) - truth)) / diff(range(truth))
  unscaled_error <- max(abs(mass(FALSE) - truth)) / diff(range(truth))

  expect_lt(scaled_error, 0.03)
  expect_gt(unscaled_error, 0.08)
  expect_lt(scaled_error, unscaled_error / 3)
})

test_that("rescaling does not fire when the bandwidth already spans a cell", {
  dt <- make_dt()
  r <- density_asymmetry(dt, n = 181, normalize = FALSE)
  expect_equal(attr(r, "bias_scale"), 1)
  expect_equal(
    r$delta,
    density_asymmetry(dt, n = 181, normalize = FALSE, rescale_narrow = FALSE)$delta
  )
})

test_that("the two signs are paired on every grid, including non-bitwise-symmetric grids", {
  dt <- make_dt(n = 200)
  for (n in c(181, 361, 1801)) {
    r <- density_asymmetry(dt, kernel_bw = 5, n = n, normalize = FALSE)
    expect_equal(nrow(r), 90)
    expect_false(any(is.na(r$delta)))
  }
})

test_that("wrapping keeps mass that sits on the boundary", {
  dt <- make_dt(boundary_frac = 0.15)
  bw <- bw.SJ(dt$bias_to_distr_corr)
  wrapped <- density_asymmetry(dt, kernel_bw = bw, normalize = FALSE)
  truncated <- density_asymmetry(dt, kernel_bw = bw, normalize = FALSE, wrap = FALSE)
  expect_gt(max(abs(wrapped$delta - truncated$delta)), 1e-3)

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
  expect_equal(raw(), raw(exclude_antipode = FALSE), tolerance = 1e-12)

  norm <- function(...) density_asymmetry(dt, kernel_bw = bw, normalize = TRUE, ...)$delta
  expect_gt(max(abs(norm() - norm(exclude_antipode = FALSE))), 1e-4)
  expect_gt(
    max(abs(raw(wrap = FALSE) - raw(wrap = FALSE, exclude_antipode = FALSE))),
    1e-4
  )
})

test_that("the estimate is equivariant to the circular space it is expressed in", {
  dt <- make_dt()
  d180 <- density_asymmetry(dt, circ_space = 180, weights_sd = 10)
  d360 <- density_asymmetry(
    dt[, .(
      abs_td_dist = abs_td_dist * 2,
      bias_to_distr_corr = bias_to_distr_corr * 2
    )],
    circ_space = 360, weights_sd = 20
  )
  setorder(d180, dist)
  setorder(d360, dist)
  expect_equal(d180$delta, d360[dist %% 2 == 0]$delta, tolerance = 1e-12)
})

test_that("the wrapped discrete density conserves mass wherever the data sit", {
  dt <- make_dt(boundary_frac = 0.15)
  for (rot in c(0, 30, 70, 120)) {
    rotated <- dt[, .(
      bias_to_distr_corr = (bias_to_distr_corr + rot + 90) %% 180 - 90
    )]
    d <- density_asymmetry_discrete(
      rotated, kernel_bw = 5, return_full_density = TRUE
    )
    expect_equal(sum(d[abs(x) < 90]$y) + d[x == -90]$y, 1, tolerance = 1e-5)
  }

  rotated <- dt[, .(
    bias_to_distr_corr = (bias_to_distr_corr + 70 + 90) %% 180 - 90
  )]
  expect_lt(
    sum(density(
      rotated$bias_to_distr_corr,
      from = -90, to = 90, n = 181, bw = 5
    )$y),
    0.97
  )
})

test_that("the discrete asymmetry is stable for narrow kernels", {
  dt <- make_dt(err_sd = 2)
  bw <- 0.05
  truth <- mean(signed_normal_mass(dt$bias_to_distr_corr, bw))

  mass <- function(rescale) {
    density_asymmetry_discrete(
      dt, kernel_bw = bw, n = 181,
      normalize = FALSE, rescale_narrow = rescale
    )$delta
  }

  scaled_error <- abs(mass(TRUE) - truth)
  unscaled_error <- abs(mass(FALSE) - truth)
  expect_lt(scaled_error, 0.05)
  expect_lt(scaled_error, unscaled_error / 3)
})

test_that("discrete wrapping and antipode exclusion behave as expected", {
  dt <- make_dt(boundary_frac = 0.15)
  bw <- bw.SJ(dt$bias_to_distr_corr)
  raw <- function(...) {
    density_asymmetry_discrete(dt, kernel_bw = bw, normalize = FALSE, ...)$delta
  }

  expect_gt(abs(raw() - raw(wrap = FALSE)), 1e-3)
  expect_equal(raw(), raw(exclude_antipode = FALSE), tolerance = 1e-12)

  quiet <- make_dt()
  expect_lt(abs(
    density_asymmetry_discrete(
      quiet, kernel_bw = 5, normalize = FALSE
    )$delta -
      density_asymmetry_discrete(
        quiet, kernel_bw = 5, normalize = FALSE, wrap = FALSE
      )$delta
  ), 1e-4)
})
