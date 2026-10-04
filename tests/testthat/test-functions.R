test_that("circular means agree across representations", {
  set.seed(1)
  x <- rnorm(500, sd = 2)

  expect_equal(circ_mean_rad(x), circ_descr(x)[["mu"]])
  expect_equal(circ_mean_rad(x), circ_mean_360(x / pi * 180) / 180 * pi)
  expect_equal(circ_mean_rad(x), weighted_circ_mean(x, rep(1, length(x))))
  expect_equal(circ_mean_rad(x), weighted_circ_mean2(x, rep(1, length(x))))
  expect_equal(circ_mean_rad(c(-pi / 2, 0, pi / 2)), 0, tolerance = 1e-15)
})

test_that("weighted means are equal", {
  set.seed(2)
  x <- rnorm(500, sd = 2)
  w <- runif(500)

  expect_equal(weighted_circ_mean(x, w), weighted_circ_mean2(x, w))
})

test_that("circular SDs agree across representations", {
  set.seed(3)
  x <- rnorm(500, sd = 2)

  expect_equal(circ_sd_rad(x), circ_descr(x)[["sigma"]])
  expect_equal(circ_sd_rad(x), circ_sd_360(x / pi * 180) / 180 * pi)
  expect_equal(circ_sd_rad(x), weighted_circ_sd(x, rep(1, length(x))))
})

test_that("Circular correlation has the expected limiting cases", {
  set.seed(4)
  x <- runif(5000, -pi, pi)

  expect_equal(circ_corr(x, x), 1, tolerance = 1e-12)
  expect_equal(circ_corr(x, x + 0.25), 1, tolerance = 1e-12)
})

test_that("Circular correlation is close to Pearson correlation for narrow data", {
  set.seed(5)
  n <- 5000
  x <- rnorm(n, sd = 0.1)
  y <- 0.5 * x + sqrt(0.75) * rnorm(n, sd = 0.1)

  expect_equal(circ_corr(x, y), cor(x, y), tolerance = 0.01)
})

test_that("conversion from circular SD to kappa works both ways", {
  test_sd_deg <- c(5, 25, 60)
  test_sd_rad <- test_sd_deg / 180 * pi

  kappa_from_deg <- vm_circ_sd_deg_to_kappa(test_sd_deg)
  kappa_from_rad <- vm_circ_sd_to_kappa(test_sd_rad)
  expect_equal(kappa_from_deg, kappa_from_rad, tolerance = 1e-10)
  expect_equal(vm_kappa_to_circ_sd_deg(kappa_from_deg), test_sd_deg, tolerance = 1e-3)
  expect_equal(vm_kappa_to_circ_sd(kappa_from_deg), test_sd_rad, tolerance = 1e-3)
})

test_that("weighted_sem implements the documented Kirchner formula", {
  x <- c(-2, 0.5, 1, 4, 7)
  w <- c(1, 2, 3, 2, 1)
  wn <- w / sum(w)
  xbar <- weighted.mean(x, wn)
  var_w <- (sum(wn * x^2) - xbar^2) / (1 - sum(wn^2))
  expected <- sqrt(var_w * sum(wn^2))

  expect_equal(weighted_sem(x, w), expected)
  expect_equal(weighted_sem(x, rep(1, length(x))), sd(x) / sqrt(length(x)))
})
