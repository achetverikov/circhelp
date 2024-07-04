test_that("means are equal", {
  x <- rnorm(500, sd = 2) # angle in radians
  expect_equal(circ_mean_rad(x), circ_descr(x)[["mu"]])
  expect_equal(circ_mean_rad(x), as.vector(circular::mean.circular(circular::circular(x))))
  expect_equal(circ_mean_rad(x), circ_mean_360(x / pi * 180) / 180 * pi)
  expect_equal(circ_mean_rad(x), weighted_circ_mean(x, rep(1, length(x))))
  expect_equal(circ_mean_rad(x), weighted_circ_mean2(x, rep(1, length(x))))
})


test_that("weighted means are equal", {
  x <- rnorm(500, sd = 2) # angle in radians
  w <- runif(500) # weights

  expect_equal(weighted_circ_mean(x, w), weighted_circ_mean2(x, w))
})

test_that("SDs are equal", {
  x <- rnorm(500, sd = 2) # angle in radians

  expect_equal(circ_sd_rad(x), circ_descr(x)[["sigma"]])
  expect_equal(circ_sd_rad(x), circ_sd_360(x / pi * 180) / 180 * pi)
  expect_equal(circ_sd_rad(x), weighted_circ_sd(x, rep(1, length(x))))
  expect_equal(circ_sd_rad(x), as.vector(circular::sd.circular(circular::circular(x))))
})

test_that("Circular correlation works properly for uniform and conditional von Mises", {
  # see Jammalamadaka & SenGupta, pp. 181-182
  kappa <- runif(1, 5, 200) # NB: for lower kappa, the results are not always within tolerance limits.
  n <- 1000000
  x <- runif(n, -pi, pi)
  y <- x + as.vector(circular::rvonmises(n, mu = circular::circular(0), kappa = kappa))
  exp_r <- a_fun(kappa)
  obs_r <- circ_corr(x, y, ill_defined = T)
  tolerance <- 1e-3
  print(sprintf("Difference between expected and observed correlation is %.6f", exp_r - obs_r))
  expect_lt(abs(exp_r - obs_r), tolerance)
})

test_that("Circular correlation close to Pearson for narrow cases", {
  # see Jammalamadaka & SenGupta
  data <- mgcv::rmvn(10000, c(0, 0), V = matrix(c(1, 0.5, 0.5, 1), ncol = 2))
  data <- data / 8
  pearson_r <- cor(data[, 1], data[, 2])
  obs_r <- circ_corr(data[, 1], data[, 2])
  tolerance <- 1e-3
  print(sprintf("Difference between expected and observed correlation is %.6f", pearson_r - obs_r))
  expect_lt(abs(pearson_r - obs_r), tolerance)
})

test_that("Circular correlation matches BAMBI::circ_cor", {
  data <- mgcv::rmvn(10000, c(0, 0), V = matrix(c(1, 0.5, 0.5, 1), ncol = 2))
  data <- data * 2
  obs_r <- circ_corr(data[, 1], data[, 2])
  exp_r <- BAMBI::circ_cor(data)[[1]]
  tolerance <- 1e-3
  print(sprintf("Difference between expected and observed correlation is %.6f", exp_r - obs_r))
  expect_lt(abs(exp_r - obs_r), tolerance)
})

test_that("conversion from circular SD to kappa works both ways", {
  test_sd_deg <- runif(1, 0, 100) # SD in degrees
  test_sd_rad <- test_sd_deg / 180 * pi

  kappa_from_deg <- vm_circ_sd_deg_to_kappa(test_sd_deg)
  kappa_from_rad <- vm_circ_sd_to_kappa(test_sd_rad)
  expect_equal(kappa_from_deg, kappa_from_rad)

  tolerance <- 1e-3
  expect_lt(abs(test_sd_deg - vm_kappa_to_circ_sd_deg(kappa_from_deg)), tolerance)
  expect_lt(abs(test_sd_rad - vm_kappa_to_circ_sd(kappa_from_deg)), tolerance)
})

test_that("Weighted sample of the mean is computed correctly", {
  n_obs <- 2000
  w <- runif(n_obs)
  #w <- rep(1, n_obs)
  w <- w/sum(w)
  x <- rnorm(n_obs, sd = 5)

  sample_means <- replicate(100000, {
    x_boot <- sample(x, n_obs, replace = T)
    mean_w <- stats::weighted.mean(x_boot, w)
    mean_w
  })

  bootstrap_sem <- sd(sample_means)
  deviation = abs((bootstrap_sem-weighted_sem(x, w))/bootstrap_sem)
  expect_lt(deviation, .05)
})
