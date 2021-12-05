test_that("means are equal", {
  x <- rnorm(500, sd = 2) # angle in radians
  expect_equal(circ_mean_rad(x), circ_descr(x)[['mu']])
  expect_equal(circ_mean_rad(x), as.vector(circular::mean.circular(circular::circular(x))))
  expect_equal(circ_mean_rad(x), circ_mean_360(x/pi*180)/180*pi)
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

  expect_equal(circ_sd_rad(x), circ_descr(x)[['sigma']])
  expect_equal(circ_sd_rad(x), circ_sd_360(x/pi*180)/180*pi)
  expect_equal(circ_sd_rad(x), as.vector(circular::sd.circular(circular::circular(x))))


})
