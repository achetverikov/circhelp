test_that("constant-variance cardinal correction returns unavailable likelihood", {
  for (period in c(180, 360)) {
    x <- seq(-period / 2 + 0.5, period / 2 - 0.5, length.out = 400)
    err <- 2 * sin(x * pi / 45) + sin(seq_along(x))
    corrected <- suppressWarnings(remove_cardinal_biases(
      err, x, space = as.character(period), var_sigma = FALSE,
      init_outliers = abs(err) > period / 4
    ))
    expect_equal(nrow(corrected), length(x))
    expect_true(all(is.finite(corrected$be_c)))
    expect_true(all(is.na(corrected$total_log_lik)))
    expect_true(all(corrected$is_outlier %in% c(0, 1)))
  }
})

test_that("unsupported spline columns produce an actionable data error", {
  x <- seq(-45, 45, length.out = 220)
  basis <- circhelp:::.pb_basis(x)
  keep <- abs(x) < 12
  B <- basis$X[keep, , drop = FALSE]
  S <- crossprod(basis$D)
  y <- sin(x[keep] / 5) + 0.2 * cos(seq_len(sum(keep)))
  expect_lt(qr(B)$rank, ncol(B))
  error <- tryCatch(circhelp:::.lambda_em_pspline(y, B, S), error = identity)
  expect_s3_class(error, "circhelp_insufficient_data")
  expect_match(conditionMessage(error), "non-outlier observations")
  expect_match(conditionMessage(error), "design rank")
  expect_equal(error$n_non_outlier, nrow(B))
  expect_equal(error$n_coefficients, ncol(B))
  expect_lt(error$design_rank, error$n_coefficients)
})

test_that("full-rank solver retains the penalized least-squares solution", {
  x <- seq(-45, 45, length.out = 220)
  basis <- circhelp:::.pb_basis(x)
  y <- sin(x / 10) + 0.1 * cos(seq_along(x))
  S <- crossprod(basis$D)
  fit <- circhelp:::.lambda_em_pspline(y, basis$X, S)
  reference <- solve(crossprod(basis$X) + fit$lambda * S,
                     crossprod(basis$X, y))
  expect_equal(as.vector(basis$X %*% fit$beta),
               as.vector(basis$X %*% reference), tolerance = 1e-7)
})

test_that("location-scale fitting preserves deficient-support errors", {
  x <- seq(-45, 45, length.out = 220)
  basis <- circhelp:::.pb_basis(x)
  keep <- abs(x) < 12
  y <- 2 * sin(x / 10) + cos(seq_along(x))
  y[!keep] <- 70
  error <- tryCatch(circhelp:::fit_weighted_locscale_normal(
    y, x, as.numeric(keep), x_basis = basis$X,
    S_penalty = crossprod(basis$D)
  ), error = identity)
  expect_s3_class(error, "circhelp_insufficient_data")
  expect_match(conditionMessage(error), "design rank")
})

test_that("too few non-outliers cannot silently select a constant spline fit", {
  x <- seq(-45, 45, length.out = 220)
  basis <- circhelp:::.pb_basis(x)
  keep <- seq_along(x) <= 8
  error <- tryCatch(circhelp:::fit_weighted_locscale_normal(
    sin(x), x, as.numeric(keep), x_basis = basis$X,
    S_penalty = crossprod(basis$D)
  ), error = identity)
  expect_s3_class(error, "circhelp_insufficient_data")
  expect_equal(error$n_non_outlier, 8)
  expect_equal(error$min_observations, ncol(basis$X) + 2L)
})

test_that("public cardinal correction identifies the deficient bin", {
  x <- seq(-89.5, 89.5, length.out = 440)
  err <- 2 * sin(x * pi / 45) + sin(seq_along(x))
  error <- tryCatch(remove_cardinal_biases(
    err, x, bias_type = "card", reassign_at_boundaries = FALSE,
    init_outliers = abs(angle_diff_90(x, 0)) > 12
  ), error = identity)
  expect_s3_class(error, "circhelp_insufficient_data")
  expect_match(conditionMessage(error), "Final cardinal-bias fit in bin")
  expect_true(!is.null(error$bin))
})
