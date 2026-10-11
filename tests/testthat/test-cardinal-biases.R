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

test_that("spline penalty identifies unsupported columns after exclusions", {
  x <- seq(-45, 45, length.out = 220)
  basis <- circhelp:::.pb_basis(x)
  keep <- abs(x) < 12
  B <- basis$X[keep, , drop = FALSE]
  S <- crossprod(basis$D)
  y <- sin(x[keep] / 5) + 0.2 * cos(seq_len(sum(keep)))
  expect_lt(qr(B)$rank, ncol(B))

  fit <- circhelp:::.lambda_em_pspline(y, B, S)
  reference <- solve(crossprod(B) + fit$lambda * S, crossprod(B, y))
  expect_true(all(is.finite(fit$beta)))
  expect_equal(as.vector(B %*% fit$beta), as.vector(B %*% reference),
               tolerance = 1e-7)
  reference_edf <- sum(diag(solve(crossprod(B) + fit$lambda * S, crossprod(B))))
  expect_equal(fit$edf, reference_edf, tolerance = 1e-7)
})

test_that("rank-deficient and full-rank spline decompositions agree", {
  x <- seq(-45, 45, length.out = 220)
  basis <- circhelp:::.pb_basis(x)
  y <- sin(x / 10) + 0.1 * cos(seq_along(x))
  S <- crossprod(basis$D)
  fast <- circhelp:::.lambda_em_pspline(y, basis$X, S)
  augmented <- circhelp:::.lambda_em_pspline_rank_deficient(
    y, basis$X, S, lambda = 10, max_iter = 50L, tol = 1e-7,
    penalty_order = 2L
  )
  expect_equal(augmented$lambda, fast$lambda, tolerance = 1e-6)
  expect_equal(augmented$edf, fast$edf, tolerance = 1e-6)
  expect_equal(as.vector(basis$X %*% augmented$beta),
               as.vector(basis$X %*% fast$beta), tolerance = 1e-6)
})

test_that("location-scale fitting survives excluded spline support", {
  x <- seq(-45, 45, length.out = 220)
  basis <- circhelp:::.pb_basis(x)
  keep <- abs(x) < 12
  y <- 2 * sin(x / 10) + cos(seq_along(x))
  y[!keep] <- 70
  fit <- circhelp:::fit_weighted_locscale_normal(
    y, x, as.numeric(keep), x_basis = basis$X,
    S_penalty = crossprod(basis$D)
  )
  expect_true(fit$ok)
  expect_true(all(is.finite(fit$pred)))
  expect_true(all(is.finite(fit$pred_sigma) & fit$pred_sigma > 0))
})
