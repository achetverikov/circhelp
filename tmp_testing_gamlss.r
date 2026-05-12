# test_pb_equivalence.R

library(gamlss)
library(gamlss.dist)

# -------------------------------------------------------------------
# 1. Plain-R helpers for default gamlss::pb()
# -------------------------------------------------------------------

.pb_basis_default <- function(x) {
  x <- as.vector(x)

  n <- length(x)
  n_unique <- length(unique(x))

  inter <- if (n < 99) 10 else 20
  inter <- min(inter, n_unique)

  degree <- 3
  order <- 2

  xl <- min(x)
  xr <- max(x)

  xmin <- xl - 0.01 * (xr - xl)
  xmax <- xr + 0.01 * (xr - xl)

  dx <- (xmax - xmin) / inter

  knots <- seq(
    from = xmin - degree * dx,
    to = xmax + degree * dx,
    by = dx
  )

  tpower <- function(x, t, p) {
    (x - t)^p * (x > t)
  }

  P <- outer(x, knots, tpower, degree)

  D_full <- diff(diag(ncol(P)), diff = degree + 1) /
    (gamma(degree + 1) * dx^degree)

  X <- (-1)^(degree + 1) * P %*% t(D_full)

  D_penalty <- diff(diag(ncol(X)), diff = order)

  list(
    X = X,
    D = D_penalty,
    inter = inter,
    degree = degree,
    order = order,
    knots = knots,
    xmin = xmin,
    xmax = xmax,
    dx = dx
  )
}

.pb_fit_default <- function(y, x, w = rep(1, length(y)),
                            lambda = 10,
                            max_iter = 50,
                            tol = 1e-7) {
  y <- as.vector(y)
  x <- as.vector(x)
  w <- as.vector(w)

  stopifnot(length(y) == length(x))
  stopifnot(length(y) == length(w))

  basis <- .pb_basis_default(x)

  X <- basis$X
  D <- basis$D

  n_eff <- sum(w != 0)
  p <- ncol(X)

  # gamlss::pb() default penalty order
  penalty_order <- 2

  qrX <- qr(sqrt(w) * X, tol = .Machine$double.eps^.8)
  R <- qr.R(qrX)
  Q <- qr.Q(qrX)
  Qy <- t(Q) %*% (sqrt(w) * y)

  regpen <- function(lambda) {
    RD <- rbind(R, sqrt(lambda) * D)
    s <- svd(RD)

    rank <- sum(s$d > max(s$d) * .Machine$double.eps^.8)

    U1 <- s$u[seq_len(p), seq_len(rank), drop = FALSE]
    y1 <- t(U1) %*% Qy

    beta <- s$v[, seq_len(rank), drop = FALSE] %*%
      (y1 / s$d[seq_len(rank)])

    H <- U1 %*% t(U1)
    edf <- sum(diag(H))

    list(
      beta = as.vector(beta),
      edf = edf
    )
  }

  sig2 <- NA_real_
  tau2 <- NA_real_

  for (iter in seq_len(max_iter)) {
    fit <- regpen(lambda)

    fitted <- as.vector(X %*% fit$beta)
    gamma <- D %*% fit$beta

    sig2 <- sum(w * (y - fitted)^2) / (n_eff - fit$edf)

    # Critical fix:
    # gamlss.pb() uses fit$edf - order, not fit$edf - nrow(D)
    tau2 <- sum(gamma^2) / (fit$edf - penalty_order)

    if (tau2 < 1e-7) {
      tau2 <- 1e-7
    }

    lambda_old <- lambda
    lambda <- sig2 / tau2

    if (lambda < 1e-7) {
      lambda <- 1e-7
    }

    if (lambda > 1e7) {
      lambda <- 1e7
    }

    if (abs(lambda - lambda_old) < tol || lambda > 1e10) {
      break
    }
  }

  fit <- regpen(lambda)
  fitted <- as.vector(X %*% fit$beta)

  list(
    fitted.values = fitted,
    residuals = y - fitted,
    beta = fit$beta,
    edf = fit$edf,
    nl.df = fit$edf - 2,
    lambda = lambda,
    sig2 = sig2,
    tau2 = tau2,
    X = X,
    D = D,
    knots = basis$knots,
    inter = basis$inter,
    iterations = iter
  )
}
# -------------------------------------------------------------------
# 2. One equivalence test against gamlss:::gamlss.pb()
# -------------------------------------------------------------------
compare_pb_once <- function(x, seed = 1, verbose = TRUE) {
  set.seed(seed)

  n <- length(x)

  y <- 0.1 * x +
    2 * sin(x / 8) +
    rnorm(n, sd = 0.5)

  w <- runif(n, 0.5, 2)

  fit_ours <- .pb_fit_default(
    y = y,
    x = x,
    w = w
  )

  fit_gamlss <- gamlss:::gamlss.pb(
    x = gamlss::pb(x),
    y = y,
    w = w
  )

  # Inspect once if needed:
  # print(names(fit_gamlss))
  # str(fit_gamlss)

  fitted_gamlss <- as.vector(fit_gamlss$fitted.values)

  lambda_gamlss <- if (!is.null(fit_gamlss$lambda)) {
    fit_gamlss$lambda
  } else {
    NA_real_
  }

  # gamlss.pb usually reports nl.df, not edf.
  # For pb(), nl.df is edf minus the constant/linear part.
  edf_gamlss <- if (!is.null(fit_gamlss$edf)) {
    fit_gamlss$edf
  } else if (!is.null(fit_gamlss$nl.df)) {
    fit_gamlss$nl.df + 2
  } else {
    NA_real_
  }

  out <- data.frame(
    n = n,
    n_unique = length(unique(x)),
    inter_ours = fit_ours$inter,

    max_abs_fitted = max(abs(fit_ours$fitted.values - fitted_gamlss)),
    mean_abs_fitted = mean(abs(fit_ours$fitted.values - fitted_gamlss)),
    cor_fitted = cor(fit_ours$fitted.values, fitted_gamlss),

    lambda_ours = fit_ours$lambda,
    lambda_gamlss = lambda_gamlss,
    lambda_diff = fit_ours$lambda - lambda_gamlss,

    edf_ours = fit_ours$edf,
    edf_gamlss = edf_gamlss,
    edf_diff = fit_ours$edf - edf_gamlss
  )

  if (verbose) {
    print(out)
  }

  invisible(list(
    summary = out,
    x = x,
    y = y,
    w = w,
    ours = fit_ours,
    gamlss = fit_gamlss
  ))
}

# -------------------------------------------------------------------
# 3. Run several cases relevant for circhelp
# -------------------------------------------------------------------

run_pb_equivalence_suite <- function() {
  cases <- list(
    dense_even = seq(-45, 45, length.out = 300),

    repeated_grid = rep(seq(-45, 45, length.out = 31), each = 10),

    small_n = seq(-45, 45, length.out = 80),

    integer_sample = sample(seq(-45, 45, by = 1), 300, replace = TRUE),

    narrow_range = seq(-10, 10, length.out = 300),

    asymmetric = sort(runif(300, -20, 45))
  )

  res <- lapply(seq_along(cases), function(i) {
    cat("\n--- Case:", names(cases)[i], "---\n")

    tmp <- compare_pb_once(
      x = cases[[i]],
      seed = i,
      verbose = TRUE
    )

    cbind(
      case = names(cases)[i],
      tmp$summary
    )
  })

  do.call(rbind, res)
}


# -------------------------------------------------------------------
# 4. Optional comparison against full gamlss()
#    This is NOT expected to be identical.
#    It checks how much the full outer gamlss() machinery differs.
# -------------------------------------------------------------------

compare_against_full_gamlss <- function(x, seed = 1) {
  set.seed(seed)

  n <- length(x)

  y <- 0.1 * x +
    2 * sin(x / 8) +
    rnorm(n, sd = 0.5)

  w <- runif(n, 0.5, 2)

  fit_ours <- .pb_fit_default(
    y = y,
    x = x,
    w = w
  )

  fit_full <- gamlss::gamlss(
    y ~ gamlss::pb(x),
    weights = w,
    family = gamlss.dist::NO,
    trace = FALSE
  )

  pred_full <- as.vector(predict(fit_full, type = "response"))

  data.frame(
    n = n,
    n_unique = length(unique(x)),
    max_abs_ours_vs_full = max(abs(fit_ours$fitted.values - pred_full)),
    mean_abs_ours_vs_full = mean(abs(fit_ours$fitted.values - pred_full)),
    cor_ours_vs_full = cor(fit_ours$fitted.values, pred_full),
    lambda_ours = fit_ours$lambda
  )
}


# -------------------------------------------------------------------
# 5. Execute tests
# -------------------------------------------------------------------

cat("\n============================================================\n")
cat("Testing .pb_fit_default() against gamlss:::gamlss.pb()\n")
cat("============================================================\n")

suite_res <- run_pb_equivalence_suite()

cat("\n============================================================\n")
cat("Summary table\n")
cat("============================================================\n")

print(suite_res)

cat("\n============================================================\n")
cat("Optional comparison against full gamlss() model\n")
cat("This is not expected to be identical.\n")
cat("============================================================\n")

full_res <- compare_against_full_gamlss(
  x = seq(-45, 45, length.out = 300),
  seed = 100
)

print(full_res)


# -------------------------------------------------------------------
# 6. Simple pass/fail checks
# -------------------------------------------------------------------

tol <- 1e-6

cat("\n============================================================\n")
cat("Pass/fail checks against gamlss:::gamlss.pb()\n")
cat("Tolerance:", tol, "\n")
cat("============================================================\n")

pass <- all(suite_res$max_abs_fitted < tol) &&
  all(abs(suite_res$lambda_diff) < tol) &&
  all(abs(suite_res$edf_diff) < tol)

if (pass) {
  cat("PASS: helper matches gamlss:::gamlss.pb() within tolerance.\n")
} else {
  cat("FAIL or approximate match only.\n")
  cat("Largest fitted-value difference:\n")
  print(max(suite_res$max_abs_fitted))

  cat("Largest lambda difference:\n")
  print(max(abs(suite_res$lambda_diff)))

  cat("Largest edf difference:\n")
  print(max(abs(suite_res$edf_diff)))
}


# -------------------------------------------------------------------
# 7. Full default-case location-scale fit:
#    y ~ pb(x)
#    sigma ~ abs(x)
#    family = NO
#
# This mimics the default GAMLSS RS-style loop for the specific case
# used in remove_cardinal_biases().
# -------------------------------------------------------------------

.fit_pb_sigma_rs_default <- function(y, x, outlier = NULL,
                                     max_iter = 50,
                                     tol = 1e-5,
                                     trace = FALSE) {
  y <- as.vector(y)
  x <- as.vector(x)

  stopifnot(length(y) == length(x))

  if (is.null(outlier)) {
    base_w <- rep(1, length(y))
  } else {
    base_w <- 1 - as.numeric(outlier)
  }

  stopifnot(length(base_w) == length(y))

  x_abs <- abs(x)
  X_sigma <- model.matrix(~ x_abs)

  # Starting values.
  # Only non-zero-weight observations should affect starts.
  keep <- base_w > 0

  mu0 <- weighted.mean(y[keep], base_w[keep])

  mu <- rep(mu0, length(y))

  sig0 <- sqrt(weighted.mean((y[keep] - mu0)^2, base_w[keep]))

  if (!is.finite(sig0) || sig0 <= 0) {
    sig0 <- stats::sd(y[keep])
  }

  if (!is.finite(sig0) || sig0 <= 0) {
    sig0 <- 1
  }

  eta_sigma <- rep(log(sig0), length(y))
  beta_sigma <- c(log(sig0), 0)

  dev_old <- Inf
  mu_fit <- NULL

  for (iter in seq_len(max_iter)) {
    # ------------------------------------------------------------
    # 1. Update mu using current sigma
    # ------------------------------------------------------------

    sigma <- exp(eta_sigma)

    # For Normal with known sigma, Fisher scoring/WLS for mu
    # uses weights proportional to 1 / sigma^2.
    w_mu <- base_w / sigma^2

    mu_fit <- .pb_fit_default(
      y = y,
      x = x,
      w = w_mu
    )

    mu <- mu_fit$fitted.values

    # ------------------------------------------------------------
    # 2. Update sigma using Fisher scoring for eta = log(sigma)
    # ------------------------------------------------------------

    sigma <- exp(eta_sigma)

    resid2_scaled <- ((y - mu)^2) / sigma^2

    # For Normal with log(sigma) link:
    # score wrt eta is ((y - mu)^2 / sigma^2 - 1)
    # expected information is 2
    z_sigma <- eta_sigma + (resid2_scaled - 1) / 2
    w_sigma <- base_w * 2

    lm_sigma <- lm.wfit(
      x = X_sigma,
      y = z_sigma,
      w = w_sigma
    )

    beta_sigma <- as.vector(lm_sigma$coefficients)

    # If a coefficient is NA because of singularity, keep old one.
    if (anyNA(beta_sigma)) {
      beta_old <- coef(lm.wfit(
        x = X_sigma[, 1, drop = FALSE],
        y = z_sigma,
        w = w_sigma
      ))

      beta_sigma <- c(as.vector(beta_old)[1], 0)
    }

    eta_sigma_new <- as.vector(X_sigma %*% beta_sigma)

    # Basic numerical guard.
    eta_sigma_new <- pmin(pmax(eta_sigma_new, -30), 30)

    sigma_new <- exp(eta_sigma_new)

    # ------------------------------------------------------------
    # 3. Check convergence via global deviance
    # ------------------------------------------------------------

    loglik <- sum(
      base_w * dnorm(y, mean = mu, sd = sigma_new, log = TRUE)
    )

    dev <- -2 * loglik

    if (trace) {
      message(
        "iter = ", iter,
        " dev = ", signif(dev, 8),
        " lambda = ", signif(mu_fit$lambda, 6),
        " edf_mu = ", signif(mu_fit$edf, 6),
        " beta_sigma = ",
        paste(signif(beta_sigma, 6), collapse = ", ")
      )
    }

    if (is.finite(dev_old)) {
      rel_change <- abs(dev_old - dev) / (abs(dev_old) + 0.1)

      if (rel_change < tol) {
        eta_sigma <- eta_sigma_new
        dev_old <- dev
        break
      }
    }

    eta_sigma <- eta_sigma_new
    dev_old <- dev
  }

  sigma <- exp(eta_sigma)

  logLik <- sum(
    base_w * dnorm(y, mean = mu, sd = sigma, log = TRUE)
  )

  list(
    pred = mu,
    fitted.values = mu,
    pred_sigma = sigma,
    sigma = sigma,
    eta_sigma = eta_sigma,
    beta_sigma = beta_sigma,
    coef_sigma = beta_sigma,
    lambda = mu_fit$lambda,
    edf_mu = mu_fit$edf,
    logLik = logLik,
    deviance = -2 * logLik,
    iterations = iter,
    converged = iter < max_iter
  )
}


# -------------------------------------------------------------------
# 8. Compare the full default helper against gamlss()
# -------------------------------------------------------------------

compare_full_default <- function(x, seed = 1, outlier_rate = 0,
                                 trace_ours = FALSE,
                                 trace_gamlss = FALSE) {
  set.seed(seed)

  n <- length(x)

  # Simulated data with smooth mean and sigma depending on abs(x).
  true_mu <- 0.1 * x + 2 * sin(x / 8)
  true_sigma <- exp(-0.2 + 0.01 * abs(x))

  y <- true_mu + rnorm(n, sd = true_sigma)

  outlier <- runif(n) < outlier_rate
  w <- 1 - as.numeric(outlier)

  fit_ours <- .fit_pb_sigma_rs_default(
    y = y,
    x = x,
    outlier = outlier,
    trace = trace_ours
  )

  fit_g <- gamlss::gamlss(
    y ~ pb(x),
    sigma.formula = ~ abs(x),
    weights = w,
    family = gamlss.dist::NO,
    trace = trace_gamlss
  )
  pred_g <- as.vector(predict(fit_g, type = "response"))
  sigma_g <- as.vector(predict(fit_g, what = "sigma", type = "response"))

  beta_sigma_g <- coef(fit_g, what = "sigma")

  out <- data.frame(
    n = length(x),
    n_unique = length(unique(x)),

    max_abs_mu = max(abs(fit_ours$pred - pred_g)),
    mean_abs_mu = mean(abs(fit_ours$pred - pred_g)),
    cor_mu = cor(fit_ours$pred, pred_g),

    max_abs_sigma = max(abs(fit_ours$pred_sigma - sigma_g)),
    mean_abs_sigma = mean(abs(fit_ours$pred_sigma - sigma_g)),
    cor_sigma = cor(fit_ours$pred_sigma, sigma_g),

    beta_sigma_ours_0 = fit_ours$beta_sigma[1],
    beta_sigma_ours_1 = fit_ours$beta_sigma[2],

    beta_sigma_g_0 = beta_sigma_g[1],
    beta_sigma_g_1 = beta_sigma_g[2],

    beta_sigma_diff_0 = fit_ours$beta_sigma[1] - beta_sigma_g[1],
    beta_sigma_diff_1 = fit_ours$beta_sigma[2] - beta_sigma_g[2],

    logLik_ours = fit_ours$logLik,
    logLik_gamlss = as.numeric(logLik(fit_g)),
    logLik_diff = fit_ours$logLik - as.numeric(logLik(fit_g)),

    iterations_ours = fit_ours$iterations,
    iterations_gamlss = fit_g$iter
  )

  invisible(list(
    summary = out,
    x = x,
    y = y,
    w = w,
    outlier = outlier,
    ours = fit_ours,
    gamlss = fit_g,
    pred_gamlss = pred_g,
    sigma_gamlss = sigma_g
  ))
}


run_full_default_suite <- function(outlier_rate = 0) {
  set.seed(123)

  cases <- list(
    dense_even = seq(-45, 45, length.out = 300),

    repeated_grid = rep(seq(-45, 45, length.out = 31), each = 10),

    small_n = seq(-45, 45, length.out = 80),

    integer_sample = sample(seq(-45, 45, by = 1), 300, replace = TRUE),

    narrow_range = seq(-10, 10, length.out = 300),

    asymmetric = sort(runif(300, -20, 45))
  )

  res <- lapply(seq_along(cases), function(i) {
    cat("\n--- Full default case:", names(cases)[i], "---\n")

    tmp <- compare_full_default(
      x = cases[[i]],
      seed = 100 + i,
      outlier_rate = outlier_rate
    )

    print(tmp$summary)

    cbind(
      case = names(cases)[i],
      tmp$summary
    )
  })

  do.call(rbind, res)
}


# -------------------------------------------------------------------
# 9. Execute full-model comparison
# -------------------------------------------------------------------

cat("\n============================================================\n")
cat("Testing .fit_pb_sigma_rs_default() against full gamlss()\n")
cat("Model: y ~ pb(x), sigma ~ abs(x), family = NO\n")
cat("============================================================\n")

full_suite_res <- run_full_default_suite(outlier_rate = 0)

cat("\n============================================================\n")
cat("Full-model summary table\n")
cat("============================================================\n")

print(full_suite_res)


# -------------------------------------------------------------------
# 10. Optional: full-model comparison with zero-weight outliers
# -------------------------------------------------------------------

cat("\n============================================================\n")
cat("Testing full model with zero-weight outliers\n")
cat("============================================================\n")

full_suite_res_outliers <- run_full_default_suite(outlier_rate = 0.05)

cat("\n============================================================\n")
cat("Full-model summary table with outliers\n")
cat("============================================================\n")

print(full_suite_res_outliers)


# -------------------------------------------------------------------
# 11. Diagnostic plotting helper for one failed/mismatching case
# -------------------------------------------------------------------

plot_full_default_comparison <- function(obj) {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  par(mfrow = c(2, 2))

  plot(
    obj$x,
    obj$y,
    main = "Data and fitted mu",
    xlab = "x",
    ylab = "y",
    pch = 16,
    cex = 0.5
  )
  lines(obj$x, obj$ours$pred, lwd = 2)
  lines(obj$x, obj$pred_gamlss, lwd = 2, lty = 2)
  legend(
    "topleft",
    legend = c("ours", "gamlss"),
    lty = c(1, 2),
    lwd = 2,
    bty = "n"
  )

  plot(
    obj$x,
    obj$ours$pred - obj$pred_gamlss,
    main = "mu difference",
    xlab = "x",
    ylab = "ours - gamlss",
    pch = 16,
    cex = 0.5
  )
  abline(h = 0, lty = 2)

  plot(
    obj$x,
    obj$ours$pred_sigma,
    main = "Fitted sigma",
    xlab = "x",
    ylab = "sigma",
    type = "l",
    lwd = 2
  )
  lines(obj$x, obj$sigma_gamlss, lwd = 2, lty = 2)
  legend(
    "topleft",
    legend = c("ours", "gamlss"),
    lty = c(1, 2),
    lwd = 2,
    bty = "n"
  )

  plot(
    obj$x,
    obj$ours$pred_sigma - obj$sigma_gamlss,
    main = "sigma difference",
    xlab = "x",
    ylab = "ours - gamlss",
    pch = 16,
    cex = 0.5
  )
  abline(h = 0, lty = 2)
}


# Example manual diagnostic:
#
# diag_obj <- compare_full_default(
#   x = seq(-45, 45, length.out = 300),
#   seed = 101
# )
# plot_full_default_comparison(diag_obj)