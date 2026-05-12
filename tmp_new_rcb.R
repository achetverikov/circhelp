get_boundary_preds <- function(
    group,
    data,
    space,
    reassign_range,
    gam_ctrl,
    poly_deg,
    angle_diff_fun,
    weights = NULL,
    gamlss_fun = gamlss::gamlss
) {
  gr_var <- outlier <- err <- x_var <- dc_var <- center_x <-
    dist_to_card <- bin_boundary_left <- bin_boundary_right <-
    bin_range <- dist_to_bin_centre <- row_i <- at_the_boundary <-
    dist_to_boundary <- dist_to_boundary_norm <- new_weight <- weight <-
    pred <- pred_sigma <- resid_at_boundaries <- NULL

  cur_df <- data[
    gr_var == group & outlier == FALSE,
    .(
      err,
      x_var,
      dc_var,
      dist_to_bin_centre = angle_diff_fun(x_var, center_x),
      weight = NULL,
      adc = abs(dist_to_card),
      center_x,
      bin_boundary_left,
      bin_boundary_right,
      bin_range
    )
  ]

  curr_bin_range <- cur_df$bin_range[[1]]
  curr_bin_center <- cur_df$center_x[[1]]

  data[, dist_to_bin_centre := angle_diff_fun(x_var, curr_bin_center)]

  data_incl_boundaries <- data[
    outlier == FALSE &
      abs(dist_to_bin_centre) < (curr_bin_range / 2 + reassign_range + 1e-12),
    .(
      row_i,
      err,
      x_var,
      dc_var,
      gr_var,
      dist_to_bin_centre,
      at_the_boundary,
      center_x
    )
  ]

  data_incl_boundaries[, dist_to_boundary :=
    abs(dist_to_bin_centre) - curr_bin_range / 2]

  data_incl_boundaries[, dist_to_boundary_norm :=
    (dist_to_boundary + reassign_range) / (2 * reassign_range)]

  if (!missing(weights) && !is.null(weights)) {
    data_incl_boundaries[, gr_var := group]
    data_incl_boundaries[
      weights,
      `:=`(weight = new_weight),
      on = .(row_i, gr_var)
    ]
  } else {
    data_incl_boundaries[, weight := ifelse(
      at_the_boundary == FALSE,
      1,
      ifelse(gr_var == group, 0.75, 0.25)
    )]
  }

  fit <- gamlss_fun(
    err ~ dist_to_bin_centre,
    ~ abs(dist_to_bin_centre),
    data = data_incl_boundaries,
    weights = weight,
    control = gam_ctrl
  )

  data_incl_boundaries[, pred :=
    predict(fit, type = "response")]

  data_incl_boundaries[, pred_sigma :=
    predict(fit, what = "sigma", type = "response")]

  data_incl_boundaries[, resid_at_boundaries := err - pred]

  data_incl_boundaries[
    ,
    .(
      row_i,
      at_the_boundary,
      x_var,
      dc_var,
      dist_to_bin_centre,
      err,
      pred,
      resid_at_boundaries,
      dist_to_boundary,
      dist_to_boundary_norm,
      weight,
      pred_sigma
    )
  ]
}


remove_cardinal_biases <- function(
    err,
    x,
    space = "180",
    bias_type = "fit",
    plots = "hide",
    poly_deg = 4,
    var_sigma = TRUE,
    var_sigma_poly_deg = 4,
    reassign_at_boundaries = TRUE,
    reassign_range = 2,
    break_points = NULL,
    init_outliers = NULL,
    debug = FALSE,
    do_plots = NULL,
    gamlss_fun = gamlss::gamlss
) {
  outlier <- dist_to_card <- dist_to_obl <- logLik <- x_var <-
    min_bp_i <- center_x <- dc_var <- gr_var <- min_boundary_i <-
    min_boundary_dist <- bin_range <- bin_boundary_left <-
    bin_boundary_right <- at_the_boundary <- row_i <- likelihood <-
    dnorm <- pred <- pred_sigma <- new_weight <- i.gr_var <-
    dist_to_bin_centre <- coef <- predict <- bias <- pred_lin <-
    be_c <- which_bin <- center_y <- outlier_f <- coef_sigma_int <-
    coef_sigma_slope <- NULL

  if (!(bias_type %in% c("fit", "card", "obl", "custom"))) {
    stop("`bias_type` should be 'fit','card', 'obl', or 'custom'")
  }

  if (bias_type == "custom" && missing(break_points)) {
    stop("If 'bias_type' is set to 'custom', you need to specify 'break_points'")
  }

  if (any(is.na(x)) || any(is.na(err))) {
    stop("There are NAs in x or err. Please remove missing values before running the function.")
  }

  if (!missing(do_plots)) {
    warning(
      "\nYou have supplied 'do_plots' argument, it is now deprecated in favor of a 'plots' argument"
    )

    if (do_plots) {
      plots <- "show"
    }
  } else {
    if (!(plots %in% c("show", "hide", "return"))) {
      stop("`plots` should be 'show','hide', or 'return'")
    }
  }

  if (space == "180") {
    x <- angle_diff_180(x, 0)
    x2 <- angle_diff_180_45(x, 0)

    obl_groups <- cut(
      x,
      breaks = seq(-90, 90, 90),
      include.lowest = TRUE
    )

    obl_bin_centers <- c(-45, 45)

    card_groups <- cut(
      x2,
      breaks = seq(-45, 180 - 45, 90),
      include.lowest = TRUE
    )

    card_bin_centers <- c(0, 90)

    angle_diff_fun <- angle_diff_180
    circ_sd_fun <- circ_sd_180
  } else if (space == "360") {
    x <- angle_diff_360(x, 0)
    x2 <- (x + 45) %% 360 - 45

    obl_groups <- cut(
      x,
      breaks = seq(-180, 180, 90),
      include.lowest = TRUE
    )

    obl_bin_centers <- seq(-135, 135, 90)

    card_groups <- cut(
      x2,
      breaks = seq(-45, 360 - 45, 90),
      include.lowest = TRUE
    )

    card_bin_centers <- seq(0, 270, 90)

    angle_diff_fun <- angle_diff_360
    circ_sd_fun <- circ_sd_360
  } else {
    stop("`space` argument should be 180 or 360.")
  }

  if (debug) {
    cat("N observation per group assuming cardinal bins: \n")
    print(table(card_groups))
    cat("N observation per group assuming oblique bins: \n")
    print(table(obl_groups))
  }

  for_fit <- data.table::data.table(
    x = x,
    x2 = x2,
    err = err,
    card_groups = card_groups,
    obl_groups = obl_groups
  )

  if (missing(init_outliers) || is.null(init_outliers)) {
    for_fit[, outlier := abs(err) > (3 * circ_sd_fun(err, na.rm = TRUE))]
  } else {
    for_fit[, outlier := init_outliers]
  }

  for_fit[, dist_to_card := angle_diff_90(x2, 0)]
  for_fit[, dist_to_obl := angle_diff_90(x, 45)]

  gam_ctrl <- gamlss::gamlss.control(trace = FALSE)

  if (debug) {
    cat("Computing bins to group the data...\n")
  }

  if (bias_type == "fit") {
    if (var_sigma) {
      sigma_formula <- ~ abs(dist_to_card)

      ll1 <- sum(
        for_fit[
          outlier == FALSE,
          logLik(gamlss_fun(
            err ~ poly(dist_to_card, var_sigma_poly_deg),
            sigma_formula,
            data = .SD,
            control = gam_ctrl
          )),
          by = .(card_groups)
        ]$V1
      )

      ll2 <- sum(
        for_fit[
          outlier == FALSE,
          logLik(gamlss_fun(
            err ~ poly(dist_to_obl, var_sigma_poly_deg),
            sigma_formula,
            data = .SD,
            control = gam_ctrl
          )),
          by = .(obl_groups)
        ]$V1
      )
    } else {
      ll1 <- sum(
        for_fit[
          outlier == FALSE,
          logLik(MASS::rlm(err ~ poly(x2, poly_deg))),
          by = .(card_groups)
        ]$V1
      )

      ll2 <- sum(
        for_fit[
          outlier == FALSE,
          logLik(MASS::rlm(err ~ poly(x, poly_deg))),
          by = .(obl_groups)
        ]$V1
      )
    }

    if (debug) {
      cat(sprintf("LL for bias type 1: %.2f, LL for bias type 2: %.2f\n", ll1, ll2))
    }

    if (ll1 >= ll2) {
      bias_type <- "card"
    } else {
      bias_type <- "obl"
    }
  }

  if (bias_type == "obl") {
    break_points <- card_bin_centers
  } else if (bias_type == "card") {
    break_points <- obl_bin_centers
  }

  for_fit$pred_sigma <- NA_real_
  for_fit$coef_sigma_int <- NA_real_
  for_fit$coef_sigma_slope <- NA_real_

  break_points <- sort(break_points)

  bin_boundaries <- c(break_points[length(break_points)], break_points)

  bin_centers <- break_points +
    (data.table::shift(
      break_points,
      1,
      fill = break_points[length(break_points)]
    ) - break_points) / 2

  bin_centers[1] <- angle_diff_fun(
    break_points[length(break_points)] +
      (break_points[1] + as.numeric(space) - break_points[length(break_points)]) / 2,
    0
  )

  bin_width <- abs(
    angle_diff_fun(
      bin_boundaries,
      data.table::shift(bin_boundaries, 1)
    )[2:length(bin_boundaries)]
  )

  bin_width[1] <- break_points[1] +
    as.numeric(space) -
    break_points[length(break_points)]

  bin_labels <- sapply(2:length(bin_boundaries), function(i) {
    sprintf("[%.2f, %.2f]", bin_boundaries[i - 1], bin_boundaries[i])
  })

  bin_labels <- factor(bin_labels, levels = bin_labels)

  for_fit[, x_var := x]

  get_bin_i <- function(x, bin_centers, bin_width) {
    within_bin <- sapply(seq_along(bin_centers), function(i) {
      abs(angle_diff_fun(x, bin_centers[i])) <= (bin_width[i] / 2)
    })

    max.col(within_bin, "first")
  }

  for_fit[, min_bp_i := get_bin_i(x, bin_centers, bin_width)]
  for_fit[, center_x := bin_centers[min_bp_i]]
  for_fit[, dc_var := angle_diff_fun(x, center_x)]
  for_fit[, gr_var := bin_labels[min_bp_i]]

  if (plots == "show" && debug == TRUE) {
    p_boundaries <- ggplot2::ggplot(
      for_fit,
      ggplot2::aes(x = .data$x, y = .data$err, color = .data$gr_var)
    ) +
      ggplot2::geom_point() +
      ggplot2::geom_vline(xintercept = angle_diff_fun(break_points, 0)) +
      ggplot2::geom_vline(
        color = "blue",
        xintercept = angle_diff_fun(bin_centers, 0)
      )

    print(p_boundaries)
  }

  for_fit[
    ,
    min_boundary_i := apply(
      sapply(break_points, function(bp) abs(angle_diff_fun(x, bp))),
      1,
      which.min
    )
  ]

  for_fit[
    ,
    min_boundary_dist := angle_diff_fun(x, break_points[min_boundary_i])
  ]

  for_fit[, bin_range := bin_width[min_bp_i]]
  for_fit[, bin_boundary_left := bin_boundaries[min_bp_i]]
  for_fit[, bin_boundary_right := bin_boundaries[min_bp_i + 1]]

  if (reassign_at_boundaries) {
    for_fit[
      ,
      at_the_boundary := (abs(min_boundary_dist) - reassign_range) < 1e-12
    ]
  } else {
    for_fit[, at_the_boundary := FALSE]
  }

  likelihoods <- numeric(0)

  if (var_sigma) {
    if (reassign_at_boundaries) {
      if (debug) {
        cat("Reassigning points at the boundaries...\n")
      }

      if (any(for_fit[, unique(bin_range)] < (2 * reassign_range))) {
        stop("Reassignment range too large compared to bin sizes")
      }

      for_fit[, row_i := .I]

      resid_at_boundaries <- for_fit[
        outlier == FALSE,
        get_boundary_preds(
          gr_var,
          copy(for_fit[outlier == FALSE]),
          space,
          reassign_range,
          gam_ctrl,
          1,
          angle_diff_fun,
          gamlss_fun = gamlss_fun
        ),
        by = .(gr_var)
      ]

      resid_at_boundaries[
        ,
        likelihood := stats::dnorm(err, pred, pred_sigma, log = FALSE)
      ]

      resid_at_boundaries[
        ,
        new_weight := ifelse(
          at_the_boundary == FALSE,
          1,
          likelihood / sum(likelihood)
        ),
        by = .(err, x_var)
      ]

      cur_weights <- resid_at_boundaries[at_the_boundary == TRUE]$new_weight
      stable_weights <- 0

      for (rep_n in 1:10) {
        weight_dt <- resid_at_boundaries[
          ,
          .(row_i, gr_var, new_weight)
        ]

        resid_at_boundaries <- resid_at_boundaries[
          ,
          get_boundary_preds(
            gr_var,
            copy(for_fit[outlier == FALSE]),
            space,
            reassign_range,
            gam_ctrl,
            ifelse(rep_n > 2, poly_deg, 1),
            angle_diff_fun,
            weights = weight_dt,
            gamlss_fun = gamlss_fun
          ),
          by = .(gr_var)
        ]

        resid_at_boundaries[
          ,
          likelihood := stats::dnorm(err, pred, pred_sigma, log = FALSE)
        ]

        resid_at_boundaries[
          ,
          new_weight := ifelse(
            at_the_boundary == FALSE,
            1,
            likelihood / sum(likelihood)
          ),
          by = .(err, x_var)
        ]

        resid_at_boundaries_c <- data.table::dcast(
          resid_at_boundaries[at_the_boundary == TRUE],
          row_i ~ gr_var,
          value.var = "resid_at_boundaries"
        )

        resid_at_boundaries_c <- resid_at_boundaries_c[
          ,
          gr_var := names(.SD)[max.col(replace(-abs(.SD), is.na(.SD), -Inf))],
          .SDcols = 2:ncol(resid_at_boundaries_c)
        ]

        prev_weights <- cur_weights
        cur_weights <- resid_at_boundaries[at_the_boundary == TRUE]$new_weight
        weights_change <- sum(abs(cur_weights - prev_weights))

        if (debug) {
          cat(sprintf(
            "Reassignment step: %i; change in weights: %.5f\n",
            rep_n,
            weights_change
          ))
        }

        for_fit[
          resid_at_boundaries_c,
          `:=`(gr_var = i.gr_var),
          on = .(row_i)
        ]

        for_fit[, center_x := bin_centers[as.numeric(gr_var)]]
        for_fit[, x_var := center_x + angle_diff_fun(x_var, center_x)]
        for_fit[, dc_var := angle_diff_fun(x, center_x)]
        for_fit[, bin_range := bin_width[as.numeric(gr_var)]]
        for_fit[, bin_boundary_left := bin_boundaries[as.numeric(gr_var)]]
        for_fit[, bin_boundary_right := bin_boundaries[as.numeric(gr_var) + 1]]

        if (rep_n > 1) {
          if (weights_change < 0.01) {
            stable_weights <- stable_weights + 1
          } else {
            stable_weights <- 0
          }

          if (stable_weights > 3) {
            if (debug) {
              cat("Reassignment stopped at stable weights\n")
            }
            break
          }
        }
      }
    }

    for_fit[, dist_to_bin_centre := angle_diff_fun(x_var, center_x)]
    for_fit[, x_var := center_x + angle_diff_fun(x_var, center_x)]
    for_fit[, dc_var := angle_diff_fun(x, center_x)]

    if (debug) {
      cat("Computing final fits...\n")
    }

    for (cg in unique(for_fit$gr_var)) {
      cur_df <- for_fit[
        gr_var == cg,
        .(
          err,
          x_var,
          dist_to_bin_centre,
          dc_var,
          outlier,
          dist_to_card
        )
      ]

      fit <- gamlss_fun(
        err ~ pb(dist_to_bin_centre),
        ~ abs(dist_to_bin_centre),
        data = cur_df,
        weights = 1 - as.numeric(cur_df$outlier),
        control = gam_ctrl
      )

      if (debug) {
        cat("Fitted GAMLSS model coefficients\n")
        print(coef(fit))
      }

      for_fit[gr_var == cg, pred := predict(fit, type = "response")]

      for_fit[
        gr_var == cg,
        pred_sigma := predict(fit, what = "sigma", type = "response")
      ]

      for_fit[gr_var == cg, bias := err * sign(pred)]

      if (debug) {
        p_pred <- ggplot2::ggplot(
          for_fit[gr_var == cg],
          ggplot2::aes(x = .data$dist_to_bin_centre, y = .data$err)
        ) +
          ggplot2::geom_point() +
          ggplot2::geom_line(ggplot2::aes(y = .data$pred))

        print(p_pred)
      }

      sigma_coef <- coef(fit, what = "sigma")

      for_fit[
        gr_var == cg,
        c("coef_sigma_int", "coef_sigma_slope") :=
          as.list(as.numeric(sigma_coef))
      ]

      likelihoods <- c(likelihoods, as.numeric(logLik(fit)))
    }
  } else {
    for_fit[
      ,
      pred := predict(
        MASS::rlm(err ~ poly(x_var, poly_deg), data = .SD[outlier == FALSE]),
        newdata = .SD[, .(x_var)]
      ),
      by = .(card_groups)
    ]

    for_fit[, pred_sigma := NA_real_]
    for_fit[, bias := err * sign(pred)]
  }

  for_fit[
    ,
    pred_lin := MASS::rlm(err ~ x_var)$fitted.values,
    by = .(gr_var)
  ]

  for_fit[, be_c := err - pred]
  for_fit[, which_bin := as.numeric(gr_var)]

  for_fit[
    ,
    center_y := predict(
      MASS::rlm(err ~ x_var),
      newdata = data.frame(x_var = center_x)
    ),
    by = .(gr_var)
  ]

  if (var_sigma) {
    for_fit[, outlier := abs(be_c) > 3 * pred_sigma]
  } else {
    for_fit[, outlier := abs(be_c) > (3 * circ_sd_fun(be_c))]
  }

  if (plots %in% c("show", "return")) {
    for_fit[
      ,
      outlier_f := factor(ifelse(outlier, "Outlier", "Non-outlier"))
    ]

    sd_val <- for_fit[, circ_sd_fun(err)]

    plots_obj <- make_plots_of_biases(for_fit, poly_deg, sd_val)

    if (plots == "show") {
      print(plots_obj)
    } else {
      return(plots_obj)
    }
  }

  for_fit[
    ,
    .(
      is_outlier = as.numeric(outlier),
      pred,
      be_c,
      which_bin,
      bias,
      bias_type,
      pred_lin,
      pred_sigma,
      coef_sigma_int,
      coef_sigma_slope,
      shifted_x = x_var,
      total_log_lik = sum(likelihoods)
    )
  ]
}

# -------------------------------------------------------------------
# Minimal plain-R replacement for the subset of gamlss::gamlss()
# needed by remove_cardinal_biases()
#
# Supports:
#   y ~ pb(x), sigma ~ abs(x)
#   y ~ x,     sigma ~ abs(x)
#   y ~ poly(x, degree), sigma ~ abs(x)
#
# Returns an object that supports:
#   predict(fit, type = "response")
#   predict(fit, what = "sigma", type = "response")
#   coef(fit)
#   coef(fit, what = "sigma")
#   logLik(fit)
# -------------------------------------------------------------------


# -------------------------------------------------------------------
# 1. Default gamlss::pb() basis
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


# -------------------------------------------------------------------
# 2. Default gamlss.pb() penalized smoother step
# -------------------------------------------------------------------

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
# 3. Detect mean model type
# -------------------------------------------------------------------

.gamlss_default_mu_type <- function(formula) {
  f_txt <- paste(deparse(formula), collapse = "")

  if (grepl("pb\\s*\\(", f_txt)) {
    return("pb")
  }

  if (grepl("poly\\s*\\(", f_txt)) {
    return("poly")
  }

  "linear"
}


.extract_response_name <- function(formula) {
  all.vars(formula[[2]])[1]
}


.extract_single_rhs_var <- function(formula) {
  vars <- all.vars(formula[[3]])
  if (length(vars) < 1) {
    stop("Could not identify RHS variable in formula.")
  }
  vars[1]
}


.extract_poly_degree <- function(formula, data = NULL, env = parent.frame()) {
  poly_call <- .find_call(formula[[3]], "poly")

  if (is.null(poly_call)) {
    stop("Could not find poly() in formula.", call. = FALSE)
  }

  if (length(poly_call) < 3) {
    stop("poly() call does not contain a degree argument.", call. = FALSE)
  }

  degree_expr <- poly_call[[3]]

  degree <- tryCatch(
    eval(degree_expr, envir = data, enclos = environment(formula)),
    error = function(e1) {
      eval(degree_expr, envir = env)
    }
  )

  degree <- as.integer(degree)

  if (length(degree) != 1 || is.na(degree)) {
    stop("poly() degree must evaluate to a single integer.", call. = FALSE)
  }

  degree
}

# -------------------------------------------------------------------
# 4. Build sigma model matrix
# -------------------------------------------------------------------

.make_sigma_matrix_default <- function(sigma.formula, data) {
  # In remove_cardinal_biases() this is always:
  #   ~ abs(some_x)
  #
  # model.matrix() handles this correctly and gives:
  #   (Intercept), abs(some_x)

  stats::model.matrix(sigma.formula, data = data)
}


# -------------------------------------------------------------------
# 5. Plain-R RS loop for:
#      y ~ pb(x), sigma ~ abs(x)
# -------------------------------------------------------------------

.fit_pb_sigma_rs_default <- function(y, x, X_sigma, weights,
                                     max_iter = 50,
                                     tol = 1e-5,
                                     trace = FALSE) {
  y <- as.vector(y)
  x <- as.vector(x)
  base_w <- as.vector(weights)

  stopifnot(length(y) == length(x))
  stopifnot(length(y) == length(base_w))
  stopifnot(nrow(X_sigma) == length(y))

  keep <- base_w > 0

  mu0 <- stats::weighted.mean(y[keep], base_w[keep])
  mu <- rep(mu0, length(y))

  sig0 <- sqrt(stats::weighted.mean((y[keep] - mu0)^2, base_w[keep]))

  if (!is.finite(sig0) || sig0 <= 0) {
    sig0 <- stats::sd(y[keep])
  }

  if (!is.finite(sig0) || sig0 <= 0) {
    sig0 <- 1
  }

  eta_sigma <- rep(log(sig0), length(y))
  beta_sigma <- rep(0, ncol(X_sigma))
  beta_sigma[1] <- log(sig0)

  dev_old <- Inf
  mu_fit <- NULL

  for (iter in seq_len(max_iter)) {
    sigma <- exp(eta_sigma)

    w_mu <- base_w / sigma^2

    mu_fit <- .pb_fit_default(
      y = y,
      x = x,
      w = w_mu
    )

    mu <- mu_fit$fitted.values

    sigma <- exp(eta_sigma)
    resid2_scaled <- ((y - mu)^2) / sigma^2

    z_sigma <- eta_sigma + (resid2_scaled - 1) / 2
    w_sigma <- base_w * 2

    lm_sigma <- stats::lm.wfit(
      x = X_sigma,
      y = z_sigma,
      w = w_sigma
    )

    beta_new <- as.vector(lm_sigma$coefficients)

    if (anyNA(beta_new)) {
      beta_new[is.na(beta_new)] <- beta_sigma[is.na(beta_new)]
    }

    beta_sigma <- beta_new
    eta_sigma_new <- as.vector(X_sigma %*% beta_sigma)

    eta_sigma_new <- pmin(pmax(eta_sigma_new, -30), 30)
    sigma_new <- exp(eta_sigma_new)

    loglik <- sum(
      base_w * stats::dnorm(y, mean = mu, sd = sigma_new, log = TRUE)
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
    base_w * stats::dnorm(y, mean = mu, sd = sigma, log = TRUE)
  )

  list(
    mu = mu,
    sigma = sigma,
    beta_mu = mu_fit$beta,
    beta_sigma = beta_sigma,
    lambda = mu_fit$lambda,
    edf_mu = mu_fit$edf,
    logLik = logLik,
    deviance = -2 * logLik,
    iterations = iter,
    converged = iter < max_iter
  )
}


# -------------------------------------------------------------------
# 5b. Fast P-spline lambda-EM: same basis + EM updates as default,
#     but precomputes QR + eigen(R^{-T} S R^{-1}) so each EM
#     iteration is O(p) instead of O(n p^2 + p^3) SVD.
#     Numerically equivalent to .pb_fit_default for the same weights.
# -------------------------------------------------------------------

.pb_fit_fast <- function(y, x, w = rep(1, length(y)),
                         lambda    = 10,
                         max_iter  = 50,
                         tol       = 1e-7) {
  y <- as.vector(y)
  x <- as.vector(x)
  w <- as.vector(w)

  basis <- .pb_basis_default(x)      # same basis as default
  X     <- basis$X
  D     <- basis$D                   # (p-2) x p difference matrix
  p     <- ncol(X)
  S     <- crossprod(D)              # p x p penalty matrix

  n_eff        <- sum(w > 0)
  penalty_order <- 2L

  # --- one-time decomposition (per RS iteration, fixed weights) ---
  sqw  <- sqrt(w)
  qr_X <- qr(X * sqw, tol = .Machine$double.eps^0.8)
  Q    <- qr.Q(qr_X)
  R    <- qr.R(qr_X)
  Rinv <- backsolve(R, diag(p))

  K   <- crossprod(Rinv, S %*% Rinv)
  K   <- (K + t(K)) / 2
  eig <- eigen(K, symmetric = TRUE)
  U   <- eig$vectors
  lk  <- pmax(eig$values, 0)

  z_w <- as.numeric(crossprod(Q, sqw * y))
  qq  <- as.numeric(crossprod(U, z_w))
  y_y <- sum(w * y^2)

  # --- EM loop: each iteration is O(p) ---
  sig2 <- NA_real_

  for (iter in seq_len(max_iter)) {
    d    <- 1 + lambda * lk
    edf  <- sum(1 / d)
    wrss <- max(y_y - sum(qq^2 * (2 / d - 1 / d^2)), 0)

    sig2 <- wrss / (n_eff - edf)
    tau2 <- sum(lk * qq^2 / d^2) / max(edf - penalty_order, 1e-6)

    if (tau2 < 1e-7) tau2 <- 1e-7

    lambda_old <- lambda
    lambda     <- sig2 / tau2
    lambda     <- max(1e-7, min(lambda, 1e7))

    if (abs(lambda - lambda_old) < tol) break
  }

  d    <- 1 + lambda * lk
  edf  <- sum(1 / d)
  beta <- backsolve(R, as.numeric(U %*% (qq / d)))

  list(
    fitted.values = as.numeric(X %*% beta),
    beta          = beta,
    edf           = edf,
    nl.df         = edf - penalty_order,
    lambda        = lambda,
    sig2          = sig2,
    iterations    = iter
  )
}


# -------------------------------------------------------------------
# 5c. RS loop identical to .fit_pb_sigma_rs_default but uses
#     .pb_fit_fast for the mu step.
# -------------------------------------------------------------------

.fit_pb_sigma_rs_fast <- function(y, x, X_sigma, weights,
                                   max_iter = 50,
                                   tol      = 1e-5,
                                   trace    = FALSE) {
  y      <- as.vector(y)
  x      <- as.vector(x)
  base_w <- as.vector(weights)

  stopifnot(length(y) == length(x))
  stopifnot(length(y) == length(base_w))
  stopifnot(nrow(X_sigma) == length(y))

  keep <- base_w > 0

  mu0 <- stats::weighted.mean(y[keep], base_w[keep])
  mu  <- rep(mu0, length(y))

  sig0 <- sqrt(stats::weighted.mean((y[keep] - mu0)^2, base_w[keep]))
  if (!is.finite(sig0) || sig0 <= 0) sig0 <- stats::sd(y[keep])
  if (!is.finite(sig0) || sig0 <= 0) sig0 <- 1

  eta_sigma  <- rep(log(sig0), length(y))
  beta_sigma <- rep(0, ncol(X_sigma))
  beta_sigma[1] <- log(sig0)

  dev_old <- Inf
  mu_fit  <- NULL

  for (iter in seq_len(max_iter)) {
    sigma <- exp(eta_sigma)
    w_mu  <- base_w / sigma^2

    mu_fit <- .pb_fit_fast(y = y, x = x, w = w_mu)   # <-- fast version
    mu     <- mu_fit$fitted.values

    sigma            <- exp(eta_sigma)
    resid2_scaled    <- ((y - mu) / sigma)^2
    z_sigma          <- eta_sigma + (resid2_scaled - 1) / 2
    w_sigma          <- base_w * 2

    lm_sigma  <- stats::lm.wfit(x = X_sigma, y = z_sigma, w = w_sigma)
    beta_new  <- as.vector(lm_sigma$coefficients)
    if (anyNA(beta_new)) beta_new[is.na(beta_new)] <- beta_sigma[is.na(beta_new)]
    beta_sigma    <- beta_new
    eta_sigma_new <- pmin(pmax(as.vector(X_sigma %*% beta_sigma), -30), 30)
    sigma_new     <- exp(eta_sigma_new)

    loglik  <- sum(base_w * stats::dnorm(y, mu, sigma_new, log = TRUE))
    dev     <- -2 * loglik

    if (trace)
      message("iter=", iter, " dev=", signif(dev, 8),
              " lambda=", signif(mu_fit$lambda, 6),
              " edf=", signif(mu_fit$edf, 5))

    if (is.finite(dev_old) &&
        abs(dev_old - dev) / (abs(dev_old) + 0.1) < tol) {
      eta_sigma <- eta_sigma_new; dev_old <- dev; break
    }

    eta_sigma <- eta_sigma_new
    dev_old   <- dev
  }

  sigma  <- exp(eta_sigma)
  logLik <- sum(base_w * stats::dnorm(y, mu, sigma, log = TRUE))

  list(
    mu         = mu,
    sigma      = sigma,
    beta_mu    = mu_fit$beta,
    beta_sigma = beta_sigma,
    lambda     = mu_fit$lambda,
    edf_mu     = mu_fit$edf,
    logLik     = logLik,
    deviance   = -2 * logLik,
    iterations = iter,
    converged  = iter < max_iter
  )
}


# -------------------------------------------------------------------
# 6. Plain-R RS loop for:
#      y ~ linear/poly, sigma ~ abs(x)
# -------------------------------------------------------------------

.fit_parametric_sigma_rs_default <- function(y, X_mu, X_sigma, weights,
                                             max_iter = 50,
                                             tol = 1e-5,
                                             trace = FALSE) {
  y <- as.vector(y)
  base_w <- as.vector(weights)

  stopifnot(nrow(X_mu) == length(y))
  stopifnot(nrow(X_sigma) == length(y))
  stopifnot(length(base_w) == length(y))

  keep <- base_w > 0

  beta_mu <- as.vector(stats::lm.wfit(
    x = X_mu[keep, , drop = FALSE],
    y = y[keep],
    w = base_w[keep]
  )$coefficients)

  beta_mu[is.na(beta_mu)] <- 0

  mu <- as.vector(X_mu %*% beta_mu)

  sig0 <- sqrt(stats::weighted.mean((y[keep] - mu[keep])^2, base_w[keep]))

  if (!is.finite(sig0) || sig0 <= 0) {
    sig0 <- stats::sd(y[keep] - mu[keep])
  }

  if (!is.finite(sig0) || sig0 <= 0) {
    sig0 <- 1
  }

  beta_sigma <- rep(0, ncol(X_sigma))
  beta_sigma[1] <- log(sig0)
  eta_sigma <- as.vector(X_sigma %*% beta_sigma)

  dev_old <- Inf

  for (iter in seq_len(max_iter)) {
    sigma <- exp(eta_sigma)

    # Mu update
    w_mu <- base_w / sigma^2

    mu_lm <- stats::lm.wfit(
      x = X_mu,
      y = y,
      w = w_mu
    )

    beta_new <- as.vector(mu_lm$coefficients)
    beta_new[is.na(beta_new)] <- beta_mu[is.na(beta_new)]
    beta_mu <- beta_new

    mu <- as.vector(X_mu %*% beta_mu)

    # Sigma update
    sigma <- exp(eta_sigma)
    resid2_scaled <- ((y - mu)^2) / sigma^2

    z_sigma <- eta_sigma + (resid2_scaled - 1) / 2
    w_sigma <- base_w * 2

    sigma_lm <- stats::lm.wfit(
      x = X_sigma,
      y = z_sigma,
      w = w_sigma
    )

    beta_sigma_new <- as.vector(sigma_lm$coefficients)
    beta_sigma_new[is.na(beta_sigma_new)] <-
      beta_sigma[is.na(beta_sigma_new)]

    beta_sigma <- beta_sigma_new
    eta_sigma_new <- as.vector(X_sigma %*% beta_sigma)

    eta_sigma_new <- pmin(pmax(eta_sigma_new, -30), 30)
    sigma_new <- exp(eta_sigma_new)

    loglik <- sum(
      base_w * stats::dnorm(y, mean = mu, sd = sigma_new, log = TRUE)
    )

    dev <- -2 * loglik

    if (trace) {
      message(
        "iter = ", iter,
        " dev = ", signif(dev, 8),
        " beta_mu = ",
        paste(signif(beta_mu, 6), collapse = ", "),
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
    base_w * stats::dnorm(y, mean = mu, sd = sigma, log = TRUE)
  )

  list(
    mu = mu,
    sigma = sigma,
    beta_mu = beta_mu,
    beta_sigma = beta_sigma,
    logLik = logLik,
    deviance = -2 * logLik,
    iterations = iter,
    converged = iter < max_iter
  )
}


# -------------------------------------------------------------------
# 7. Main drop-in gamlss replacement
# -------------------------------------------------------------------


# -------------------------------------------------------------------
# Minimal gamlss-like wrapper for remove_cardinal_biases()
# -------------------------------------------------------------------

gamlss_default_fun <- function(formula,
                               sigma.formula = NULL,
                               data,
                               weights = NULL,
                               control = NULL,
                               engine = c("r", "tmb"),
                               DLL = "circhelp",
                               ...) {
  engine <- match.arg(engine)

  if (is.null(sigma.formula)) {
    stop("sigma.formula must be supplied.", call. = FALSE)
  }

  if (missing(data) || is.null(data)) {
    stop("data must be supplied.", call. = FALSE)
  }

  data <- as.data.frame(data)

  y_name <- .extract_response_name(formula)

  if (!y_name %in% names(data)) {
    stop("Response variable not found in data: ", y_name, call. = FALSE)
  }

  y <- data[[y_name]]

  if (is.null(weights)) {
    weights <- rep(1, length(y))
  } else {
    weights <- eval(substitute(weights), data, parent.frame())
    weights <- as.vector(weights)
  }

  if (length(weights) != length(y)) {
    stop("weights must have the same length as the response.", call. = FALSE)
  }

  mu_type <- .gamlss_default_mu_type(formula)

  X_sigma <- .make_sigma_matrix_default(
    sigma.formula = sigma.formula,
    data = data
  )

  if (mu_type == "pb") {
    x_name <- .extract_single_rhs_var(formula)

    if (!x_name %in% names(data)) {
      stop("RHS variable not found in data: ", x_name, call. = FALSE)
    }

    x <- data[[x_name]]

    # Keep the validated R implementation for pb().
    # This reproduces the default gamlss pb + sigma model closely.
    fit <- .fit_pb_sigma_rs_default(
      y = y,
      x = x,
      X_sigma = X_sigma,
      weights = weights,
      trace = FALSE
    )

    beta_mu <- fit$beta_mu

  } else if (mu_type == "poly") {
    degree <- .extract_poly_degree(
      formula = formula,
      data = data,
      env = parent.frame()
    )

    x_name <- .extract_single_rhs_var(formula)

    if (!x_name %in% names(data)) {
      stop("RHS variable not found in data: ", x_name, call. = FALSE)
    }

    X_mu <- stats::model.matrix(
      stats::as.formula(sprintf("~ poly(%s, %d)", x_name, degree)),
      data = data
    )

    if (engine == "tmb") {
      fit <- .fit_parametric_sigma_tmb(
        y = y,
        X_mu = X_mu,
        X_sigma = X_sigma,
        weights = weights,
        DLL = DLL
      )
    } else {
      fit <- .fit_parametric_sigma_rs_default(
        y = y,
        X_mu = X_mu,
        X_sigma = X_sigma,
        weights = weights,
        trace = FALSE
      )
    }

    beta_mu <- fit$beta_mu

  } else {
    X_mu <- stats::model.matrix(formula, data = data)

    if (engine == "tmb") {
      fit <- .fit_parametric_sigma_tmb(
        y = y,
        X_mu = X_mu,
        X_sigma = X_sigma,
        weights = weights,
        DLL = DLL
      )
    } else {
      fit <- .fit_parametric_sigma_rs_default(
        y = y,
        X_mu = X_mu,
        X_sigma = X_sigma,
        weights = weights,
        trace = FALSE
      )
    }

    beta_mu <- fit$beta_mu
  }

  out <- list(
    fitted.values = fit$mu,
    mu = fit$mu,
    sigma = fit$sigma,

    coefficients = beta_mu,
    sigma.coefficients = fit$beta_sigma,

    logLik = fit$logLik,
    deviance = fit$deviance,

    weights = weights,
    formula = formula,
    sigma.formula = sigma.formula,
    data = data,
    call = match.call(),

    mu_type = mu_type,
    engine = engine,
    iterations = fit$iterations,
    converged = fit$converged,

    fit = fit
  )

  class(out) <- "gamlss_default_fit"

  out
}

# -------------------------------------------------------------------
# Fast drop-in: same as gamlss_default_fun but uses .pb_fit_fast
# (eigendecomp EM) for the pb() mu step.
# -------------------------------------------------------------------

gamlss_fast_fun <- function(formula,
                             sigma.formula = NULL,
                             data,
                             weights = NULL,
                             control = NULL,
                             ...) {
  if (is.null(sigma.formula)) stop("sigma.formula must be supplied.", call. = FALSE)
  if (missing(data) || is.null(data)) stop("data must be supplied.", call. = FALSE)

  data   <- as.data.frame(data)
  y_name <- .extract_response_name(formula)
  if (!y_name %in% names(data)) stop("Response variable not found: ", y_name, call. = FALSE)

  y <- data[[y_name]]

  if (is.null(weights)) {
    weights <- rep(1, length(y))
  } else {
    weights <- as.vector(eval(substitute(weights), data, parent.frame()))
  }

  mu_type <- .gamlss_default_mu_type(formula)
  X_sigma <- .make_sigma_matrix_default(sigma.formula = sigma.formula, data = data)

  if (mu_type == "pb") {
    x_name <- .extract_single_rhs_var(formula)
    if (!x_name %in% names(data)) stop("RHS variable not found: ", x_name, call. = FALSE)

    fit <- .fit_pb_sigma_rs_fast(
      y       = y,
      x       = data[[x_name]],
      X_sigma = X_sigma,
      weights = weights,
      trace   = isTRUE(control$trace)
    )
    beta_mu <- fit$beta_mu

  } else if (mu_type == "poly") {
    degree <- .extract_poly_degree(formula = formula, data = data, env = parent.frame())
    x_name <- .extract_single_rhs_var(formula)
    if (!x_name %in% names(data)) stop("RHS variable not found: ", x_name, call. = FALSE)
    X_mu <- stats::model.matrix(
      stats::as.formula(sprintf("~ poly(%s, %d)", x_name, degree)), data = data)
    fit     <- .fit_parametric_sigma_rs_default(y = y, X_mu = X_mu,
                 X_sigma = X_sigma, weights = weights, trace = FALSE)
    beta_mu <- fit$beta_mu

  } else {
    X_mu <- stats::model.matrix(formula, data = data)
    fit     <- .fit_parametric_sigma_rs_default(y = y, X_mu = X_mu,
                 X_sigma = X_sigma, weights = weights, trace = FALSE)
    beta_mu <- fit$beta_mu
  }

  out <- list(
    fitted.values      = fit$mu,
    mu                 = fit$mu,
    sigma              = fit$sigma,
    coefficients       = beta_mu,
    sigma.coefficients = fit$beta_sigma,
    logLik             = fit$logLik,
    deviance           = fit$deviance,
    weights            = weights,
    formula            = formula,
    sigma.formula      = sigma.formula,
    data               = data,
    call               = match.call(),
    mu_type            = mu_type,
    engine             = "fast",
    iterations         = fit$iterations,
    converged          = fit$converged,
    fit                = fit
  )
  class(out) <- "gamlss_default_fit"
  out
}


# TMB::compile("src/normal_ls.cpp")
# dyn.load(TMB::dynlib("src/normal_ls"))
# -------------------------------------------------------------------
# TMB-backed convenience wrapper
# -------------------------------------------------------------------

gamlss_tmb_fun <- function(formula,
                           sigma.formula = NULL,
                           data,
                           weights = NULL,
                           control = NULL,
                           ...) {
  gamlss_default_fun(
    formula = formula,
    sigma.formula = sigma.formula,
    data = data,
    weights = weights,
    control = control,
    engine = "tmb",
    DLL = "normal_ls",
    ...
  )
}


.fit_normal_ls_tmb <- function(y,
                               X_mu,
                               X_sigma,
                               weights  = NULL,
                               start    = NULL,
                               D_mu     = matrix(0, 0, ncol(X_mu)),
                               lambda_mu = 0,
                               DLL      = "circhelp",
                               control  = list()) {
  if (!requireNamespace("TMB", quietly = TRUE)) {
    stop("Package 'TMB' is required for engine = 'tmb'.", call. = FALSE)
  }

  y <- as.vector(y)
  X_mu <- as.matrix(X_mu)
  X_sigma <- as.matrix(X_sigma)

  if (is.null(weights)) {
    weights <- rep(1, length(y))
  }

  weights <- as.vector(weights)

  stopifnot(length(y) == nrow(X_mu))
  stopifnot(length(y) == nrow(X_sigma))
  stopifnot(length(y) == length(weights))

  keep <- weights > 0

  if (!any(keep)) {
    stop("All weights are zero.", call. = FALSE)
  }

  if (is.null(start)) {
    beta_mu_start <- as.vector(stats::lm.wfit(
      x = X_mu[keep, , drop = FALSE],
      y = y[keep],
      w = weights[keep]
    )$coefficients)

    beta_mu_start[is.na(beta_mu_start)] <- 0

    mu_start <- as.vector(X_mu %*% beta_mu_start)

    sig0 <- sqrt(stats::weighted.mean(
      (y[keep] - mu_start[keep])^2,
      weights[keep]
    ))

    if (!is.finite(sig0) || sig0 <= 0) {
      sig0 <- stats::sd(y[keep] - mu_start[keep])
    }

    if (!is.finite(sig0) || sig0 <= 0) {
      sig0 <- 1
    }

    beta_sigma_start <- rep(0, ncol(X_sigma))
    beta_sigma_start[1] <- log(sig0)

    start <- list(
      beta_mu = beta_mu_start,
      beta_sigma = beta_sigma_start
    )
  }

  obj <- TMB::MakeADFun(
    data = list(
      y        = y,
      X_mu     = X_mu,
      X_sigma  = X_sigma,
      weights  = weights,
      D_mu     = D_mu,
      lambda_mu = lambda_mu
    ),
    parameters = start,
    DLL = DLL,
    silent = TRUE
  )

  opt <- stats::nlminb(
    start = obj$par,
    objective = obj$fn,
    gradient = obj$gr,
    control = utils::modifyList(
      list(
        eval.max = 1000,
        iter.max = 1000,
        rel.tol = 1e-10,
        x.tol = 1e-10
      ),
      control
    )
  )

  par <- obj$env$parList(opt$par)

  mu <- as.vector(X_mu %*% par$beta_mu)
  sigma <- exp(as.vector(X_sigma %*% par$beta_sigma))

  logLik <- sum(
    weights * stats::dnorm(y, mean = mu, sd = sigma, log = TRUE)
  )

  list(
    mu = mu,
    sigma = sigma,
    beta_mu = par$beta_mu,
    beta_sigma = par$beta_sigma,
    logLik = logLik,
    deviance = -2 * logLik,
    iterations = opt$iterations,
    converged = opt$convergence == 0,
    opt = opt
  )
}


.fit_parametric_sigma_tmb <- function(y,
                                      X_mu,
                                      X_sigma,
                                      weights,
                                      DLL = "circhelp") {
  .fit_normal_ls_tmb(
    y = y,
    X_mu = X_mu,
    X_sigma = X_sigma,
    weights = weights,
    DLL = DLL
  )
}
# -------------------------------------------------------------------
# 8. Methods needed by remove_cardinal_biases()
# -------------------------------------------------------------------

predict.gamlss_default_fit <- function(object,
                                       newdata = NULL,
                                       what = c("mu", "sigma"),
                                       type = c("response", "link"),
                                       ...) {
  what <- match.arg(what)
  type <- match.arg(type)

  if (!is.null(newdata)) {
    stop("newdata prediction is not implemented for gamlss_default_fit.")
  }

  if (what == "mu") {
    return(as.vector(object$mu))
  }

  if (what == "sigma") {
    if (type == "link") {
      return(log(as.vector(object$sigma)))
    }

    return(as.vector(object$sigma))
  }
}


coef.gamlss_default_fit <- function(object,
                                    what = c("mu", "sigma"),
                                    ...) {
  what <- match.arg(what)

  if (what == "mu") {
    return(object$coefficients)
  }

  object$sigma.coefficients
}


logLik.gamlss_default_fit <- function(object, ...) {
  val <- object$logLik

  attr(val, "df") <- length(object$coefficients) +
    length(object$sigma.coefficients)

  attr(val, "nobs") <- sum(object$weights > 0)

  class(val) <- "logLik"

  val
}


.find_call <- function(expr, fun_name) {
  if (!is.call(expr)) {
    return(NULL)
  }

  # Handles poly(...), stats::poly(...), gamlss::pb(...), etc.
  call_head <- expr[[1]]

  is_target <- FALSE

  if (is.symbol(call_head)) {
    is_target <- identical(as.character(call_head), fun_name)
  } else if (is.call(call_head) && identical(as.character(call_head[[1]]), "::")) {
    is_target <- identical(as.character(call_head[[3]]), fun_name)
  }

  if (is_target) {
    return(expr)
  }

  for (i in seq_along(expr)[-1]) {
    out <- .find_call(expr[[i]], fun_name)
    if (!is.null(out)) {
      return(out)
    }
  }

  NULL
}


# ===================================================================
# REML P-spline alternative for gamlss_fun
# ===================================================================

# -------------------------------------------------------------------
# R1. B-spline basis + 2nd-difference penalty (matches gamlss::pb())
# -------------------------------------------------------------------

.make_pspline_reml <- function(x, ndx = NULL, deg = 3L, ord = 2L) {
  x <- as.vector(x)
  n_dist <- length(unique(x))
  if (is.null(ndx)) ndx <- if (length(x) < 100L) 10L else 20L
  ndx <- min(ndx, n_dist)
  xl <- min(x); xr <- max(x)
  rng <- xr - xl
  xl  <- xl - 0.01 * rng
  xr  <- xr + 0.01 * rng
  dx  <- (xr - xl) / ndx
  knots <- seq(xl - deg * dx, xr + deg * dx, by = dx)
  B <- splines::splineDesign(knots, x, ord = deg + 1L, outer.ok = TRUE)
  K <- ncol(B)
  D <- diff(diag(K), differences = ord)
  S <- crossprod(D)
  list(X = B, S = S)
}

# -------------------------------------------------------------------
# R2. REML criterion for lambda selection
#     criterion: WRSS(lam) + sum(log(1 + lam*lk)) - rk*log(lam)
# -------------------------------------------------------------------

.gcv_pspline_reml <- function(y_w, B, S_full,
                               lambda_range = c(-15, 15)) {
  p    <- ncol(B)
  qr_B <- qr(B)
  Q    <- qr.Q(qr_B)
  R    <- qr.R(qr_B)
  Rinv <- backsolve(R, diag(p))
  K    <- crossprod(Rinv, S_full %*% Rinv)
  K    <- (K + t(K)) / 2
  eig  <- eigen(K, symmetric = TRUE)
  U    <- eig$vectors
  lk   <- pmax(eig$values, 0)
  rk   <- sum(lk > .Machine$double.eps^0.5)
  z_w  <- as.numeric(crossprod(Q, y_w))
  qq   <- as.numeric(crossprod(U, z_w))
  y_y  <- sum(y_w^2)

  reml_fn <- function(log_lam) {
    lam  <- exp(log_lam)
    d    <- 1 + lam * lk
    wrss <- max(y_y - sum(qq^2 * (2 / d - 1 / d^2)), 0)
    wrss + sum(log(d)) - rk * log_lam
  }

  opt  <- stats::optimize(reml_fn, interval = lambda_range)
  lam  <- exp(opt$minimum)
  d    <- 1 + lam * lk
  beta <- backsolve(R, as.numeric(U %*% (qq / d)))
  edf  <- sum(1 / d)
  list(lambda = lam, edf = edf, beta = beta)
}

# -------------------------------------------------------------------
# R3. RS loop: REML P-spline mu + nlminb sigma
# -------------------------------------------------------------------

.fit_pb_sigma_reml <- function(y, x, X_sigma, weights,
                                max_iter = 20,
                                tol = 1e-5,
                                trace = FALSE) {
  y      <- as.vector(y)
  x      <- as.vector(x)
  base_w <- as.vector(weights)

  stopifnot(length(y) == length(x))
  stopifnot(length(y) == length(base_w))
  stopifnot(nrow(X_sigma) == length(y))

  # Identify structural zeros BEFORE any weight flooring
  pos <- base_w > 0

  sp   <- .make_pspline_reml(x)
  B    <- sp$X
  S    <- sp$S

  # Initialise sigma parameters
  keep  <- pos
  mu0   <- stats::weighted.mean(y[keep], base_w[keep])
  sig0  <- sqrt(stats::weighted.mean((y[keep] - mu0)^2, base_w[keep]))
  if (!is.finite(sig0) || sig0 <= 0) sig0 <- stats::sd(y[keep])
  if (!is.finite(sig0) || sig0 <= 0) sig0 <- 1

  gamma <- rep(0, ncol(X_sigma))
  gamma[1] <- log(sig0)

  mu       <- rep(mu0, length(y))
  beta_mu  <- NULL
  dev_old  <- Inf

  for (iter in seq_len(max_iter)) {
    sigma_vec <- as.numeric(exp(X_sigma %*% gamma))
    w_mu      <- base_w / sigma_vec^2

    # Weighted-sqrt transform for REML (only pos rows)
    sqw       <- sqrt(w_mu[pos])
    B_gcv     <- B[pos, , drop = FALSE] * sqw
    y_gcv     <- y[pos] * sqw

    gcv_res <- tryCatch(
      .gcv_pspline_reml(y_gcv, B_gcv, S),
      error = function(e) NULL
    )

    if (is.null(gcv_res) || !all(is.finite(gcv_res$beta))) break

    beta_mu <- gcv_res$beta
    mu      <- as.numeric(B %*% beta_mu)
    resid   <- y - mu

    sigma_nll <- function(par) {
      sig <- exp(as.numeric(X_sigma %*% par))
      if (any(sig <= 0)) return(.Machine$double.xmax)
      sum(base_w * (log(sig) + 0.5 * (resid / sig)^2))
    }

    opt_g <- try(
      stats::nlminb(gamma, sigma_nll,
                    control = list(iter.max = 100, eval.max = 200)),
      silent = TRUE
    )

    if (!inherits(opt_g, "try-error") && all(is.finite(opt_g$par))) {
      gamma_new <- opt_g$par
    } else {
      gamma_new <- gamma
    }

    sigma_new <- as.numeric(exp(X_sigma %*% gamma_new))
    loglik    <- sum(base_w * stats::dnorm(y, mu, sigma_new, log = TRUE))
    dev       <- -2 * loglik

    if (trace) {
      message("iter = ", iter,
              " dev = ", signif(dev, 8),
              " lambda = ", signif(gcv_res$lambda, 6),
              " edf = ", signif(gcv_res$edf, 5),
              " gamma = ", paste(signif(gamma_new, 5), collapse = ", "))
    }

    if (is.finite(dev_old)) {
      if (abs(dev_old - dev) / (abs(dev_old) + 0.1) < tol) {
        gamma   <- gamma_new
        dev_old <- dev
        break
      }
    }

    gamma   <- gamma_new
    dev_old <- dev
  }

  sigma  <- as.numeric(exp(X_sigma %*% gamma))
  logLik <- sum(base_w * stats::dnorm(y, mu, sigma, log = TRUE))

  list(
    mu         = mu,
    sigma      = sigma,
    beta_mu    = beta_mu,
    beta_sigma = gamma,
    lambda     = if (!is.null(gcv_res)) gcv_res$lambda else NA_real_,
    edf_mu     = if (!is.null(gcv_res)) gcv_res$edf    else NA_real_,
    logLik     = logLik,
    deviance   = -2 * logLik,
    iterations = iter,
    converged  = iter < max_iter
  )
}

# -------------------------------------------------------------------
# R4. Public gamlss_fun interface
# -------------------------------------------------------------------

gamlss_reml_fun <- function(formula,
                             sigma.formula = NULL,
                             data,
                             weights = NULL,
                             control = NULL,
                             ...) {
  if (is.null(sigma.formula)) {
    stop("sigma.formula must be supplied.", call. = FALSE)
  }

  if (missing(data) || is.null(data)) {
    stop("data must be supplied.", call. = FALSE)
  }

  data   <- as.data.frame(data)
  y_name <- .extract_response_name(formula)

  if (!y_name %in% names(data)) {
    stop("Response variable not found in data: ", y_name, call. = FALSE)
  }

  y <- data[[y_name]]

  if (is.null(weights)) {
    weights <- rep(1, length(y))
  } else {
    weights <- eval(substitute(weights), data, parent.frame())
    weights <- as.vector(weights)
  }

  if (length(weights) != length(y)) {
    stop("weights must have the same length as the response.", call. = FALSE)
  }

  mu_type <- .gamlss_default_mu_type(formula)

  X_sigma <- .make_sigma_matrix_default(
    sigma.formula = sigma.formula,
    data          = data
  )

  if (mu_type == "pb") {
    x_name <- .extract_single_rhs_var(formula)

    if (!x_name %in% names(data)) {
      stop("RHS variable not found in data: ", x_name, call. = FALSE)
    }

    x <- data[[x_name]]

    fit <- .fit_pb_sigma_reml(
      y        = y,
      x        = x,
      X_sigma  = X_sigma,
      weights  = weights,
      trace    = isTRUE(control$trace)
    )

    beta_mu <- fit$beta_mu

  } else if (mu_type == "poly") {
    degree <- .extract_poly_degree(
      formula = formula,
      data    = data,
      env     = parent.frame()
    )

    x_name <- .extract_single_rhs_var(formula)

    if (!x_name %in% names(data)) {
      stop("RHS variable not found in data: ", x_name, call. = FALSE)
    }

    X_mu <- stats::model.matrix(
      stats::as.formula(sprintf("~ poly(%s, %d)", x_name, degree)),
      data = data
    )

    fit     <- .fit_parametric_sigma_rs_default(
      y       = y,
      X_mu    = X_mu,
      X_sigma = X_sigma,
      weights = weights,
      trace   = isTRUE(control$trace)
    )
    beta_mu <- fit$beta_mu

  } else {
    X_mu <- stats::model.matrix(formula, data = data)

    fit     <- .fit_parametric_sigma_rs_default(
      y       = y,
      X_mu    = X_mu,
      X_sigma = X_sigma,
      weights = weights,
      trace   = isTRUE(control$trace)
    )
    beta_mu <- fit$beta_mu
  }

  out <- list(
    fitted.values      = fit$mu,
    mu                 = fit$mu,
    sigma              = fit$sigma,
    coefficients       = beta_mu,
    sigma.coefficients = fit$beta_sigma,
    logLik             = fit$logLik,
    deviance           = fit$deviance,
    weights            = weights,
    formula            = formula,
    sigma.formula      = sigma.formula,
    data               = data,
    call               = match.call(),
    mu_type            = mu_type,
    engine             = "reml",
    iterations         = fit$iterations,
    converged          = fit$converged,
    fit                = fit
  )

  class(out) <- "gamlss_default_fit"

  out
}