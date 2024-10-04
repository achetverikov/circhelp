#' Circular mean
#'
#' @param x vector of values
#' @param na.rm a logical value indicating whether NA values should be removed before the computation proceeds
#'
#' @return mean of values in the vector
#' @export
#'
#' @examples
#' x <- runif(1000, -pi, pi)
#' mean(x)
#' circ_mean_rad(x)
#'
#' @describeIn circ_mean_rad circular mean in 2pi space
circ_mean_rad <- function(x, na.rm = FALSE) {
  sinr <- sum(sin(x), na.rm = na.rm)
  cosr <- sum(cos(x), na.rm = na.rm)
  circmean <- atan2(sinr, cosr)
  circmean
}

#' @describeIn circ_mean_rad circular mean in 180° space (e.g., line orientation)
#' @export
circ_mean_180 <- function(x, na.rm = FALSE) {
  circ_mean_rad(x / 90 * pi, na.rm = na.rm) / pi * 90
}

#' @describeIn circ_mean_rad circular mean in 360° space
#' @export
circ_mean_360 <- function(x, na.rm = FALSE) {
  circ_mean_rad(x / 180 * pi, na.rm = na.rm) / pi * 180
}

#' Differences between angles in different circular spaces
#'
#' @param a first angle
#' @param b second angle
#' @details By default, all functions return values in ± half-range space (e.g., -pi to pi for 2pi radian space used by `angle_diff_rad()`) but `angle_diff_180_45()` and `angle_diff_360_90()` return values in \[-1/4 range, 3/4 range\] space
#'
#' @return difference between a and b
#' @export
#'
#' @examples
#' angle_diff_180(5, 175)
#' angle_diff_360(5, 175)
#' angle_diff_90(5, 175)
#' angle_diff_rad(5, 175)
#'
#' angle_diff_360(300, 0)
#' angle_diff_360_90(300, 0)

#' @describeIn angle_diff_rad angle difference in radians
#' @export
angle_diff_rad <- function(a, b) {
  c <- a - b
  (c + pi) %% (2 * pi) - pi
}

#' @describeIn angle_diff_rad angle difference in 360 degree space
#' @export

angle_diff_360 <- function(a, b) {
  angle_diff_rad(a / 180 * pi, b / 180 * pi) / pi * 180
}

#' @describeIn angle_diff_rad angle difference in 180 degree space (e.g., line orientation)
#' @export
#'
angle_diff_180 <- function(a, b) {
  angle_diff_rad(a / 90 * pi, b / 90 * pi) / pi * 90
}

#' @describeIn angle_diff_rad angle difference in 90 degree space
#' @export

angle_diff_90 <- function(a, b) {
  angle_diff_rad(a / 45 * pi, b / 45 * pi) / pi * 45
}

#' @describeIn angle_diff_rad angle difference in 180 degree space from -45 to 135
#' @export

angle_diff_180_45 <- function(a, b) {
  c <- a - b
  (c + 45) %% 180 - 45
}

#' @describeIn angle_diff_rad angle difference in 360 degree space from -90 to 270
#' @export

angle_diff_360_90 <- function(a, b) {
  c <- a - b
  (c + 90) %% 360 - 90
}

#' Circular correlation coefficient
#'
#' Computes a circular correlation coefficient as defined in Jammalamadaka & SenGupta (2001).
#' @param a first variable
#' @param b second variable
#' @param ill_defined is one of the variables mean is not well-defined (e.g., it is uniformly distributed)?
#' @param mu fix the mean parameter of both vectors to a certain value
#' @param na.rm a logical value indicating whether NA values should be removed before the computation proceeds
#'
#' @return correlation coefficient
#' @references {
#' Jammalamadaka, S. R., & SenGupta, A. (2001). Topics in Circular Statistics. WORLD SCIENTIFIC. \doi{10.1142/4031}
#' }
#' @export
#'
#' @examples
#' requireNamespace("mgcv")
#' data <- mgcv::rmvn(10000, c(0, 0), V = matrix(c(1, 0.5, 0.5, 1), ncol = 2))
#' circ_corr(data[, 1], data[, 2])
circ_corr <- function(a, b, ill_defined = FALSE, mu = NULL, na.rm = FALSE) {
  if (na.rm) {
    a <- a[!is.na(a)]
    b <- b[!is.na(b)]
  }
  mu_a <- circ_mean_rad(a)
  mu_b <- circ_mean_rad(b)
  if (!is.null(mu)) {
    mu_a <- mu_b <- mu
  } else if (ill_defined) {
    mean_diff <- circ_mean_rad(a - b)
    mean_sum <- circ_mean_rad(a + b)
    mu_a <- (mean_diff + mean_sum) / 2
    mu_b <- (mean_sum - mean_diff) / 2
  }
  sin_a <- sin(a - mu_a)
  sin_b <- sin(b - mu_b)
  rho <- sum(sin_a * sin_b) / sqrt(sum(sin_a * sin_a) * sum(sin_b * sin_b))
  rho
}

#' Circular-linear correlation
#'
#' \loadmathjax
#' Implementation of the circular-linear correlation measure introduced by Mardia (1976) and Johnson and Wehrly (1977) as cited in Jammalamadaka & Sengupta (2001).
#'
#' @param circ_x circular variable
#' @param lin_x linear variable
#' @param na.rm a logical value indicating whether NA values should be removed before the computation proceeds
#'
#' @details This measure is computed as \mjsdeqn{r^2 = (r_{xc}^2+r_{xs}^2-2 r_{xc} r_{xs}r_{cs})/(1-r_{cs}^2)} where \mjseqn{r_{xc} = corr(x, cos(\alpha))}, \mjseqn{r_{xs} = corr(x, sin(\alpha))}, \mjseqn{r_{cs} = corr(cos(\alpha), sin(\alpha))}, and \mjseqn{\alpha} and \mjseqn{x} are the circular and linear variables, respectively.
#'
#'
#' @return circular-linear correlation measure
#' @export
#' @references {
#' Jammalamadaka, S. R., & SenGupta, A. (2001). Topics in Circular Statistics. WORLD SCIENTIFIC. \doi{10.1142/4031}
#' }
#' @examples
#'
#' x <- rnorm(50)
#' a <- as.vector(circular::rvonmises(50, 0, 5))
#' circ_lin_corr(x + a, x)
circ_lin_corr <- function(circ_x, lin_x, na.rm = FALSE) {
  if (na.rm) {
    circ_x <- circ_x[!is.na(circ_x)]
    lin_x <- lin_x[!is.na(lin_x)]
  }
  cos_a <- cos(circ_x)
  sin_a <- sin(circ_x)
  r_xcos <- stats::cor(lin_x, cos_a)
  r_xsin <- stats::cor(lin_x, sin_a)
  r_cossin <- stats::cor(cos_a, sin_a)

  r_squared <- ((r_xcos^2) + (r_xsin^2) - (2 * r_xcos * r_xsin * r_cossin)) / (1 - (r_cossin^2))

  sqrt(r_squared)
}
#' Weighted circular parameters
#'
#' @param x vector of values (in radians)
#' @param w vector of weights
#' @param na.rm a logical value indicating whether NA values should be removed before the computation proceeds
#'
#' @return weighted mean of values in the vector
#' @export
#'
#' @examples
#' x <- rnorm(1000, 0, 0.5)
#' w <- runif(1000, 0, 1)
#' weighted.mean(x, w)
#' weighted_circ_mean(x, w)
#'
#' @describeIn weighted_circ_mean weighted circular mean

weighted_circ_mean <- function(x, w, na.rm = FALSE) {
  if (length(w) != length(x)) {
    stop("Weights (w) should have the same length as values (x)")
  }

  sum_w <- sum(w, na.rm = na.rm)
  atan2(sum(w * sin(x), na.rm = na.rm) / sum_w, sum(w * cos(x), na.rm = na.rm) / sum_w)
}

#' @describeIn weighted_circ_mean an alternative way to compute weighted circular mean (the results are the same)
#' @export
weighted_circ_mean2 <- function(x, w, na.rm = FALSE) {
  if (length(w) != length(x)) {
    stop("Weights (w) should have the same length as values (x)")
  }
  z <- exp(1i * x)
  Arg(sum(w * z, na.rm = na.rm) / sum(w, na.rm = na.rm))
}

#' @describeIn weighted_circ_mean weighted circular SD
#' @export

weighted_circ_sd <- function(x, w, na.rm = FALSE) {
  sum_w <- sum(w, na.rm = na.rm)

  r <- sqrt((sum(w * sin(x), na.rm = na.rm) / sum_w)^2 + (sum(w * cos(x), na.rm = na.rm) / sum_w)^2)
  sqrt(-2 * log(r))
}

#' @describeIn weighted_circ_mean weighted mean resultant length
#' @export

weighted_circ_rho <- function(x, w, na.rm = FALSE) {
  sum_w <- sum(w, na.rm = na.rm)

  r <- sqrt((sum(w * sin(x), na.rm = na.rm) / sum_w)^2 + (sum(w * cos(x), na.rm = na.rm) / sum_w)^2)
  r
}

#' Get angle value in \[-pi, pi\] space
#'
#' @param x angle
#'
#' @return angle in \[-pi, pi\] space
#' @export
#'
#' @examples
#' correct_angle_rad(4 * pi)
#'
correct_angle_rad <- function(x) {
  angle_diff_rad(x, 0)
}

#' Circular standard deviation
#'
#' @param x vector of angles
#' @param na.rm a logical value indicating whether NA values should be removed before the computation proceeds
#'
#' @return standard deviation of values in the vector
#' @export
#'
#' @examples
#' circ_sd_rad(rnorm(50))
#' circ_sd_180(rnorm(50))
#'
#' @describeIn circ_sd_rad SD of angles in radians
circ_sd_rad <- function(x, na.rm = FALSE) {
  if (na.rm) {
    x <- x[!is.na(x)]
  }

  r <- sum(exp(1i * x))

  # mean resultant length
  rho <- abs(r) / length(x)

  sqrt(-2 * log(rho))
}

#' @describeIn circ_sd_rad SD of angles in 360 degree space
#' @export

circ_sd_360 <- function(x, na.rm = FALSE) {
  circ_sd_rad(x / 180 * pi, na.rm = na.rm) / pi * 180
}
#' @describeIn circ_sd_rad SD of angles in 180 degree space
#' @export
circ_sd_180 <- function(x, na.rm = FALSE) {
  circ_sd_rad(x / 90 * pi, na.rm = na.rm) / pi * 90
}

circ_dist <- function(x, y) {
  angle(exp(1i * x) / exp(1i * y))
}

angle <- function(x) {
  atan2(Im(x), Re(x))
}

#' A set of descriptive statistics for circular data
#'
#' @param x vector of angles
#' @param w weights for the values in the vector
#' @param d correction for the bias for data with known spacing
#' @param na.rm a logical value indicating whether NA values should be removed before the computation proceeds
#'
#' @return a list with descriptive statistics
#' \itemize{
#'  \item mu - mean
#'  \item sigma - standard deviation
#'  \item skew_pewsey - skewness as defined by Pewsey
#'  \item skew_fischer - skewness as defined by Fischer
#'  \item rho - mean resultant length
#'  \item skew_rel_to_zero - skewness relative to zero
#' }
#' @export
#'
#' @examples
#' x <- c(rnorm(50, 0, 0.5), rnorm(20, 1, 0.5))
#' circ_descr(x)
#'
circ_descr <- function(x, w = NULL, d = NULL, na.rm = FALSE) {
  if (na.rm) {
    x <- x[!is.na(x)]
  }
  if (is.null(w)) {
    w <- rep(1, length(x))
  }
  # compute weighted sum of cos and sin of angles
  r <- sum(w * exp(1i * x))
  mu <- angle(r)

  # mean resultant length
  rho <- abs(r) / sum(w)
  # for data with known spacing, apply correction factor to correct for bias
  # in the estimation of r (see Zar, p. 601, equ. 26.16)
  if (!is.null(d)) {
    c <- d / 2 / sin(d / 2)
    rho <- c * rho
  }

  # 2nd moment
  alpha_rel <- angle(exp(1i * x) / exp(1i * mu))

  p <- 2
  cbar <- mean(cos(p * alpha_rel) * w)
  sbar <- mean(sin(p * alpha_rel) * w)
  mp <- cbar + 1i * sbar

  rho2 <- abs(mp) # shouldn't rho2 be adjusted for d as well?
  mu2 <- angle(mp)

  b <- sum(w * sin(2 * circ_dist(x, mu))) / sum(w)
  b0 <- rho2 * sin(circ_dist(mu2, 2 * mu)) / (1 - rho)^(3 / 2)

  b_rel_to_zero <- sum(w * sin(2 * circ_dist(x, 0))) / sum(w)
  s0 <- sqrt(-2 * log(rho))
  list(mu = mu, sigma = s0, skew_pewsey = b, skew_fischer = b0, rho = rho, skew_rel_to_zero = b_rel_to_zero)
}

#' Remove cardinal biases
#'
#' @param err a vector of errors, deviations of response from the true stimuli
#' @param x a vector of true stimuli in degrees (see space)
#' @param space circular space to use (a string: `180` or `360`)
#' @param bias_type bias type to use (`fit`, `card`, `obl`, or `custom`, see details)
#' @param plots a string `hide`, `show`, or `return` to hide, show, or return plots (default: `hide`)
#' @param do_plots deprecated, use the parameter `plots` instead
#' @param poly_deg degree of the fitted polynomials for each bin (default: 4)
#' @param var_sigma allow standard deviation (width) of the fitted response distribution to vary as a function of distance to the nearest cardinal (default: True)
#' @param var_sigma_poly_deg degree of the fitted polynomials for each bin for the first approximation for the response distribution to select the best fitting model (default: 4)
#' @param reassign_at_boundaries select the bin for the observations at the boundaries between bins based on the best-fitting polynomial (default: True)
#' @param reassign_range maximum distance to the boundary at which reassignment can occur (default: 2 degrees)
#' @param break_points can be used to assign custom break points instead of cardinal/oblique ones with `bias_type` set to `custom` (default: NULL)
#' @param init_outliers a vector determining which errors are initially assumed to be outliers (default: NULL)
#' @param debug print some extra info (default: False)
#'
#' @details
#' If the `bias_type` is set to `fit`, the function computes the cardinal biases in the following way:
#' \enumerate{
#' \item Create two sets of bins, splitting the stimuli vector into bins centered at cardinal and at oblique directions.
#' \item For each set of bins, fit a nth-degree polynomial for the responses in each bin, optionally allowing the distribution of responses to vary in width as a function of distance to the nearest cardinal (regardless of whether the bins are centered at the cardinal or at the oblique, the width of the response distribution usually increases as the distance to cardinals increase).
#' \item Choose the best-fitting model between the one using cardinal and the one using oblique bins.
#' \item Compute the residuals of the best-fitting model - that's your bias-corrected error - and the biases (see below).
#' }
#' The bias is computed by flipping the sign of errors when the average predicted error is negative, so, that, for example, if on average the responses are shifted clockwise relative to the true values, the trial-by-trial error would count as bias when it is also shifted clockwise.
#'
#' If `bias_type` is set to `obl` or `card`, only one set of bins is used, centred at cardinal or oblique angles, respectively.
#'
#' For additional examples see the help vignette:
#' \code{vignette("cardinal_biases", package = "circhelp")}
#'
#' @return If `plots=='return'`, returns the three plots showing the biases
#' (combined together with [patchwork::wrap_plots()]). Otherwise, returns a list with the following elements:
#' \itemize{
#' \item is_outlier - 0 for outliers (defined as `±3*pred_sigma` for the model with varying sigma or as `±3\*SD` for the simple model)
#' \item pred predicted error
#' \item be_c error corrected for biases (`be_c = observed error - pred`)
#' \item which_bin the numeric ID of the bin that the stimulus belong to
#' \item bias the bias computed as described above
#' \item bias_typ bias type (cardinal or oblique)
#' \item pred_lin predicted error for a simple linear model for comparison
#' \item pred_sigma predicted SD of the error distribution
#' \item coef_sigma_int, coef_sigma_slope intercept and slope for the sigma prediction
#'
#' }
#' @references {
#' \itemize{
#' \item Chetverikov, A., & Jehee, J. F. M. (2023). Motion direction is represented as a bimodal probability distribution in the human visual cortex. Nature Communications, 14(7634). \doi{10.1038/s41467-023-43251-w}
#' \item van Bergen, R. S., Ma, W. J., Pratte, M. S., & Jehee, J. F. M. (2015). Sensory uncertainty decoded from visual cortex predicts behavior. Nature Neuroscience, 18(12), 1728-1730. \doi{10.1038/nn.4150}
#' }
#' }
#' @export
#' @import data.table gamlss
#' @importFrom MASS rlm
#' @examples
#'
#' # Data in orientation domain from Pascucci et al. (2019, PLOS Bio),
#' # https://doi.org/10.5281/zenodo.2544946
#'
#' ex_data <- Pascucci_et_al_2019_data[observer == 4, ]
#' remove_cardinal_biases(ex_data$err, ex_data$orientation, plots = "show")
#'
#' # Data in motion domain from Bae & Luck (2018, Neuroimage),
#' # https://osf.io/2h6w9/
#' ex_data_bae <- Bae_Luck_2018_data[subject_Num == unique(subject_Num)[5], ]
#' remove_cardinal_biases(ex_data_bae$err, ex_data_bae$TargetDirection,
#'   space = "360", plots = "show"
#' )
#'
#' # Using a stricter initial outlier boundary
#'
#' remove_cardinal_biases(ex_data_bae$err, ex_data_bae$TargetDirection,
#'   space = "360", plots = "show",
#'   init_outliers = abs(ex_data_bae$err) > 60
#' )
#'
#' # We can also use just one bin by setting `bias_type` to custom
#' # and setting the `break_points` at the ends of the range for x
#'
#' remove_cardinal_biases(ex_data_bae$err, ex_data_bae$TargetDirection,
#'   space = "360", bias_type = "custom",
#'   break_points = c(-180, 180), plots = "show",
#'   reassign_at_boundaries = FALSE, poly_deg = 8,
#'   init_outliers = abs(ex_data_bae$err) > 60
#' )
#'
remove_cardinal_biases <- function(err, x, space = "180", bias_type = "fit", plots = "hide", poly_deg = 4, var_sigma = TRUE, var_sigma_poly_deg = 4, reassign_at_boundaries = TRUE, reassign_range = 2, break_points = NULL, init_outliers = NULL, debug = FALSE, do_plots = NULL) {
  outlier <- dist_to_card <- dist_to_obl <- logLik <- x_var <- min_bp_i <- center_x <- dc_var <- gr_var <- min_boundary_i <- min_boundary_dist <- bin_range <- bin_boundary_left <- bin_boundary_right <- at_the_boundary <- row_i <- likelihood <- dnorm <- pred <- pred_sigma <- new_weight <- i.gr_var <- dist_to_bin_centre <- coef <- predict <- bias <- pred_lin <- be_c <- which_bin <- center_y <- outlier_f <- coef_sigma_int <- . <- coef_sigma_slope <- NULL # due to NSE notes in R CMD check

  if (!(bias_type %in% c("fit", "card", "obl", "custom"))) {
    stop("`bias_type` should be 'fit','card', 'obl', or 'custom'")
  }
  if (bias_type == "custom" & missing(break_points)) {
    stop("If 'bias_type' is set to 'custom', you need to specify 'break_points'")
  }

  if (any(is.na(x)) | any(is.na(err))) {
    stop("There are NAs in x or err. Please remove missing values before running the function.")
  }


  if (!missing(do_plots)) {
    warnings("\nYou have supplied 'do_plots' argument, it is now deprecated in favor of a 'plots' argument")
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
    obl_groups <- cut(x, breaks = seq(-90, 90, 90), include.lowest = TRUE)
    obl_bin_centers <- c(-45, 45)
    card_groups <- cut(x2, breaks = seq(-45, 180 - 45, 90), include.lowest = TRUE)
    card_bin_centers <- c(0, 90)
    angle_diff_fun <- angle_diff_180
    circ_sd_fun <- circ_sd_180
  } else if (space == "360") {
    x <- angle_diff_360(x, 0)
    x2 <- (x + 45) %% 360 - 45
    obl_groups <- cut(x, breaks = seq(-180, 180, 90), include.lowest = TRUE)
    obl_bin_centers <- seq(-135, 135, 90)
    card_groups <- cut(x2, breaks = seq(-45, 360 - 45, 90), include.lowest = TRUE)
    card_bin_centers <- seq(0, 270, 90)
    angle_diff_fun <- angle_diff_360
    circ_sd_fun <- circ_sd_360
  } else {
    stop("`space` argument should be 180 or 360.")
  }

  if (debug) {
    cat("N observation per group assuming cardinal bins: \n")
    cat(table(card_groups))
    cat("N observation per group assuming oblique bins: \n")
    cat(table(obl_groups))
  }
  for_fit <- data.table(
    x = x,
    x2 = x2,
    err, card_groups, obl_groups
  )
  if (missing(init_outliers)) {
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
      sigma_formula <- "~abs(dist_to_card)" # assumes that uncertainty changes linearly as a function of distance to cardinals regardless of the bias direction
      ll1 <- sum(for_fit[outlier == FALSE, logLik(gamlss::gamlss(err ~ poly(dist_to_card, var_sigma_poly_deg), sigma_formula, .SD, control = gam_ctrl)), by = .(card_groups)]$V1)
      ll2 <- sum(for_fit[outlier == FALSE, logLik(gamlss::gamlss(err ~ poly(dist_to_obl, var_sigma_poly_deg), sigma_formula, .SD, control = gam_ctrl)), by = .(obl_groups)]$V1)
    } else {
      ll1 <- sum(for_fit[outlier == FALSE, logLik(MASS::rlm(err ~ poly(x2, poly_deg))), by = .(card_groups)]$V1)
      ll2 <- sum(for_fit[outlier == FALSE, logLik(MASS::rlm(err ~ poly(x, poly_deg))), by = .(obl_groups)]$V1)
    }
    if (debug) {
      cat(sprintf("LL for bias type 1: %.2f, LL for bias type 2: %.2f", ll1, ll2))
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
  bin_centers <- break_points + (shift(break_points, 1, fill = break_points[length(break_points)]) - break_points) / 2
  bin_centers[1] <- angle_diff_fun(break_points[length(break_points)] + (break_points[1] + as.numeric(space) - break_points[length(break_points)]) / 2, 0)
  bin_width <- abs(angle_diff_fun(bin_boundaries, shift(bin_boundaries, 1))[2:length(bin_boundaries)])
  bin_width[1] <- break_points[1] + as.numeric(space) - break_points[length(break_points)]
  bin_labels <- sapply(2:(length(bin_boundaries)), \(i) {
    sprintf("[%.2f, %.2f]", bin_boundaries[i - 1], bin_boundaries[i])
  })
  bin_labels <- factor(bin_labels, levels = bin_labels)

  for_fit[, x_var := x]
  get_bin_i <- function(x, bin_centers, bin_width) {
    within_bin <- sapply(1:length(bin_centers), \(i) abs(angle_diff_fun(x, bin_centers[i])) <= (bin_width[i] / 2))
    max.col(within_bin, "first")
  }
  for_fit[, min_bp_i := get_bin_i(x, bin_centers, bin_width)]
  for_fit[, center_x := bin_centers[min_bp_i]]
  for_fit[, dc_var := angle_diff_fun(x, center_x)]
  for_fit[, gr_var := bin_labels[min_bp_i]]

  if (plots == "show" & debug == TRUE) {
    p_boundaries <- ggplot(for_fit, aes(x = x, y = err, color = gr_var)) +
      geom_point() +
      geom_vline(xintercept = angle_diff_fun(break_points, 0)) +
      geom_vline(color = "blue", xintercept = angle_diff_fun(bin_centers, 0))
    print(p_boundaries)
  }
  for_fit[, min_boundary_i := apply(sapply(break_points, \(bp) abs(angle_diff_fun(x, bp))), 1, which.min)]

  for_fit[, min_boundary_dist := angle_diff_fun(x, break_points[min_boundary_i])]
  for_fit[, bin_range := bin_width[min_bp_i]]
  for_fit[, bin_boundary_left := bin_boundaries[min_bp_i]]
  for_fit[, bin_boundary_right := bin_boundaries[min_bp_i + 1]]
  if (reassign_at_boundaries) {
    for_fit[, at_the_boundary := (abs(min_boundary_dist) - reassign_range) < (1e-12)]
  }

  if (var_sigma) {
    # get predictions
    if (reassign_at_boundaries) {
      if (debug) cat("Reassigning points at the boundaries...")
      if (any(for_fit[, unique(bin_range)] < (2 * reassign_range))) {
        stop("Reassignment range too large compared to bin sizes")
      }
      for_fit[, row_i := 1:.N]
      resid_at_boundaries <- for_fit[outlier == FALSE,
        get_boundary_preds(gr_var, copy(for_fit[outlier == FALSE]), space, reassign_range, gam_ctrl, ifelse(rep_n > 2, poly_deg, 1), angle_diff_fun),
        by = .(gr_var)
      ]
      resid_at_boundaries[, likelihood := dnorm(err, pred, pred_sigma, log = FALSE)]
      resid_at_boundaries[, new_weight := ifelse(at_the_boundary == FALSE, 1, likelihood / sum(likelihood)), by = .(err, x_var)]
      cur_weights <- resid_at_boundaries[at_the_boundary == TRUE, ]$new_weight
      stable_weights <- 0
      for (rep_n in 1:10) {
        weight_dt <- resid_at_boundaries[, .(row_i, gr_var, new_weight)]
        resid_at_boundaries <- resid_at_boundaries[,
          get_boundary_preds(gr_var, copy(for_fit[outlier == FALSE]), space, reassign_range, gam_ctrl, ifelse(rep_n > 2, poly_deg, 1), angle_diff_fun, weights = weight_dt),
          by = .(gr_var)
        ]
        resid_at_boundaries[, likelihood := dnorm(err, pred, pred_sigma, log = FALSE)]
        resid_at_boundaries[, new_weight := ifelse(at_the_boundary == FALSE, 1, likelihood / sum(likelihood)), by = .(err, x_var)]

        resid_at_boundaries_c <- dcast(resid_at_boundaries[at_the_boundary == TRUE], row_i ~ gr_var, value.var = "resid_at_boundaries")
        resid_at_boundaries_c <- resid_at_boundaries_c[, gr_var := names(.SD)[max.col(replace(-abs(.SD), is.na(.SD), -Inf))], .SDcols = 2:ncol(resid_at_boundaries_c)]
        prev_weights <- cur_weights

        cur_weights <- resid_at_boundaries[at_the_boundary == TRUE, ]$new_weight
        weights_change <- sum(abs(cur_weights - prev_weights))
        if (debug) {
          cat(sprintf("Reassignment step: %i; change in weights: %.5f", rep_n, weights_change))
        }
        for_fit[resid_at_boundaries_c, `:=`(gr_var = i.gr_var), on = .(row_i)]
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
              cat("Reassignment stopped at stable weights")
            }
            break
          }
        }
      }
    }
    for_fit[, dist_to_bin_centre := angle_diff_fun(x_var, center_x)]
    for_fit[, x_var := center_x + angle_diff_fun(x_var, center_x)]
    for_fit[, dc_var := angle_diff_fun(x, center_x)]
    if (debug) cat("Computing final fits...")

    likelihoods <- c()
    for (cg in unique(for_fit$gr_var)) {
      cur_df <- for_fit[gr_var == cg, .(err, x_var, dist_to_bin_centre, dc_var, outlier, dist_to_card)]
      fit <- gamlss::gamlss(err ~ pb(dist_to_bin_centre),
        ~ abs(dist_to_bin_centre),
        data = cur_df,
        weights = 1 - as.numeric(cur_df$outlier),
        control = gam_ctrl
      )

      if (debug) {
        cat("Fitted GAMLSS model coefficients\n")
        cat(coef(fit))
      }

      for_fit[gr_var == cg, pred := predict(fit, type = "response")]
      for_fit[gr_var == cg, pred_sigma := predict(fit, what = "sigma", type = "response")]

      for_fit[gr_var == cg, bias := err * sign(pred)]

      if (debug) {
        p_pred <- ggplot(for_fit[gr_var == cg], aes(x = dist_to_bin_centre, y = err)) +
          geom_point() +
          geom_line(aes(y = pred))
        print(p_pred)
      }
      for_fit[gr_var == cg, c("coef_sigma_int", "coef_sigma_slope") := data.frame(t(coef(fit, what = "sigma")))]
      likelihoods <- c(likelihoods, logLik(fit))
    }
  } else {
    for_fit[, pred := predict(MASS::rlm(err ~ poly(x_var, poly_deg), .SD[outlier == FALSE]), newdata = .SD[, .(x_var)]), by = .(card_groups)]
  }

  for_fit[, pred_lin := MASS::rlm(err ~ x_var)$fitted.values, by = .(gr_var)]

  for_fit[, be_c := err - pred]
  for_fit[, which_bin := as.numeric(gr_var)]
  for_fit[, center_y := predict(MASS::rlm(err ~ x_var),
    newdata = data.frame(x_var = center_x)
  ), by = .(gr_var)]

  if (var_sigma) {
    for_fit[, outlier := abs(be_c) > 3 * pred_sigma]
  } else {
    for_fit[, outlier := abs(be_c) > (3 * circ_sd_fun(be_c))]
  }

  if (plots %in% c("show", "return")) {
    for_fit[, outlier_f := factor(ifelse(outlier, "Outlier", "Non-outlier"))]
    sd_val <- for_fit[, circ_sd_fun(err)]
    plots_obj <- make_plots_of_biases(for_fit, poly_deg, sd_val)
    if (plots == "show") {
      print(plots_obj)
    } else {
      return(plots_obj)
    }
  }
  return(for_fit[, .(is_outlier = as.numeric(outlier), pred, be_c, which_bin, bias, bias_type, pred_lin, pred_sigma, coef_sigma_int, coef_sigma_slope, shifted_x = x_var, total_log_lik = sum(likelihoods))])
}

#' Remove cardinal biases for data with orientation (color, motion, ...) set in discrete steps
#'
#' @param err a vector of errors, deviations of response from the true stimuli
#' @param x a vector of true stimuli in degrees (see space)
#' @param space circular space to use (a string: `180` or `360`)
#' @param init_outliers a vector determining which errors are initially assumed to be outliers (default: NULL)
#'
#' @return returns a data.table with the following columns:
#' \itemize{
#' \item is_outlier - 0 for outliers (defined as ±3*predicted SD, where SD and mean are computed after excluding initially assumed outliers)
#' \item be_c error corrected for biases (`be_c = observed error - pred`)
#' }
#' @export
#'
remove_cardinal_biases_discrete <- function(err, x, space, init_outliers = NULL) {
  outlier <- be_c <- mean_err <- is_outlier <- . <- NULL # due to NSE notes in R CMD check

  if (space %in% c("180", "360")) {
    angle_diff_fun <- get(paste0("angle_diff_", space))
    circ_sd_fun <- get(paste0("circ_sd_", space))
    circ_mean_fun <- get(paste0("circ_mean_", space))
  } else {
    stop("`space` argument should be 180 or 360.")
  }

  data <- data.table(err, x)
  if (missing(init_outliers)) {
    data[, outlier := abs(err) > (3 * circ_sd_fun(err))]
  } else {
    data[, outlier := init_outliers]
  }
  data[, c("mean_err") := .(circ_mean_fun(err[outlier == FALSE])), by = x]
  data[, be_c := angle_diff_fun(err, mean_err), by = x]
  data[, is_outlier := abs(be_c) > (3 * circ_sd_fun(be_c[outlier == FALSE])), by = x]
  data[, .(be_c, is_outlier)]
}

#' Plots biases using the data from `remove_cardinal_biases`
#'
#' @param data data prepared by  `remove_cardinal_biases`
#' @param poly_deg the degree of polynomial
#' @param sd_val sd used to determine the outliers
#'
#' @return does not return anything
#'
#' @import ggplot2
#' @importFrom patchwork wrap_plots
#' @keywords internal
#'
make_plots_of_biases <- function(data, poly_deg, sd_val) {
  requireNamespace("patchwork")
  requireNamespace("ggplot2")

  outlier <- x_var <- gr_var <- err <- outlier_f <- pred <- pred_sigma <- be_c <- bias <- NULL # due to NSE notes in R CMD check

  common_plot_pars <- list(
    scale_x_continuous(breaks = seq(-180, 360, 90)),
    labs(x = "Orientation", shape = NULL, color = "Bin"),
    theme(legend.position = c(1, 1), legend.justification = c(1, 1)),
    guides(color = guide_legend(override.aes = list(size = 1)))
  )
  alpha <- 100 / data[, .N] * length(unique(data$gr_var))

  p1 <- ggplot(data = data[outlier == FALSE], aes(x = x_var, color = gr_var)) +
    geom_point(data = data, aes(y = err, shape = outlier_f), alpha = alpha) +
    geom_line(aes(y = pred), size = 1) +
    geom_line(aes(y = pred + 3 * pred_sigma)) +
    geom_line(aes(y = pred - 3 * pred_sigma)) +
    geom_hline(yintercept = c(-1, 1) * data[, 3 * sd_val], linetype = 2) +
    common_plot_pars +
    labs(y = "Error")

  p1a <- ggplot(data = data[outlier == FALSE], aes(x = x_var, color = gr_var)) +
    geom_point(data = data, aes(y = be_c, shape = outlier_f), alpha = alpha) +
    geom_hline(yintercept = c(-1, 1) * data[, 3 * circ_sd_180(be_c)], linetype = 2) +
    geom_line(data = data[outlier == FALSE], aes(y = 3 * pred_sigma)) +
    geom_line(data = data[outlier == FALSE], aes(y = -3 * pred_sigma)) +
    common_plot_pars +
    labs(y = "Bias-corrected error")

  p1b <- ggplot(data, aes(x = x_var, y = bias, color = gr_var, shape = outlier_f)) +
    geom_point(alpha = alpha) +
    common_plot_pars +
    geom_hline(yintercept = c(-1, 1) * data[, 3 * sd_val], linetype = 2) +
    theme(legend.position = "none") +
    labs(y = "Bias")

  p1 <- p1 + labs(shape = NULL) + guides(color = "none")
  patchwork::wrap_plots(p1, p1a, p1b, guides = "collect")
}

#' Pad circular data on both ends
#'
#' @param data data.table to pad
#' @param circ_var circular variable
#' @param circ_borders range of the circular variable
#' @param circ_part padding proportion
#' @param verbose print extra info
#'
#' @details Pads the data by adding a part of the data (default: 1/6th) from one end to another end. Useful to roughly account for circularity when using non-circular methods.
#' @return a padded data.table
#' @export
#'
#' @import data.table
#' @examples
#'
#' dt <- data.table::data.table(x = runif(1000, -90, 90), y = rnorm(1000))
#' pad_circ(dt, "x", verbose = TRUE)
#'
pad_circ <- function(data, circ_var, circ_borders = c(-90, 90), circ_part = 1 / 6, verbose = FALSE) {
  circ_range <- max(circ_borders) - min(circ_borders)

  data1 <- copy(data[get(circ_var) < (circ_borders[1] + circ_range * circ_part), ])
  data1[, (circ_var) := get(circ_var) + circ_range]

  data2 <- copy(data[get(circ_var) > (circ_borders[2] - circ_range * circ_part), ])
  data2[, (circ_var) := get(circ_var) - circ_range]
  if (verbose) {
    cat(sprintf("Rows in original DT: %i, padded on the left: %i, padded on the right: %i", data[, .N], data1[, .N], data2[, .N]))
  }

  rbind(data, data1, data2)
}

#' Get polynomial predictions for values at the boundaries
#'
#' A helper function for [remove_cardinal_biases()].
#'
#' @param group group (bin) id
#' @param data dataset
#' @param space see [remove_cardinal_biases()]
#' @param reassign_range see [remove_cardinal_biases()]
#' @param gam_ctrl control object for gam models
#' @param poly_deg see [remove_cardinal_biases()]
#' @param angle_diff_fun a function to compute difference between angles
#'
#' @return a data.table with predicted values
#'
#' @import data.table
#' @keywords internal
#'

get_boundary_preds <- function(group, data, space, reassign_range, gam_ctrl, poly_deg, angle_diff_fun, weights = NULL) {
  gr_var <- outlier <- err <- x_var <- dc_var <- center_x <- dist_to_card <- bin_boundary_left <- bin_boundary_right <- bin_range <- dist_to_bin_centre <- row_i <- at_the_boundary <- dist_to_boundary <- dist_to_boundary_norm <- new_weight <- weight <- pred <- . <- predict <- pred_sigma <- resid_at_boundaries <- NULL # due to NSE notes in R CMD check
  cur_df <- data[gr_var == group & outlier == FALSE, .(err, x_var, dc_var,
    dist_to_bin_centre = angle_diff_fun(x_var, center_x), weight = NULL, adc = abs(dist_to_card), center_x, bin_boundary_left, bin_boundary_right, bin_range
  )]

  curr_bin_range <- cur_df$bin_range[[1]]
  curr_bin_center <- cur_df$center_x[[1]]
  boundary1 <- cur_df$bin_boundary_left[[1]]
  boundary2 <- cur_df$bin_boundary_right[[1]]

  data[, dist_to_bin_centre := angle_diff_fun(x_var, curr_bin_center)]
  data_incl_boundaries <- data[
    outlier == FALSE &
      abs(dist_to_bin_centre) < (curr_bin_range / 2 + reassign_range + 1e-12),
    .(
      row_i, err, x_var, dc_var, gr_var,
      dist_to_bin_centre,
      at_the_boundary,
      center_x
    )
  ]
  data_incl_boundaries[, dist_to_boundary := (abs(dist_to_bin_centre) - curr_bin_range / 2)]
  data_incl_boundaries[, dist_to_boundary_norm := (dist_to_boundary + reassign_range) / (2 * reassign_range)]

  if (!missing(weights)) {
    data_incl_boundaries[, gr_var := group]
    data_incl_boundaries[weights, `:=`(weight = new_weight), on = .(row_i, gr_var)]
  } else {
    data_incl_boundaries[, weight := ifelse(at_the_boundary == FALSE, 1, ifelse(gr_var == group, 0.75, 0.25))]
  }


  fit <- gamlss::gamlss(err ~ dist_to_bin_centre,
    ~ abs(dist_to_bin_centre),
    data = data_incl_boundaries,
    weights = weight,
    control = gam_ctrl
  )

  data_incl_boundaries[, pred := predict(fit, type = "response")]
  data_incl_boundaries[, pred_sigma := predict(fit, what = "sigma", type = "response")]

  data_incl_boundaries[, resid_at_boundaries := err - pred]
  data_incl_boundaries[, .(row_i, at_the_boundary, x_var, dc_var, dist_to_bin_centre, err, pred, resid_at_boundaries, dist_to_boundary, dist_to_boundary_norm, weight, pred_sigma)]
}


a_fun <- function(x) {
  besselI(x, 1, expon.scaled = TRUE) / besselI(x, 0, expon.scaled = TRUE)
}

inverse <- function(f, lower = 1e-16, upper = 1000) {
  function(y) stats::uniroot((function(x) f(x) - y), lower = lower, upper = upper, extendInt = "yes")[[1]]
}

#' Conversion between the circular SD and kappa of von Mises
#'
#' @param kappa von Mises kappa parameter
#' @param sd_deg circular SD of von Mises (degrees)
#' @param sd circular SD of von Mises (radians)
#'
#' @return `vm_kappa_to_circ_sd` and `vm_kappa_to_circ_sd_deg` return circular SD (in radians or degrees, respectively) corresponding to a given kappa. `vm_circ_sd_to_kappa` and `vm_circ_sd_deg_to_kappa` return kappa corresponding to a given circular SD (in radians or degrees, respectively).
#'
#' @export
#'
#' @examples
#'
#' vm_kappa <- 5
#' vm_sd <- vm_kappa_to_circ_sd(vm_kappa)
#'
#' vm_circ_sd_to_kappa(vm_sd)
#'
#' x <- circular::rvonmises(10000, mu = circular::circular(0), kappa = vm_kappa)
#'
#' sprintf("Expected SD: %.2f, actual SD: %.2f", vm_sd, circ_sd_rad(x))
#'
vm_kappa_to_circ_sd <- function(kappa) {
  sqrt(-2 * log(a_fun(kappa)))
}


#' @export
#' @describeIn vm_kappa_to_circ_sd get circular SD (in degrees) from kappa
vm_kappa_to_circ_sd_deg <- function(kappa) {
  vm_kappa_to_circ_sd(kappa) / pi * 180
}

#' @export
#' @describeIn vm_kappa_to_circ_sd get kappa from circular SD (in radians)

vm_circ_sd_to_kappa <- function(sd) {
  vm_circ_sd_inverse <- inverse(vm_kappa_to_circ_sd)
  sapply(sd, function(x) tryCatch(vm_circ_sd_inverse(x), error = function(e) paste("Can't convert sigma = ", x, " to kappa, error ", e)))
}


#' @export
#' @describeIn vm_kappa_to_circ_sd get kappa from circular SD (in degrees)

vm_circ_sd_deg_to_kappa <- function(sd_deg) {
  vm_circ_sd_to_kappa(sd_deg / 180 * pi)
}



#' Weighted standard error of the mean (SEM_w)
#'
#' Computes the variance of a weighted mean following the definitions given by Kirchner (2006).
#' @param x variable to compute the SEM for
#' @param w weights
#' @param na.rm should NAs be removed
#'
#' @details
#' James Kirchner describes two different cases when the weighted variance is computed. The code here implements Case I where "one wants to give more weight to some points than to others, because they are considered to be more important" and "the weights differ but the uncertainties associated with the individual xi are assumed to be the same" (Kirchner, 2006, p. 1). The formula used is:
#' \mjsdeqn{SEM_w = \sqrt{\left(\sum_{i = 1}^{N} (w_{i} x_i^2)-\bar{x}^2\right)\frac{\sum_{i = 1}^{N} w_i^2}{1-\sum_{i = 1}^{N} w_i^2}} }
#' The expected error is within 5% of the bootstrapped SEM (at larger sample sizes).
#'
#' @return weighted standard error of the mean
#' @references {
#' \itemize{
#' \item Kirchner, J. 2006. Data Analysis Toolkit #12: Weighted averages and their uncertainties. \url{https://seismo.berkeley.edu/~kirchner/Toolkits/Toolkit_12.pdf}.  Retrieved on 04.07.2024.
#' \item Bevington, P. R. 1969. Data Reduction and Error Analysis for the Physical Sciences. McGraw-Hill, 336 pp.
#'
#' }
#' }
#' @export
#' @examples
#' set.seed(1)
#' n_obs <- 200
#' w <- runif(n_obs)
#' w <- w / sum(w)
#' x <- rnorm(n_obs, sd = 5)
#' weighted_sem(x, w)
weighted_sem <- function(x, w, na.rm = FALSE) {
  if (na.rm) {
    w <- w[i <- !is.na(x)]
    x <- x[i]
  }
  n <- length(w)
  w <- w / sum(w)
  x_w_bar <- stats::weighted.mean(x, w, na.rm = na.rm)
  # case I
  var_w_1 <- (sum(w * x^2) - x_w_bar^2) / (1 - sum(w^2))
  sem_w_1 <- sqrt(var_w_1 * sum(w^2))
  # Kirchner also provides case II that assumes equal importance but different variances
  # var_w_2 <- (sum(w*x^2)-x_w_bar^2)*n/(n-1)
  # sem_w_2 <- sqrt(var_w_2/n)

  return(sem_w_1)
}

#' Compute predictions for circular LOESS
#' @param object a circular LOESS object
#' @param newdata a data.frame with a variable x on which the predictions are computed
#' @param ... other arguments (passed to circ_loess)
#'
#'
#' @return a data.frame with predictions
#' @method predict circ_loess
#' @export
#' @keywords internal
#'
predict.circ_loess <- function(object, newdata, ...) {
  res <- circ_loess(angle = object$angle, y = object$y, xseq = newdata$x, circ_space = object$circ_space, span = object$span, ...)
  list(fit = data.frame(y = res$y_est, ymin = res$y_est - 1.96 * res$y_se, ymax = res$y_est + 1.96 * res$y_se), se.fit = res$y_se)
}




#' An implementation of circular-linear locally-weighted regression (LOESS)
#' \loadmathjax
#'
#' Provides an locally-weighted average when the independent variable is circular and depended variable is linear. Mainly to use with ggplot2.
#'
#' @param formula the formula, e.g., y ~ x
#' @param data data to use
#' @param angle a vector of angles (not used if a formula is provided)
#' @param y dependent variable vector (not used if a formula is provided)
#' @param xseq a grid to compute predictions on (optional, the default is to use 500 points spanning the circle)
#' @param circ_space circular space to use (90, 180, 360, or 2*pi)
#' @param span a span to adjust the degree of smoothing
#' @param ... other arguments (ignored)
#'
#' @details Weights for the regression are computed as
#' \mjsdeqn{w = (1-(d/d_{max})^3)^3}
#' where _d_ is the angular difference between the point at which the estimate is computed and the angles in the data, and \mjseqn{d_{max}} is the maximum possible distance. If `span` is above 1, all points are included and \mjseqn{d_{max} = {circ\_space}/(4*span)}. Otherwise, a proportion \mjseqn{\alpha} of the points included based on their distance to the point at which the estimate is computed and \mjseqn{d_{max}} is the corresponding maximal distance.

#' @return an object of `circ_loess` class with the following parameters:
#' * `angle` the angles in the data
#' * `y` the dependent variable vales in the data
#' * `xseq` the grid on which the loess values are estimated
#' * `y_est` the estimated loess values
#' * `y_se` standard errors
#' * `w` weights
#' * `circ_space` circular space used
#' * `span` span used
#'
#' @seealso [stats::loess()]
#'
#' @export
#'
#' @examples
#' p <- ggplot(Pascucci_et_al_2019_data, aes(x = orientation, y = err)) +
#'   geom_point(alpha = 0.05) +
#'   labs(x = "Orientation, deg.", y = "Error, deg.")
#' p1 <- p + geom_smooth(method = "loess") + ggtitle("Standard LOESS")
#' p2 <- p + geom_smooth(method = "circ_loess", method.args = list(circ_space = 180, span = 0.5)) +
#'   ggtitle("Circular LOESS, span = 0.5")
#' p3 <- p + geom_smooth(method = "circ_loess", method.args = list(circ_space = 180, span = 0.2)) +
#'   ggtitle("Circular LOESS, span = 0.2")
#' (p1 + p2 + p3)
#'
circ_loess <- function(formula = NULL, data = NULL, angle = NULL, y = NULL, xseq = NULL, circ_space = NULL, span = 0.75, ...) {
  if (!is.null(formula)) {
    M <- stats::model.frame(formula, data)
    angle <- M[, 2]
    y <- M[, 1]
  }

  if (is.null(circ_space)) {
    message("circular space is not set, tryin to guess based on the data (prone to errors)...")
    if (max(abs(angle)) > 90) {
      circ_space <- 360
      message("circ_loess assuming 360 deg. space")
    } else if (max(abs(angle)) > 45) {
      circ_space <- 180
      message("circ_loess assuming 180 deg. space")
    } else if (max(abs(angle)) > 15) {
      circ_space <- 90
      message("circ_loess assuming 90 deg. space")
    } else {
      circ_space <- 2 * pi
      message("circ_loess assuming 2pi space")
    }
  }
  if (circ_space %in% c(90, 180, 360)) {
    diff_fun <- get(paste0("angle_diff_", circ_space))
  } else if (circ_space == 2 * pi) {
    diff_fun <- angle_diff_rad
  } else {
    stop("Unknown circ_space value. Should be 90, 180, 360, or 2*pi.")
  }

  range <- c(-0.5, 0.5) * circ_space

  angle <- diff_fun(angle, 0)

  if (is.null(xseq)) {
    xseq <- seq(min(angle), max(angle), length.out = 500)
  }

  y_est <- sapply(xseq, function(x) {
    dist <- abs(diff_fun(x, angle))
    if (span < 1) {
      included_obs <- dist <= stats::quantile(dist, span)
      angle <- angle[included_obs]
      dist <- dist[included_obs]
      y <- y[included_obs]
      max_dist <- max(dist)
    } else {
      max_dist <- diff(range) / 2 * span
    }


    w <- (1 - (dist / max_dist)^3)^3

    list(stats::weighted.mean(y, w), weighted_sem(y, w)^0.5, w)
  })
  structure(list(
    angle = angle, y = y, xseq = xseq, y_est = unlist(y_est[1, ]), circ_space = circ_space, span = span,
    y_se = unlist(y_est[2, ]), w = unlist(y_est[3, ])
  ), class = "circ_loess")
}


#' Compute asymmetry in weighted probability density
#'
#' This function calculates the asymmetry in the probability density of a given variable (usually errors) relative to another variable (usually dissimilarity) using kernel density estimation. The asymmetry is computed for each x-axis value, and the result can be averaged or returned for each value individually.
#'
#' @param dt data.table with the data.
#' @param circ_space Circular space, which can be 180 or 360 (default: 180).
#' @param weights_sd Standard deviation of the Gaussian window to use across `xvar` (default: 10).
#' @param kernel_bw Bandwidth for the kernel density estimator across `yvar`. If NULL, it is computed using [stats::bw.SJ()] (default: NULL).
#' @param xvar X-axis variable, such as dissimilarity between items (default: "abs_td_dist").
#' @param yvar Y-axis variable, normally errors (default: "bias_to_distr_corr").
#' @param by A vector of grouping variable names (default: an empty vector).
#' @param n The number of steps for the x-axis variable at which the density is computed (default: 181).
#' @param average If TRUE, the asymmetry is averaged for each x-value (default: TRUE).
#' @param return_full_density If TRUE, returns the full data.table with density computed at each point (default: FALSE).
#' @return A data.table with the grouping variables, `dist` - the values of X-axis variable at which the density is computed, and `delta` - the difference (asymmetry) in probability density for positive and negative values of `yvar`; or the full density data if `return_full_density` is TRUE.
#' @export
#' @importFrom stats as.formula bw.SJ density weights
#' @import data.table
#'
#' @examples
#'
#' data(Pascucci_et_al_2019_data)
#' ex_data <- Pascucci_et_al_2019_data
#' ex_data[, err := angle_diff_180(reported, orientation)] # response errors
#' ex_data[, prev_ori := shift(orientation), by = observer] # orientation on previous trial
#'
#' # determine the shift in orientations between trials
#' ex_data[, diff_in_ori := angle_diff_180(prev_ori, orientation)]
#' ex_data[, abs_diff_in_ori := abs(diff_in_ori)]
#' ex_data[, err_rel_to_prev_targ := ifelse(diff_in_ori < 0, -err, err)]
#'
#' err_dens <- density_asymmetry(ex_data[!is.na(err_rel_to_prev_targ)],
#'   circ_space = 180, weights_sd = 10, xvar = "abs_diff_in_ori",
#'   yvar = "err_rel_to_prev_targ", by = c("observer")
#' )
#'
#' ggplot(err_dens, aes(x = dist, y = delta)) +
#'   geom_line(stat = "summary", fun = mean) +
#'   labs(y = "Asymmetry in error probability density, %", x = "Absolute orientation difference, °")
#'
density_asymmetry <- function(dt, circ_space = 180, weights_sd = 10, kernel_bw = NULL, xvar = "abs_td_dist", yvar = "bias_to_distr_corr", by = c(), n = 181, average = T, return_full_density = F) {
  x_val <- x <- x_sign <- delta <- `1` <- `-1` <- total <- ratio <- . <- NULL # due to NSE notes in R CMD check

  if (!(circ_space %in% c(180, 360))) {
    stop("`circ_space` should be 180 or 360")
  }

  max_diss <- circ_space / 2
  if (is.null(kernel_bw)) {
    kernel_bw <- bw.SJ(dt[, get(yvar)])
  }

  res <- rbindlist(sapply(1:max_diss, \(i) {
    dt[, weights := dnorm(get(xvar), mean = i, sd = weights_sd)]

    res <- dt[, density(get(yvar),
      from = -max_diss, to = max_diss, n = n, bw = kernel_bw,
      weights = weights / sum(weights)
    )[c("x", "y")], by = by]
    res[, x_val := abs(x)]
    res[, x_sign := sign(x)]
    res$dist <- i

    if (return_full_density) {
      return(res)
    }

    res <- dcast(res,
      as.formula(paste(paste(by, collapse = "+"), "+dist+x_val~x_sign")),
      value.var = "y"
    )
    res[, delta := `1` - `-1`]
    res[, total := `1` + `-1`]
    res[, ratio := `1` / `-1`]

    res
  }, simplify = F))

  if (return_full_density) {
    return(res)
  }

  if (average) {
    res <- res[!is.na(delta), .(delta = sum(delta)), by = c("dist", by)]
  }

  attr(res, "kernel_bw") <- kernel_bw

  res
}

#' Compute asymmetry in discrete weighted probability density
#'
#' This function calculates the asymmetry in the probability density of a given variable using kernel density estimation. Unlike `density_asymmetry`, it does not take an `xvar` and assumes that the asymmetry is calculated for the whole dataset or subsets defined by `by` (which could also include a discrete x variable if needed) .
#'
#' @param dt data.table with the data.
#' @param yvar Y-axis variable, normally errors (default: "bias_to_distr_corr").
#' @param circ_space Circular space, which can be 180 or 360 (default: 180).
#' @param kernel_bw Bandwidth for the kernel density estimator across `yvar`. If NULL, it is computed using [stats::bw.SJ()] (default: NULL).
#' @param by A vector of grouping variable names (default: an empty vector).
#' @param n The number of steps for the density computation (default: 181).
#' @param average If TRUE, the asymmetry is averaged (default: TRUE).
#' @param return_full_density If TRUE, returns the full data.table with density computed at each point (default: FALSE).
#' @return A data.table with the grouping variables and `delta` - the difference (asymmetry) in probability density for positive and negative values of `yvar`; or the full density data if `return_full_density` is TRUE.
#' @export
#' @importFrom stats as.formula bw.SJ density weights
#' @import data.table
#'
#' @examples
#'
#' data(Pascucci_et_al_2019_data)
#' ex_data <- Pascucci_et_al_2019_data
#' ex_data[, err := angle_diff_180(reported, orientation)] # response errors
#' ex_data[, prev_ori := shift(orientation), by = observer] # orientation on previous trial
#'
#' # determine the shift in orientations between trials
#' ex_data[, diff_in_ori := angle_diff_180(prev_ori, orientation)]
#' ex_data[, abs_diff_in_ori := abs(diff_in_ori)]
#' ex_data[, err_rel_to_prev_targ := ifelse(diff_in_ori < 0, -err, err)]
#' ex_data[, similarity_discrete := ifelse(abs_diff_in_ori>45, 'Dissimilar', 'Similar')]
#' err_dens_discrete <- density_asymmetry_discrete(ex_data[!is.na(err_rel_to_prev_targ)],
#'   circ_space = 180, yvar = "err_rel_to_prev_targ", by = c("observer","similarity_discrete")
#' )
#'
#' ggplot(err_dens_discrete, aes(x = similarity_discrete, y = delta))+
#' geom_violin() + geom_point() +
#'   labs(y = "Asymmetry in error probability density, %", x = "Previous target")
#'


density_asymmetry_discrete <- function(dt, yvar = "bias_to_distr_corr", circ_space = 180, kernel_bw = NULL, by = c(), n = 181, average = T, return_full_density = F) {
  x_val <- x <- x_sign <- delta <- `1` <- `-1` <- total <- ratio <- . <- NULL # due to NSE notes in R CMD check

  if (!(circ_space %in% c(180, 360))) {
    stop("`circ_space` should be 180 or 360")
  }

  max_diss <- circ_space / 2
  if (is.null(kernel_bw)) {
    kernel_bw <- bw.SJ(dt[, get(yvar)])
  }

  res <- dt[, density(get(yvar),
                      from = -max_diss, to = max_diss, n = n, bw = kernel_bw,
  )[c("x", "y")], by = by]
  res[, x_val := abs(x)]
  res[, x_sign := sign(x)]

  if (return_full_density) {
    return(res)
  }

  res <- dcast(res,
               as.formula(paste(paste(by, collapse = "+"), "+x_val~x_sign")),
               value.var = "y"
  )
  res[, delta := `1` - `-1`]
  res[, total := `1` + `-1`]
  res[, ratio := `1` / `-1`]


  if (average) {
    res <- res[!is.na(delta), .(delta = sum(delta)), by = by]
  }

  attr(res, "kernel_bw") <- kernel_bw

  res
}
