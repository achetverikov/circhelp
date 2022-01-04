
#' Circular mean
#'
#' @param x vector of values
#'
#' @return mean of values in the vector
#' @export
#'
#' @examples
#' x <- runif(1000,-pi,pi)
#' mean(x)
#' circ_mean_rad(x)
#'
#' @describeIn circ_mean_rad circular mean in 2pi space
circ_mean_rad<-function (x){
  # from CircStats package
  sinr <- sum(sin(x))
  cosr <- sum(cos(x))
  circmean <- atan2(sinr, cosr)
  circmean
}

#' @describeIn circ_mean_rad circular mean in 180° space (e.g., line orientation)
#' @export
circ_mean_180<-function(x){
  circ_mean_rad(x/90*pi)/pi*90
}

#' @describeIn circ_mean_rad circular mean in 360° space
#' @export
circ_mean_360<-function(x){
  circ_mean_rad(x/180*pi)/pi*180
}

#' Differences between angles in different circular spaces
#'
#' @param a first angle
#' @param b second angle
#' @details By default, all functions return values in ± half-range space (e.g., -pi to pi for 2pi radian space used by angle_diff_rad) but angle_diff_180_45 and angle_diff_360_90 return values in \[-1/4 range, 3/4 range\] space
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
angle_diff_rad<-function(a,b){
  c = a - b
  (c+pi)%%(2*pi) - pi
}

#' @describeIn angle_diff_rad angle difference in 360 degree space
#' @export

angle_diff_360<-function(a,b){
  angle_diff_rad(a/180*pi, b/180*pi)/pi*180
}

#' @describeIn angle_diff_rad angle difference in 180 degree space (e.g., line orientation)
#' @export
#'
angle_diff_180<-function(a,b){
  angle_diff_rad(a/90*pi, b/90*pi)/pi*90
}

#' @describeIn angle_diff_rad angle difference in 90 degree space
#' @export

angle_diff_90<-function(a,b){
  angle_diff_rad(a/45*pi, b/45*pi)/pi*45
}

#' @describeIn angle_diff_rad angle difference in 180 degree space from -45 to 135
#' @export

angle_diff_180_45<-function(a,b){
  c = a - b
  (c+45)%%180 - 45
}

#' @describeIn angle_diff_rad angle difference in 360 degree space from -90 to 270
#' @export

angle_diff_360_90<-function(a,b){
  c = a - b
  (c+90)%%360 - 90
}

#' Circular correlation coefficient
#'
#' Computes a circular correlation coefficient as defined in Jammalamadaka & SenGupta (2001).
#' @param a first variable
#' @param b second variable
#' @param ill_defined is one of the variables mean is not well-defined (e.g., it is uniformly distributed)?
#' @param mu fix the mean parameter of both vectors to a certain value
#'
#' @return correlation value
#' @references {
#' Jammalamadaka, S. R., & SenGupta, A. (2001). Topics in Circular Statistics. WORLD SCIENTIFIC. https://doi.org/10.1142/4031
#' }
#' @export
#'
#' @examples
#' data <- rmvn(10000, c(0,0), V = matrix(c(1,0.5,0.5,1), ncol = 2))
#' circ_corr(data[,1], data[,2])


circ_corr <- function(a, b, ill_defined = FALSE, mu = NULL){
  mu_a <- circ_mean_rad(a)
  mu_b <- circ_mean_rad(b)
  if (!is.null(mu)){
    mu_a = mu_b = mu
  } else if (ill_defined){
    mean_diff <- circ_mean_rad(a-b)
    mean_sum <- circ_mean_rad(a+b)
    mu_a <- (mean_diff + mean_sum)/2
    mu_b <- (mean_sum - mean_diff)/2
  }
  sin_a <- sin(a - mu_a)
  sin_b <- sin(b - mu_b)
  rho <- sum(sin_a*sin_b)/sqrt(sum(sin_a*sin_a)*sum(sin_b*sin_b))
  rho
}

#' Weighted circular parameters
#'
#' @param x vector of values (in radians)
#' @param w vector of weights
#'
#' @return weighted mean of values in the vector
#' @export
#'
#' @examples
#' x <- rnorm(1000,0, 0.5)
#' w <- runif(1000, 0, 1)
#' weighted.mean(x, w)
#' weighted_circ_mean(x, w)
#'
#' @describeIn weighted_circ_mean weighted circular mean

weighted_circ_mean<-function(x, w){
  if (length(w)!=length(x))
    stop('Weights (w) should have the same length as values (x)')

  sum_w<-sum(w)
  atan2(sum(w*sin(x))/sum_w, sum(w*cos(x))/sum_w)

}

#' @describeIn weighted_circ_mean an alternative way to compute weighted circular mean (the results are the same)
#' @export
weighted_circ_mean2 <- function(x, w){
  if (length(w)!=length(x))
    stop('Weights (w) should have the same length as values (x)')
  z = exp(1i*x)
  Arg(sum(w*z)/sum(w))
}

#' @describeIn weighted_circ_mean weighted circular SD
#' @export

weighted_circ_sd<-function(x, w){
  sum_w<-sum(w)

  r <- sqrt((sum(w*sin(x))/sum_w)^2+(sum(w*cos(x))/sum_w)^2)
  sqrt(-2*log(r))
}

#' @describeIn weighted_circ_mean weighted mean resultant length
#' @export

weighted_circ_rho<-function(x, w){
  sum_w<-sum(w)

  r <- sqrt((sum(w*sin(x))/sum_w)^2+(sum(w*cos(x))/sum_w)^2)
  r
}


#' Correct angle value
#'
#' @param x angle
#'
#' @return angle in \[-pi, pi\] space
#' @export
#'
#' @examples
#' correct_angle_rad(4*pi)
#'
correct_angle_rad <- function(x) {
  angle_diff_rad(x, 0)
}

#' Circular standard deviation
#'
#' @param x vector of angles
#'
#' @return standard deviation of values in the vector
#' @export
#'
#' @examples
#' circ_sd_rad(rnorm(50))
#' circ_sd_180(rnorm(50))
#'
#' @describeIn circ_sd_rad SD of angles in radians
circ_sd_rad <- function(x){
  r = sum(exp(1i*x));

  # mean resultant length
  rho = abs(r)/length(x);

  sqrt(-2*log(rho))
}

#' @describeIn circ_sd_rad SD of angles in 360 degree space
#' @export

circ_sd_360 <- function(x){
  circ_sd_rad(x/180*pi)/pi*180
}
#' @describeIn circ_sd_rad SD of angles in 180 degree space
#' @export
circ_sd_180 <- function(x){
  circ_sd_rad(x/90*pi)/pi*90
}

circ_dist <- function(x, y){
  angle(exp(1i*x)/exp(1i*y))
}

angle <- function(x){
  atan2(Im(x), Re(x))
}

#' A set of descriptive statistics for circular data
#'
#' @param x vector of angles
#' @param w weights for the values in the vector
#' @param d correction for the bias for data with known spacing
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
#' x <- c(rnorm(50,0,0.5), rnorm(20,1,0.5))
#' circ_descr(x)
#'
circ_descr <- function(x, w = NULL, d = NULL){
  if (is.null(w))
    w = rep(1, length(x))
  # compute weighted sum of cos and sin of angles
  r = sum(w*exp(1i*x));
  mu = angle(r);

  # mean resultant length
  rho = abs(r)/sum(w);
  # for data with known spacing, apply correction factor to correct for bias
  # in the estimation of r (see Zar, p. 601, equ. 26.16)
  if (!is.null(d)){
    c = d/2/sin(d/2);
    rho = c*rho;
  }

  # 2nd moment
  alpha_rel = angle(exp(1i*x)/exp(1i*mu))

  p = 2
  cbar = mean(cos(p*alpha_rel)*w)
  sbar = mean(sin(p*alpha_rel)*w)
  mp = cbar + 1i*sbar

  rho2 = abs(mp) # shouldn't rho2 be adjusted for d as well?
  mu2 = angle(mp)

  b = sum(w*sin(2*circ_dist(x,mu)))/sum(w)
  b0 = rho2*sin(circ_dist(mu2,2*mu))/(1-rho)^(3/2)

  b_rel_to_zero = sum(w*sin(2*circ_dist(x,0)))/sum(w)
  s0 = sqrt(-2*log(rho))
  list(mu = mu, sigma = s0, skew_pewsey = b, skew_fischer = b0, rho = rho, skew_rel_to_zero = b_rel_to_zero)
}

#' Remove cardinal biases
#'
#' @param err a vector of errors, deviations of response from the true stimuli
#' @param x a vector of true stimuli in degrees (see space)
#' @param space circular space to use (a string: '180' or '360')
#' @param bias_type bias type to use ('fit', 'card', or 'obl', see details)
#' @param do_plot show plots (default: False)
#' @param poly_deg degree of the fitted polynomials for each bin (default: 4)
#' @param var_sigma allow standard deviation (width) of the fitted response distribution to vary as a function of distance to the nearest cardinal (default: True)
#' @param var_sigma_poly_deg degree of the fitted polynomials for each bin for the first approximation for the response distribution to select the best fitting model (default: 4)
#' @param reassign_at_boundaries select the bin for the observations at the boundaries between bins based on the best-fitting polynomial (default: True)
#' @param reassign_range maximum distance to the boundary at which reassignment can occur (default: 2 degrees)
#' @param debug print some extra info (default: False)
#'
#' @details
#' If the bias_type is set to 'fit', the function computes the cardinal biases in the following way:
#' \enumerate{
#' \item Create two sets of bins, splitting the stimuli vector into bins centered at cardinal and at oblique directions.
#' \item For each set of bins, fit a nth-degree polynomial for the responses in each bin, optionally allowing the distribution of responses to vary in width as a function of distance to the nearest cardinal (regardless of whether the bins are centered at the cardinal or at the oblique, the width of the response distribution usually increases as the distance to cardinals increase).
#' \item Choose the best-fitting model between the one using cardinal and the one using oblique bins.
#' \item Compute the residuals of the best-fitting model - that's your bias-corrected error - and the biases (see below).
#' }
#' The bias is computed by flipping the sign of errors when the average predicted error is negative, so, that, for example, if on average the responses are shifted clockwise relative to the true values, the trial-by-trial error would count as bias when it is also shifted clockwise.
#'
#' If bias_type is set to 'obl' or 'card', only one set of bins is used, centred at cardinal or oblique angles, respectively.
#'
#'
#' @return a list with the follwing elements:
#' \itemize{
#' \item is_outlier - 0 for outliers (defined as ±3*pred_sigma for the model with varying sigma or as ±3\*SD for the simple model)
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
#' @export
#'
#' @examples
#'
#' # Data in orientation domain from Pascucci et al. (2019, PLOS Bio),
#' # https://doi.org/10.5281/zenodo.2544946
#'
#' data <- fread('https://zenodo.org/record/2544946/files/Experiment2_rawdata.csv?download=1')
#' data[, err:=angle_diff_180(reported, orientation)]
#' data[observer==4, remove_cardinal_biases(err, orientation, do_plot = TRUE)]
#'
#' # Data in motion domain from Bae & Luck (2018, Neuroimage), https://osf.io/2h6w9/
#' data_motion <- fread('https://osf.io/4m2kb/download')
#' data_motion[, err:=angle_diff_360(RespAngle, TargetDirection)]
#' data_motion[subject_Num == unique(subject_Num)[5],
#' remove_cardinal_biases(err, TargetDirection, space = '360', do_plot = TRUE)]
#'
remove_cardinal_biases <- function(err, x, space = '180', bias_type = 'fit', do_plot = FALSE, poly_deg = 4,  var_sigma = TRUE, var_sigma_poly_deg = 4, reassign_at_boundaries = TRUE, reassign_range = 2, debug = FALSE){
  require(MASS)
  require(mgcv)
  require(ggplot2)
  require(gamlss)
  require(data.table)

  if (!(bias_type %in% c('fit','card','obl')) ){
    stop("`bias_type` should be 'fit','card', or 'obl'" )
  }

  if (space == '180'){
    x <- angle_diff_180(x,0)
    x2 <- angle_diff_180_45(x,0)
    obl_groups <- cut(x, breaks = seq(-90, 90, 90), include.lowest = T)
    obl_bin_centers <- c(-45, 45)
    card_groups <- cut(x2, breaks = seq(-45, 180-45, 90), include.lowest = T)
    card_bin_centers <- c(0, 90)
    angle_diff_fun <- angle_diff_180
    circ_sd_fun <- circ_sd_180

  } else if (space == '360'){
    x <- angle_diff_360(x,0)
    x2 <- (x+45)%%360 - 45
    obl_groups <- cut(x, breaks = seq(-180, 180, 90), include.lowest = T)
    obl_bin_centers <- seq(-135, 135, 90)
    card_groups <- cut(x2, breaks = seq(-45, 360-45, 90), include.lowest = T)
    card_bin_centers <- seq(0, 270, 90)
    angle_diff_fun <- angle_diff_360
    circ_sd_fun <- circ_sd_360
  } else stop('`space` argument should be 180 or 360.')

  if (debug){
    print(table(card_groups))
    print(table(obl_groups))
  }
  for_fit<-data.table(x = x,
                      x2 = x2,
                      err, card_groups, obl_groups)
  for_fit[,outlier:=abs(err)>3*circ_sd_fun(err)]
  for_fit[,dist_to_card:=angle_diff_90(x2, 0)]
  for_fit[,dist_to_obl:=angle_diff_90(x, 45)]
  gam_ctrl <- gamlss.control(trace = F)

  if (bias_type == 'fit'){
    if (var_sigma){
      sigma_formula <- '~abs(dist_to_card)' # assumes that uncertainty changes linearly as a function of distance to cardinals regardless of the bias direction
      ll1<-sum(for_fit[outlier==F,logLik(gamlss(err~poly(dist_to_card, var_sigma_poly_deg), sigma_formula,.SD, control = gam_ctrl)), by=.(card_groups)]$V1)
      ll2<-sum(for_fit[outlier==F,logLik(gamlss(err~poly(dist_to_obl, var_sigma_poly_deg), sigma_formula, .SD, control = gam_ctrl)), by=.(obl_groups)]$V1)
    }
    else {
      ll1<-sum(for_fit[outlier==F,logLik(rlm(err~poly(x2,poly_deg))), by=.(card_groups)]$V1)
      ll2<-sum(for_fit[outlier==F,logLik(rlm(err~poly(x,poly_deg))), by=.(obl_groups)]$V1)
    }
    if (debug)
      print(sprintf('LL for bias type 1: %.2f, LL for bias type 2: %.2f', ll1, ll2))

    if (ll1>=ll2){
      bias_type = 'card'
    } else bias_type = 'obl'

  }

  for_fit$pred_sigma<-NA_real_
  for_fit$coef_sigma_int<-NA_real_
  for_fit$coef_sigma_slope<-NA_real_
  if (bias_type == 'card'){
    for_fit[,x_var:=x2]
    for_fit[,dc_var:=dist_to_card]
    for_fit[,gr_var:=card_groups]
    for_fit[,at_the_boundary:=(abs(dist_to_obl)-reassign_range)<(1e-12)]

    bin_centers <- card_bin_centers
  } else {
    for_fit[,x_var := x]
    for_fit[,dc_var := dist_to_obl]
    for_fit[,gr_var := obl_groups]
    for_fit[,at_the_boundary:=(abs(dist_to_card)-reassign_range)<(1e-12)]

    bin_centers <- obl_bin_centers

  }
  for_fit[,center_x:=bin_centers[as.numeric(gr_var)]]
  if (var_sigma){
    # get predictions
    resid_at_boundaries <- for_fit[,
                                   get_boundary_preds(gr_var, for_fit, space, reassign_range, gam_ctrl, poly_deg, angle_diff_fun), by = .(gr_var)]
    resid_at_boundaries[,resid:=angle_diff_fun(err, pred)]
    resid_at_boundaries <- dcast(resid_at_boundaries, x_var+dc_var+err~gr_var,
                                 value.var = 'resid')
    resid_at_boundaries <- resid_at_boundaries[,gr_var:=names(.SD)[max.col(-abs(.SD))], .SDcols = 4:5]
    for_fit[resid_at_boundaries, `:=` (gr_var = i.gr_var), on = .(x_var, dc_var, err)]
    for_fit[,center_x:=bin_centers[as.numeric(gr_var)]]
    for_fit[,dist_to_bin_centre:=angle_diff_fun(x_var, center_x)]
    for_fit[,x_var:=center_x+dist_to_bin_centre]
    for (cg in unique(for_fit$gr_var)){
      cur_df <- for_fit[gr_var==cg,.(err, x_var, dist_to_bin_centre, dc_var, outlier, dist_to_card)]
      fit <- gamlss(err~poly(dist_to_bin_centre, poly_deg),
                    ~abs(dist_to_card),
                    data = cur_df,
                    weights = 1-as.numeric(cur_df$outlier),
                    control = gam_ctrl)
      if (debug){
        print(coef(fit))
      }

      for_fit[gr_var==cg, pred:=predict(fit, type='response')]
      for_fit[gr_var==cg, pred_sigma:=predict(fit, what = 'sigma', type='response')]

      for_fit[gr_var==cg, bias:=err*sign(pred)]

      if (debug){
        print(ggplot(for_fit[gr_var==cg], aes(x = dist_to_bin_centre, y = err))+geom_point()+geom_line(aes(y = pred)))
      }
      for_fit[gr_var==cg, c('coef_sigma_int','coef_sigma_slope'):= data.frame(t(coef(fit, what = 'sigma')))]
    }
  } else {
    for_fit[,pred:=predict(rlm(err~poly(x_var,poly_deg),.SD[outlier==F]),newdata=.SD[,.(x_var)]),  by=.(card_groups)]
  }

  for_fit[,pred_lin:=rlm(err~x_var)$fitted.values, by=.(gr_var)]

  for_fit[,be_c:=err-pred]
  for_fit[,which_bin:=as.numeric(gr_var)]
  for_fit[,center_y:=predict(rlm(err~x_var),
                             newdata = data.frame(x_var = center_x)), by=.(gr_var)]

  if (var_sigma){
    for_fit[,outlier:=abs(be_c)>3*pred_sigma]
  }
  else {
    for_fit[,outlier:=abs(be_c)>(3*circ_sd_fun(be_c))]
  }

  # for_fit[outlier==F,pred:=gls(err~x2, weights = varFunc(~sqrt(abs(dist_to_card))), data = .SD)$fitted, by=.(card_groups)]
  # for_fit[outlier==F,weights:=attr(gls(err~x2, weights = varFunc(~sqrt(abs(dist_to_card))), data = .SD)$model$varStruct, 'weights'), by=.(card_groups)]
  #for_fit[outlier==F,pred_rlm:=rlm(err~poly(x2,1), data = .SD)$fitted, by=.(card_groups)]
  if (do_plot){
    for_fit[,outlier_f := factor(ifelse(outlier, 'Outlier','Non-outlier'))]
    sd_val <- data[,circ_sd_fun(err)]
    make_plots_of_biases(for_fit, poly_deg, sd_val)
  }
  for_fit[,.(is_outlier = as.numeric(outlier), pred, be_c, which_bin, bias, bias_type, pred_lin, pred_sigma, coef_sigma_int, coef_sigma_slope)]
}

make_plots_of_biases <- function(data, poly_deg, sd_val){
  require(patchwork)
  require(ggplot2)

  common_plot_pars <- list(scale_x_continuous(breaks = seq(-180,360,90)),
                           labs(x = 'Orientation', shape = NULL, color = 'Bin'),
                           theme(legend.position = c(1,1), legend.justification = c(1,1)),
                          guides(color = guide_legend(override.aes = list(size = 1))))
  alpha = 100/data[,.N]*length(unique(data$gr_var))

  p1<-ggplot(data=data[outlier==F],aes(x = x_var, color = gr_var))+
    geom_point(data = data, aes( y = err, shape = outlier_f), alpha = alpha)+
    geom_line(aes(y = pred), size = 1)+
    #geom_line(aes(y = pred_lin), linetype = 2, size = 1)+
    geom_line(aes(y=pred+3*pred_sigma))+
    geom_line(aes(y=pred-3*pred_sigma))+
    geom_hline(yintercept = c(-1,1)*data[,3*sd_val], linetype = 2)+
    common_plot_pars+
    labs(y = 'Error')

  p1a<-ggplot(data=data[outlier==F],aes(x = x_var, color = gr_var))+
    geom_point(data = data, aes( y = be_c, shape = outlier_f), alpha = alpha)+
    geom_hline(yintercept = c(-1,1)*data[,3*circ_sd_180(be_c)], linetype = 2)+
    geom_line(data=data[outlier==F], aes(y=3*pred_sigma))+
    geom_line(data=data[outlier==F], aes(y=-3*pred_sigma))+
    common_plot_pars+labs(y = 'Bias-corrected error')

  p1b<-ggplot(data, aes(x = x_var, y = bias,  color = gr_var, shape = outlier_f))+
    geom_point(alpha = alpha)+
    common_plot_pars+
    geom_hline(yintercept = c(-1,1)*data[,3*sd_val], linetype = 2)+theme(legend.position = 'none')+labs(y = 'Bias')

  p1<-p1+labs(shape = NULL)+guides(color = 'none')
  print(p1+p1a+p1b+plot_layout(guides = 'collect'))
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
#' @examples
#'
#' require(data.table)
#' dt <- data.table(x = runif(1000,-90,90), y = rnorm(1000))
#' pad_circ(dt, 'x', verbose = TRUE)
#'
pad_circ <- function(data, circ_var, circ_borders=c(-90,90), circ_part = 1/6, verbose = FALSE){
  require(data.table)
  circ_range <- max(circ_borders)-min(circ_borders)

  data1<-copy(data[get(circ_var)<(circ_borders[1]+circ_range*circ_part),])
  data1[,(circ_var):=get(circ_var)+ circ_range]

  data2<-copy(data[get(circ_var)>(circ_borders[2]-circ_range*circ_part),])
  data2[,(circ_var):=get(circ_var)- circ_range]
  if (verbose)
    print(sprintf('Rows in original DT: %i, padded on the left: %i, padded on the right: %i', data[,.N],data1[,.N], data2[,.N]))

  rbind(data,data1,data2)
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
#' @keywords internal
#'
get_boundary_preds <- function(group, data, space, reassign_range, gam_ctrl, poly_deg, angle_diff_fun){
  cur_df <- data[gr_var==group&outlier==F,.(err, x_var, dc_var,
                                            dist_to_bin_centre = angle_diff_fun(x_var, center_x), weight = 1-as.numeric(outlier), adc = abs(dist_to_card), center_x)]

  fit <- gamlss(err~poly(dist_to_bin_centre, poly_deg),
                ~ adc,
                data = cur_df,
                control = gam_ctrl)
  # boundary predictions
  bin_range <- as.numeric(space)/2
  boundary1 <- cur_df$center_x[[1]] - bin_range/2
  boundary2 <- cur_df$center_x[[1]]  + bin_range/2
  data_at_boundaries <- data[outlier == F & at_the_boundary == T,
        .(err, x_var, dc_var,
          dist_to_bin_centre = angle_diff_fun(x_var, cur_df$center_x[1]),
          weight = 1 - as.numeric(outlier), adc = abs(dist_to_card), center_x)]

  data_at_boundaries <- unique(data_at_boundaries)

  pred_at_boundaries <- predict(fit, type ='response',
                                newdata = data_at_boundaries,
                                data = cur_df)
  resid_at_boundaries <- resid(fit)
  data_at_boundaries[, .(x_var, dc_var, dist_to_bin_centre, err, pred = pred_at_boundaries)]
}


a_fun <- function(x) {
  besselI(x, 1, expon.scaled = T)/besselI(x, 0, expon.scaled = T)
}

inverse <- function (f, lower = 1e-16, upper = 1000) {
  function (y) uniroot((function (x) f(x) - y), lower = lower, upper = upper, extendInt = 'yes')[[1]]
}

vm_circ_sd <- function(kappa) {
  sqrt(-2*log(a_fun(kappa)))
}

vm_circ_sd_inverse_deg <-  function(sd_deg) {
  vm_circ_sd_inverse <- inverse(vm_circ_sd)
  sapply(sd_deg, function(x) tryCatch(vm_circ_sd_inverse(x/180*pi), error = function(e) paste("Can't convert sigma = ",x, " to kappa, error ",e)))
}

vm_circ_sd_deg <- function(kappa) {
  sqrt(-2*log(a_fun(kappa)))/pi*180
}
