## ----include = FALSE----------------------------------------------------------
library(patchwork)
library(circhelp)
library(data.table)
library(ggplot2)

if (requireNamespace("ragg", quietly = TRUE)) {
  knitr::opts_chunk$set(dev = "ragg_png")
}

fig_width <- 3.6
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  warning = FALSE,
  fig.width = fig_width, fig.height = fig_width, dpi = 320, res = 320,
  out.width = "40%", out.height = "40%"
)

## ----setup--------------------------------------------------------------------
# load the data
# data <- fread('https://zenodo.org/record/2544946/files/Experiment2_rawdata.csv?download=1')
data <- Pascucci_et_al_2019_data
data[, err := angle_diff_180(reported, orientation)] # response errors
data[, prev_ori := shift(orientation), by = observer] # orientation on previous trial
data[, diff_in_ori := angle_diff_180(prev_ori, orientation)] # shift in orientations between trials

## ----serial_dependence--------------------------------------------------------
ggplot(pad_circ(data, "diff_in_ori"), aes(x = diff_in_ori, y = err)) +
  geom_line(aes(group = observer), stat = "smooth", size = 0.4, color = "black", alpha = 0.2, method = "loess") +
  geom_smooth(se = T, method = "loess") +
  coord_cartesian(xlim = c(-90, 90)) +
  scale_x_continuous(breaks = seq(-90, 90, 90)) +
  labs(y = "Error, °", x = "Orientation difference, °")

## ----card_biases_ex-----------------------------------------------------------
ggplot(data[observer == 4, ], aes(x = angle_diff_180(orientation, 0), y = err)) +
  geom_point() +
  coord_cartesian(xlim = c(-90, 90)) +
  scale_x_continuous(breaks = seq(-90, 90, 90)) +
  labs(y = "Error, °", x = "Orientation, °")

## ----correct_biases_ex, fig.width = fig_width*3, out.width='100%'-------------
ex_subj_data <- data[observer == 4, ]
res <- remove_cardinal_biases(ex_subj_data$err, ex_subj_data$orientation, plots = "show")

## -----------------------------------------------------------------------------
data[, err_corrected := remove_cardinal_biases(err, orientation)$be_c, by = observer]
data[, err_mean_corrected := angle_diff_180(err, circ_mean_180(err)), by = observer]

## ----fig.width = fig_width*3, out.width='100%'--------------------------------
datam <- melt(data[!is.na(diff_in_ori)], id.vars = c("diff_in_ori", "observer"), measure.vars = c("err", "err_corrected", "err_mean_corrected"))
datam[, variablef := factor(variable, levels = c("err", "err_mean_corrected", "err_corrected"), labels = c("Raw error", "Mean-corrected", "Mean and cardinal bias removed"))]
datam[, err_rel_to_prev_targ := ifelse(diff_in_ori < 0, -value, value)]

ggplot(pad_circ(datam, "diff_in_ori"), aes(x = diff_in_ori, y = value)) +
  geom_line(aes(group = observer), stat = "smooth", size = 0.4, color = "black", alpha = 0.2, method = "loess") +
  geom_smooth(se = TRUE, method = "loess") +
  facet_grid(~variablef) +
  coord_cartesian(xlim = c(-90, 90), ylim = c(-5, 5)) +
  scale_x_continuous(breaks = seq(-90, 90, 90)) +
  labs(y = "Error, °", x = "Orientation difference, °")

## ----fig.width = fig_width*3, out.width='100%'--------------------------------
ggplot(datam, aes(
  x = abs(diff_in_ori),
  y = err_rel_to_prev_targ,
  color = variablef
)) +
  geom_hline(linetype = 2, yintercept = 0) +
  geom_smooth(se = TRUE, method = "loess") +
  facet_grid(~variablef) +
  coord_cartesian(xlim = c(0, 90)) +
  theme(legend.position = "top") +
  labs(color = NULL, y = "Bias towards previous targets", x = "Absolute orientation difference, °") +
  scale_x_continuous(breaks = seq(0, 90, 30))

## ----density, fig.width=fig_width*2, out.width='66%'--------------------------
data[, abs_diff_in_ori := abs(diff_in_ori)]
data[, err_rel_to_prev_targ := ifelse(diff_in_ori < 0, -err_corrected, err_corrected)]

err_dens <- density_asymmetry(data[!is.na(err_rel_to_prev_targ)], circ_space = 180, weights_sd = 10, 
                              xvar = "abs_diff_in_ori", yvar = "err_rel_to_prev_targ", by = c("observer"))

data[,abs_diff_in_ori_bin:=apastats::adv_cut(abs_diff_in_ori, cuts = seq(0,90,10), include_oob = T, num_labels = F, labels_at_means = T)]

mean_err <- apastats::get_superb_ci(data[!is.na(err_rel_to_prev_targ)], wid = 'observer', within = 'abs_diff_in_ori_bin', value_var = 'err_rel_to_prev_targ')
setDT(mean_err)
mean_err[,abs_diff_in_ori_bin:=as.numeric(abs_diff_in_ori_bin)]

((ggplot(err_dens[, mean_cl_normal(delta * 100), by = .(dist)], aes(x = dist, y = y, ymin = ymin, ymax = ymax)) +
  geom_line() +
  geom_ribbon(alpha = 0.1) +
  labs(y = "Asymmetry in error probability density, %")) +
 (ggplot(mean_err, 
         aes(x = abs_diff_in_ori_bin, y = center, ymin = center+lowerwidth, ymax = center+upperwidth)) +
  geom_line()+
  geom_pointrange() +
  labs(y = "Mean error, deg."))) & labs(x = "Absolute orientation difference, °")

