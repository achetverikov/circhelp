# tmp_test_new_rcb.R
# 4-way engine comparison: gamlss | gamlss_default_fun | gamlss_reml_fun | gamlss_tmb_fun

suppressPackageStartupMessages({
  library(circhelp)
  library(data.table)
  library(gamlss)
  library(gamlss.dist)
})
source("tmp_new_rcb.R")

# -------------------------------------------------------------------
# Preconditions
# -------------------------------------------------------------------

for (fn in c("gamlss_default_fun", "gamlss_fast_fun", "gamlss_reml_fun",
             "gamlss_tmb_fun", "remove_cardinal_biases")) {
  if (!exists(fn)) stop(fn, " not found — source tmp_new_rcb.R first.")
}

# -------------------------------------------------------------------
# TMB dynamic library check
# -------------------------------------------------------------------

.tmb_dll_path    <- "src/normal_ls.so"
.tmb_dll_name    <- "normal_ls"          # name R registers after dyn.load()

.tmb_available <- tryCatch({
  if (!requireNamespace("TMB", quietly = TRUE)) stop("TMB not installed")
  if (!file.exists(.tmb_dll_path))              stop("DLL file not found")

  # Load only if not already registered (avoids "already loaded" warnings)
  loaded_names <- names(getLoadedDLLs())
  if (!.tmb_dll_name %in% loaded_names) {
    dyn.load(.tmb_dll_path)
  }

  # Confirm TMB can now see the DLL by its registered name
  .tmb_dll_name %in% names(getLoadedDLLs())
}, error = function(e) FALSE)

if (.tmb_available) {
  cat("TMB DLL available as '", .tmb_dll_name, "' (", .tmb_dll_path, ")\n", sep = "")

  # Smoke-test: fit a tiny model to confirm the template compiles & runs correctly
  .tmb_smoke <- tryCatch({
    set.seed(0)
    n   <- 20L
    X1  <- cbind(1, rnorm(n))
    X2  <- cbind(1, abs(rnorm(n)))
    y0  <- X1 %*% c(0, 1) + rnorm(n, 0, 0.5)
    obj <- TMB::MakeADFun(
      data       = list(y = y0, X_mu = X1, X_sigma = X2,
                        weights = rep(1, n),
                        D_mu = matrix(0, 0, 2), lambda_mu = 0),
      parameters = list(beta_mu = c(0, 0), beta_sigma = c(0, 0)),
      DLL        = .tmb_dll_name,
      silent     = TRUE
    )
    opt <- stats::nlminb(obj$par, obj$fn, obj$gr)
    stopifnot(opt$convergence == 0)
    cat("TMB smoke-test passed (", .tmb_dll_name, ")\n", sep = "")
  }, error = function(e) {
    stop("TMB smoke-test FAILED for '", .tmb_dll_name, "':\n  ", conditionMessage(e),
         "\nRecompile with: TMB::compile('src/normal_ls.cpp'); dyn.load(TMB::dynlib('src/normal_ls'))",
         call. = FALSE)
  })
} else {
  warning(
    "TMB engine will be skipped: '", .tmb_dll_path, "' not found or failed to load.\n",
    "  Compile with: TMB::compile('src/normal_ls.cpp'); dyn.load(TMB::dynlib('src/normal_ls'))",
    call. = FALSE
  )
}

# -------------------------------------------------------------------
# Engines
# -------------------------------------------------------------------

ENGINES <- c(
  list(
    gamlss  = gamlss::gamlss,
    default = gamlss_default_fun,
    fast    = gamlss_fast_fun,
    reml    = gamlss_reml_fun
  ),
  if (.tmb_available) list(tmb = gamlss_tmb_fun) else list()
)

# -------------------------------------------------------------------
# Helpers
# -------------------------------------------------------------------

compare_two <- function(ref, cmp, label = "", engine = "") {
  ref <- as.data.table(ref)
  cmp <- as.data.table(cmp)
  stopifnot(nrow(ref) == nrow(cmp))

  bin_diff     <- ref$which_bin != cmp$which_bin
  bin_diff[is.na(bin_diff)] <- FALSE
  outlier_diff <- ref$is_outlier != cmp$is_outlier
  outlier_diff[is.na(outlier_diff)] <- FALSE

  data.table(
    dataset          = label,
    engine           = engine,
    n                = nrow(ref),
    same_bins        = !any(bin_diff),
    n_bin_diff       = sum(bin_diff),
    same_outliers    = !any(outlier_diff),
    n_outlier_diff   = sum(outlier_diff),
    max_abs_pred     = max(abs(ref$pred        - cmp$pred),        na.rm = TRUE),
    mean_abs_pred    = mean(abs(ref$pred       - cmp$pred),        na.rm = TRUE),
    max_abs_sigma    = max(abs(ref$pred_sigma  - cmp$pred_sigma),  na.rm = TRUE),
    mean_abs_sigma   = mean(abs(ref$pred_sigma - cmp$pred_sigma),  na.rm = TRUE),
    max_abs_be_c     = max(abs(ref$be_c        - cmp$be_c),        na.rm = TRUE),
    llik_ref         = unique(ref$total_log_lik)[1],
    llik_cmp         = unique(cmp$total_log_lik)[1],
    llik_diff        = unique(cmp$total_log_lik)[1] - unique(ref$total_log_lik)[1]
  )
}

# -------------------------------------------------------------------
# Core 4-way runner
# -------------------------------------------------------------------

run_4way <- function(data,
                     err_col,
                     x_col,
                     label,
                     space               = "180",
                     init_outliers       = NULL,
                     bias_type           = "fit",
                     reassign_at_boundaries = TRUE,
                     reassign_range      = 2,
                     break_points        = NULL,
                     poly_deg            = 4,
                     var_sigma_poly_deg  = 4,
                     debug               = FALSE) {
  data <- as.data.table(data)

  # Drop rows with NA in the response or predictor
  keep <- !is.na(data[[err_col]]) & !is.na(data[[x_col]])
  if (any(!keep)) {
    message("  [", label, "] dropping ", sum(!keep), " row(s) with NA in ",
            err_col, " or ", x_col)
    data <- data[keep]
  }

  err  <- data[[err_col]]
  x    <- data[[x_col]]

  cat("\n============================================================\n")
  cat("Dataset:", label, "  N:", length(err),
      "  space:", space, "  bias_type:", bias_type, "\n")
  cat("============================================================\n")

  fits  <- vector("list", length(ENGINES))
  times <- numeric(length(ENGINES))
  names(fits) <- names(times) <- names(ENGINES)

  for (eng in names(ENGINES)) {
    cat("  Running", eng, "... ")
    t <- system.time(
      fits[[eng]] <- tryCatch(
        remove_cardinal_biases(
          err                    = err,
          x                      = x,
          space                  = space,
          bias_type              = bias_type,
          plots                  = "hide",
          poly_deg               = poly_deg,
          var_sigma              = TRUE,
          var_sigma_poly_deg     = var_sigma_poly_deg,
          reassign_at_boundaries = reassign_at_boundaries,
          reassign_range         = reassign_range,
          break_points           = break_points,
          init_outliers          = init_outliers,
          debug                  = debug,
          gamlss_fun             = ENGINES[[eng]]
        ),
        error = function(e) { message("ERROR: ", conditionMessage(e)); NULL }
      )
    )
    times[[eng]] <- unname(t["elapsed"])
    cat(round(times[[eng]], 2), "s\n")
  }

  # Comparisons vs gamlss baseline
  ref_fit <- as.data.table(fits[["gamlss"]])

  cmp_rows <- rbindlist(lapply(setdiff(names(ENGINES), "gamlss"), function(eng) {
    if (is.null(fits[[eng]])) {
      return(data.table(dataset = label, engine = eng, n = nrow(ref_fit),
                        error = TRUE))
    }
    row <- compare_two(ref_fit, as.data.table(fits[[eng]]),
                       label = label, engine = eng)
    row[, time_gamlss  := times[["gamlss"]]]
    row[, time_engine  := times[[eng]]]
    row[, speedup      := times[["gamlss"]] / times[[eng]]]
    row
  }), fill = TRUE)

  # Print concise table
  print(cmp_rows[, .(engine, n_bin_diff, n_outlier_diff,
                      max_abs_pred, max_abs_sigma,
                      llik_diff, time_gamlss, time_engine, speedup)])

  invisible(list(
    label  = label,
    fits   = fits,
    times  = times,
    cmp    = cmp_rows
  ))
}

# -------------------------------------------------------------------
# Datasets
# -------------------------------------------------------------------

data("Pascucci_et_al_2019_data", package = "circhelp")
data("Bae_Luck_2018_data",       package = "circhelp")

pascucci_dt <- as.data.table(Pascucci_et_al_2019_data)
bae_dt      <- as.data.table(Bae_Luck_2018_data)

# -------------------------------------------------------------------
# Run 4-way comparisons — all observers / subjects per dataset
# -------------------------------------------------------------------

pascucci_obs <- sort(unique(pascucci_dt$observer))
bae_subj     <- sort(unique(bae_dt$subject_Num))

# Dataset 1: Pascucci et al. 2019 — orientation (180-degree space)
res_pascucci <- lapply(pascucci_obs, function(obs) {
  run_4way(
    data       = pascucci_dt[observer == obs],
    err_col    = "err",
    x_col      = "orientation",
    label      = paste0("Pascucci_obs_", obs),
    space      = "180",
    bias_type  = "fit",
    reassign_at_boundaries = TRUE,
    reassign_range         = 2
  )
})

# Dataset 2: Bae & Luck 2018 — motion direction (360-degree space)
res_bae <- lapply(bae_subj, function(subj) {
  run_4way(
    data       = bae_dt[subject_Num == subj],
    err_col    = "err",
    x_col      = "TargetDirection",
    label      = paste0("Bae_subj_", subj),
    space      = "360",
    bias_type  = "fit",
    reassign_at_boundaries = TRUE,
    reassign_range         = 2
  )
})

all_results <- c(res_pascucci, res_bae)

# -------------------------------------------------------------------
# Combined summary
# -------------------------------------------------------------------

all_cmp <- rbindlist(lapply(all_results, `[[`, "cmp"), fill = TRUE)

cat("\n============================================================\n")
cat("Per-dataset comparison vs gamlss baseline\n")
cat("============================================================\n")
print(all_cmp[, .(
  dataset, engine,
  n_bin_diff, n_outlier_diff,
  max_abs_pred, mean_abs_pred,
  max_abs_sigma, mean_abs_sigma,
  llik_diff,
  time_gamlss, time_engine, speedup
)])

# -------------------------------------------------------------------
# Divergence summary: aggregated across datasets per engine
# -------------------------------------------------------------------

cat("\n============================================================\n")
cat("Divergence summary (aggregated across all datasets)\n")
cat("  all metrics vs gamlss baseline\n")
cat("============================================================\n")

div_summary <- all_cmp[, .(
  n_datasets        = .N,

  # Bin / outlier assignment agreement
  pct_any_bin_diff  = mean(n_bin_diff > 0) * 100,
  mean_n_bin_diff   = mean(n_bin_diff),
  max_n_bin_diff    = max(n_bin_diff),

  pct_any_out_diff  = mean(n_outlier_diff > 0) * 100,
  mean_n_out_diff   = mean(n_outlier_diff),

  # Prediction accuracy vs gamlss
  mean_max_abs_pred  = mean(max_abs_pred,  na.rm = TRUE),
  mean_mean_abs_pred = mean(mean_abs_pred, na.rm = TRUE),
  worst_max_abs_pred = max(max_abs_pred,   na.rm = TRUE),

  # Sigma accuracy vs gamlss
  mean_max_abs_sigma  = mean(max_abs_sigma,  na.rm = TRUE),
  mean_mean_abs_sigma = mean(mean_abs_sigma, na.rm = TRUE),
  worst_max_abs_sigma = max(max_abs_sigma,   na.rm = TRUE),

  # Log-likelihood difference vs gamlss (positive = engine beats gamlss)
  mean_llik_diff = mean(llik_diff, na.rm = TRUE),
  min_llik_diff  = min(llik_diff,  na.rm = TRUE),
  max_llik_diff  = max(llik_diff,  na.rm = TRUE),

  # Speed
  mean_speedup  = mean(speedup, na.rm = TRUE),
  worst_speedup = min(speedup,  na.rm = TRUE)
), by = engine]

print(div_summary)

# One table: metrics as rows, engines as columns, values averaged across observers
metric_avg <- melt(
  all_cmp[, .(engine, max_abs_pred, mean_abs_pred,
               max_abs_sigma, mean_abs_sigma,
               n_bin_diff, n_outlier_diff,
               llik_diff, speedup)],
  id.vars     = "engine",
  variable.name = "metric"
)[, .(mean_value = mean(value, na.rm = TRUE)), by = .(metric, engine)]

cat("\n--- Mean across observers (columns = engines) ---\n")
print(dcast(metric_avg, metric ~ engine, value.var = "mean_value"))

# -------------------------------------------------------------------
# Timing summary across engines
# -------------------------------------------------------------------

cat("\n============================================================\n")
cat("Timing summary (seconds)\n")
cat("============================================================\n")

timing_dt <- rbindlist(lapply(all_results, function(res) {
  rbindlist(lapply(names(res$times), function(eng) {
    data.table(dataset = res$label, engine = eng,
               elapsed = res$times[[eng]])
  }))
}))

timing_wide <- dcast(timing_dt, dataset ~ engine, value.var = "elapsed")
print(timing_wide)

cat("\nMean elapsed per engine:\n")
print(timing_dt[, .(mean_s = mean(elapsed), total_s = sum(elapsed)), by = engine])

# -------------------------------------------------------------------
# Save
# -------------------------------------------------------------------

saveRDS(
  list(all_cmp = all_cmp, timing = timing_dt, all_results = all_results),
  file = "remove_cardinal_biases_4way_comparison.rds"
)
data.table::fwrite(all_cmp, file = "remove_cardinal_biases_4way_summary.csv")

cat("\nSaved:\n")
cat("  remove_cardinal_biases_4way_comparison.rds\n")
cat("  remove_cardinal_biases_4way_summary.csv\n")
