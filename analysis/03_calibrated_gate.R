library(here)
library(lavaan)
source(here("R", "measures.R"))
source(here("R", "simulate.R"))
source(here("R", "parcels.R"))
source(here("R", "gate_harness.R"))

n_subjects <- 322
configs3 <- list(c(200, 200), c(200, 100), c(100, 100))
cfg_name <- function(cfg) paste0(cfg[1], "/", cfg[2])
seed_eval_main <- 20260711
seed_eval_spot <- 777
r_cal <- 20
r_eval <- 30
r_spot <- 15

dir.create(here("figures"), showWarnings = FALSE)

message("=== CALIBRATION (SIMULATION, seeds 1000+i) ===")
cal_grid <- seq(0, 0.6, by = 0.1)
cal_rows <- list()
counter <- 0
for (cfg in configs3) {
  message("calibration config ", cfg_name(cfg), " ...")
  for (rho_s in cal_grid) {
    for (r in seq_len(r_cal)) {
      counter <- counter + 1
      row <- gate_replicate(
        n_subjects, list(level = 0.4, sens = rho_s, dprime = 0.3),
        cfg[1], cfg[2], 1000 + counter
      )
      row$population_rho_sens <- rho_s
      row$config <- cfg_name(cfg)
      cal_rows[[counter]] <- row
    }
  }
}
cal <- do.call(rbind, cal_rows)
cal$phase <- "calibration"

calibrations <- lapply(configs3, function(cfg) {
  d <- cal[cal$config == cfg_name(cfg) & cal$converged & cal$admissible, ]
  lin <- lm(rho_sens ~ realized_sens, data = d)
  quad <- lm(rho_sens ~ realized_sens + I(realized_sens^2), data = d)
  qs <- summary(quad)$coefficients
  list(
    config = cfg_name(cfg),
    n_fit = nrow(d),
    a = unname(coef(lin)[1]),
    b = unname(coef(lin)[2]),
    V = vcov(lin),
    resid_sd = summary(lin)$sigma,
    quad_coef = qs["I(realized_sens^2)", "Estimate"],
    quad_p = qs["I(realized_sens^2)", "Pr(>|t|)"],
    lin_data = d[, c("realized_sens", "rho_sens")]
  )
})
names(calibrations) <- vapply(configs3, cfg_name, "")

cal_coefs <- do.call(rbind, lapply(calibrations, function(cl) data.frame(
  data_type = "SIMULATION", config = cl$config, n_fit = cl$n_fit,
  intercept_a = cl$a, slope_b = cl$b,
  se_a = sqrt(cl$V[1, 1]), se_b = sqrt(cl$V[2, 2]),
  resid_sd = cl$resid_sd, quad_coef = cl$quad_coef, quad_p = cl$quad_p
)))
message("\nCalibration mappings (est = a + b * truth):")
print(cal_coefs, digits = 3, row.names = FALSE)

png(here("figures", "step0b2_SIMULATION_calibration_residuals.png"),
    width = 1200, height = 420)
par(mfrow = c(1, 3), mar = c(4.5, 4.5, 3, 1))
for (cl in calibrations) {
  d <- cl$lin_data
  res <- d$rho_sens - (cl$a + cl$b * d$realized_sens)
  plot(d$realized_sens, res,
       xlab = "realized planted rho_sens", ylab = "linear-fit residual",
       main = paste0("SIMULATION calibration, ", cl$config,
                     " (quad p = ", signif(cl$quad_p, 2), ")"),
       pch = 16, col = "#00000055")
  abline(h = 0, lty = 2)
  lo <- lowess(d$realized_sens, res)
  lines(lo, col = "red", lwd = 2)
}
dev.off()

correct_sens <- function(rs, rs_se, cl) {
  est <- (rs - cl$a) / cl$b
  grad <- c(-1 / cl$b, -est / cl$b)
  var_cal <- as.numeric(t(grad) %*% cl$V %*% grad)
  se <- sqrt(rs_se^2 / cl$b^2 + var_cal)
  list(est = est, se = se)
}

apply_correction <- function(row, cl) {
  cs <- correct_sens(row$rho_sens, row$rho_sens_se, cl)
  gap_c <- row$rho_level - cs$est
  gap_var <- row$rho_level_se^2 + cs$se^2 - 2 * row$cov_rs_rl / cl$b
  gap_se <- sqrt(pmax(gap_var, 0))
  row$rho_sens_c <- cs$est
  row$rho_sens_c_se <- cs$se
  row$rho_sens_c_lo <- cs$est - 1.96 * cs$se
  row$rho_sens_c_hi <- cs$est + 1.96 * cs$se
  row$gap_c <- gap_c
  row$gap_c_lo <- gap_c - 1.96 * gap_se
  row$gap_c_hi <- gap_c + 1.96 * gap_se
  row
}

scenarios <- list(
  GAP = list(level = 0.6, sens = 0.2, dprime = 0),
  NULLGAP = list(level = 0.4, sens = 0.4, dprime = 0),
  ZEROSENS = list(level = 0.4, sens = 0.0, dprime = 0),
  DCONF = list(level = 0.0, sens = 0.0, dprime = 0.4)
)

run_eval_suite <- function(base_seed, n_rep, label) {
  rows <- list()
  idx <- 0
  for (sc in names(scenarios)) {
    for (cfg in configs3) {
      message(sprintf("[%s] %s %s ...", label, sc, cfg_name(cfg)))
      cl <- calibrations[[cfg_name(cfg)]]
      for (r in seq_len(n_rep)) {
        idx <- idx + 1
        row <- gate_replicate(n_subjects, scenarios[[sc]], cfg[1], cfg[2],
                              base_seed + idx * 1009)
        row$scenario <- sc
        row$config <- cfg_name(cfg)
        row$rep <- r
        rows[[idx]] <- apply_correction(row, cl)
      }
    }
  }
  out <- do.call(rbind, rows)
  out$suite <- label
  out
}

eval_main <- run_eval_suite(seed_eval_main, r_eval, "eval_main")
eval_spot <- run_eval_suite(seed_eval_spot, r_spot, "eval_spot")

summarize_cell <- function(d) {
  ok <- d$converged & d$admissible
  realized_gap <- d$realized_level - d$realized_sens
  width_unc <- d$rho_sens_hi - d$rho_sens_lo
  width_cor <- d$rho_sens_c_hi - d$rho_sens_c_lo
  data.frame(
    data_type = "SIMULATION", scenario = d$scenario[1], config = d$config[1],
    n_rep = nrow(d), pct_bad = 100 * mean(!ok),
    med_rs_c = median(d$rho_sens_c[ok]),
    med_rs_c_bias = median((d$rho_sens_c - d$realized_sens)[ok]),
    rs_c_coverage = sum(ok & d$rho_sens_c_lo <= d$realized_sens &
                          d$realized_sens <= d$rho_sens_c_hi),
    med_gap_c = median(d$gap_c[ok]),
    med_gap_c_bias = median((d$gap_c - realized_gap)[ok]),
    gap_c_positive = sum(ok & d$gap_c > 0),
    gap_c_ci_zero = sum(ok & d$gap_c_lo <= 0 & 0 <= d$gap_c_hi),
    med_rs_u_bias = median((d$rho_sens - d$realized_sens)[ok]),
    med_gap_u_bias = median((d$gap - realized_gap)[ok]),
    med_width_ratio = median((width_cor / width_unc)[ok]),
    n_narrower = sum(ok & width_cor < width_unc)
  )
}

summarize_suite <- function(res) {
  cells <- split(res, list(res$scenario, res$config), drop = TRUE)
  out <- do.call(rbind, lapply(cells, summarize_cell))
  rownames(out) <- NULL
  out[order(out$config, out$scenario), ]
}

sum_main <- summarize_suite(eval_main)
sum_spot <- summarize_suite(eval_spot)

message("\n=== EVALUATION SUMMARY (SIMULATION, fresh seed ", seed_eval_main,
        ") ===")
print(sum_main, digits = 3, row.names = FALSE)

config_criteria <- function(smry, n_rep) {
  lapply(setNames(nm = unique(smry$config)), function(cf) {
    g <- function(sc) smry[smry$scenario == sc & smry$config == cf, ]
    gap <- g("GAP"); nul <- g("NULLGAP"); zer <- g("ZEROSENS"); dcf <- g("DCONF")
    all_cells <- rbind(gap, nul, zer, dcf)
    out <- list(
      config = cf,
      a_gap_bias = abs(gap$med_rs_c_bias) <= 0.07,
      a_gap_coverage = gap$rs_c_coverage >= ceiling(25 / 30 * n_rep),
      a_gap_positive = gap$gap_c_positive >= ceiling(27 / 30 * n_rep),
      b_null_gap_bias = abs(nul$med_gap_c_bias) <= 0.07,
      b_null_gap_ci_zero = nul$gap_c_ci_zero >= ceiling(25 / 30 * n_rep),
      c_zerosens_median = abs(zer$med_rs_c) <= 0.07,
      d_dconf_spurious = dcf$med_rs_c < 0.10,
      e_convergence = all(all_cells$pct_bad <= 20),
      f_ci_honesty = all(all_cells$med_width_ratio >= 1) &
        sum(all_cells$n_narrower) == 0
    )
    out$all_pass <- all(unlist(out[setdiff(names(out), "config")]))
    out
  })
}

crit_main <- config_criteria(sum_main, r_eval)
crit_spot <- config_criteria(sum_spot, r_spot)

pass <- vapply(crit_main, `[[`, TRUE, "all_pass")
branch <- if (!pass[["200/200"]]) {
  "GATE FAILED OUTRIGHT: revert to 0b'-validated UNCORRECTED estimator, 200/200-only design + floor reframe"
} else if (pass[["200/100"]] && pass[["100/100"]]) {
  "GATE PASSED, FULL: all 6 domain pairs confirmatory, corrected d'-partialed latent estimator"
} else if (pass[["200/100"]]) {
  "GATE PASSED WITH SCOPE RESTRICTION: 5 pairs confirmatory; Calories x GDP (100/100) bias-band only"
} else {
  "GATE PASSED, NARROW: memory x vision (200/200) sole confirmatory pair; other pairs descriptive with bias bands"
}
spot_pass <- vapply(crit_spot, `[[`, TRUE, "all_pass")

message("\n=== PRE-SPECIFIED CRITERIA (corrected estimator, main eval seed) ===")
str(crit_main)
message("Config pass map (main): ",
        paste(names(pass), pass, sep = "=", collapse = ", "))
message("Config pass map (spot-check seed ", seed_eval_spot, "): ",
        paste(names(spot_pass), spot_pass, sep = "=", collapse = ", "))
message("\nDECISION TREE BRANCH: ", branch)

write.csv(cal, here("results", "step0b2_SIMULATION_calibration_replicates.csv"),
          row.names = FALSE)
write.csv(cal_coefs, here("results", "step0b2_SIMULATION_calibration_coefs.csv"),
          row.names = FALSE)
write.csv(rbind(eval_main, eval_spot),
          here("results", "step0b2_SIMULATION_eval_replicates.csv"),
          row.names = FALSE)
write.csv(rbind(cbind(sum_main, suite = "eval_main"),
                cbind(sum_spot, suite = "eval_spot")),
          here("results", "step0b2_SIMULATION_eval_summary.csv"),
          row.names = FALSE)

meta <- list(
  data_type = "SIMULATION",
  script = "analysis/03_calibrated_gate.R",
  run_utc = format(Sys.time(), "%Y-%m-%d %H:%M:%S", tz = "UTC"),
  r_version = R.version.string,
  lavaan_version = as.character(packageVersion("lavaan")),
  seeds = list(calibration = "1000+i (i = 1..420)",
               eval_main = seed_eval_main, eval_spot = seed_eval_spot),
  n_subjects = n_subjects,
  r_cal = r_cal, r_eval = r_eval, r_spot = r_spot,
  calibration_coefs = lapply(calibrations, function(cl) {
    cl[c("config", "n_fit", "a", "b", "resid_sd", "quad_coef", "quad_p")]
  }),
  criteria_main = crit_main,
  criteria_spot = crit_spot,
  decision_tree_branch = branch
)
jsonlite::write_json(meta, here("results", "step0b2_SIMULATION_meta.json"),
                     auto_unbox = TRUE, pretty = TRUE, digits = 8)

message("\nStep 0b'' complete. Outputs in results/step0b2_SIMULATION_*.")
