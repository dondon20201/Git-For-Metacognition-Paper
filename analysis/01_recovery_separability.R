library(here)
source(here("R", "measures.R"))
source(here("R", "simulate.R"))

dir.create(here("results"), showWarnings = FALSE)

seed_main <- 20260710
seed_check <- 42424242
n_subjects <- 322
n_boot <- 2000

step0a <- read.csv(here("results", "step0a_basic_summary.csv"))
acc_range <- range(step0a$mean_accuracy)
d_prime_mean <- 1.2
message(sprintf(
  "Empirical per-domain accuracy range (step 0a): %.3f-%.3f; simulator d' mean %.1f -> accuracy %.3f",
  acc_range[1], acc_range[2], d_prime_mean, pnorm(d_prime_mean / 2)
))

set.seed(1)
for (i in 1:20) {
  conf <- sample(1:7, 60, replace = TRUE)
  acc <- rbinom(60, 1, 0.6)
  if (sum(acc) %in% c(0, 60)) next
  cc <- conf[acc == 1]
  ce <- conf[acc == 0]
  pairs <- outer(cc, ce, `-`)
  bf_auroc <- (sum(pairs > 0) + 0.5 * sum(pairs == 0)) / length(pairs)
  bf_gamma <- (sum(pairs > 0) - sum(pairs < 0)) /
    (sum(pairs > 0) + sum(pairs < 0))
  stopifnot(
    abs(auroc2(conf, acc) - bf_auroc) < 1e-12,
    abs(gk_gamma(conf, acc) - bf_gamma) < 1e-12
  )
}
message("Estimator unit checks vs brute force: PASS (20 random datasets)")

set.seed(seed_main)
n_cell <- 200000

cell_stats <- function(d_prime, lambda, beta) {
  tr <- simulate_trials(n_cell, d_prime, lambda, beta)
  c(
    mean_conf = mean(tr$confidence),
    accuracy = mean(tr$accuracy),
    auroc2 = auroc2(tr$confidence, tr$accuracy),
    gamma = gk_gamma(tr$confidence, tr$accuracy)
  )
}

beta_grid <- seq(-1.5, 1.5, by = 0.5)
check_level <- data.frame(
  varied = "beta", d_prime = 1.2, lambda = 0.5, beta = beta_grid,
  t(vapply(beta_grid, function(b) cell_stats(1.2, 0.5, b), numeric(4)))
)

lambda_grid <- seq(0.1, 0.9, by = 0.1)
check_sens <- data.frame(
  varied = "lambda", d_prime = 1.2, lambda = lambda_grid, beta = 0,
  t(vapply(lambda_grid, function(l) cell_stats(1.2, l, 0), numeric(4)))
)

dprime_grid <- c(0.8, 1.2, 1.6)
check_acc <- data.frame(
  varied = "d_prime", d_prime = dprime_grid, lambda = 0.5, beta = 0,
  t(vapply(dprime_grid, function(d) cell_stats(d, 0.5, 0), numeric(4)))
)
check_acc$accuracy_target <- pnorm(check_acc$d_prime / 2)

selfcheck <- rbind(
  cbind(check_level, accuracy_target = pnorm(0.6)),
  cbind(check_sens, accuracy_target = pnorm(0.6)),
  check_acc
)
selfcheck$data_type <- "SIMULATION"

selfcheck_flags <- list(
  level_moves_mean_conf_monotonically = all(diff(check_level$mean_conf) > 0),
  level_leaves_accuracy_flat_range = diff(range(check_level$accuracy)),
  level_auroc2_coupling_range = diff(range(check_level$auroc2)),
  sens_moves_auroc2_monotonically = all(diff(check_sens$auroc2) > 0),
  sens_leaves_mean_conf_flat_range = diff(range(check_sens$mean_conf)),
  sens_leaves_accuracy_flat_range = diff(range(check_sens$accuracy)),
  max_abs_accuracy_vs_target = max(abs(selfcheck$accuracy - selfcheck$accuracy_target))
)
message("\n=== SELF-CHECK (SIMULATION) ===")
str(selfcheck_flags)
write.csv(selfcheck, here("results", "step0b_SIMULATION_selfcheck.csv"),
          row.names = FALSE)

recovery_leakage <- function(seed, n_trials) {
  set.seed(seed + n_trials)
  cohort <- run_cohort(draw_traits(n_subjects, d_prime_mean), n_trials)
  row <- function(kind, est_col, planted_col) {
    bc <- boot_cor(cohort[[est_col]], cohort[[planted_col]], n_boot)
    data.frame(
      data_type = "SIMULATION", n_trials = n_trials, kind = kind,
      estimated = est_col, planted = planted_col,
      r = bc["r"], ci_lo = bc["ci_lo"], ci_hi = bc["ci_hi"]
    )
  }
  rbind(
    row("recovery", "all.mean_conf", "beta"),
    row("recovery", "all.auroc2", "lambda"),
    row("recovery", "all.gamma", "lambda"),
    row("leakage", "all.mean_conf", "lambda"),
    row("leakage", "all.auroc2", "beta"),
    row("leakage", "all.gamma", "beta"),
    row("confound_dprime", "all.mean_conf", "d_prime"),
    row("confound_dprime", "all.auroc2", "d_prime"),
    row("confound_dprime", "all.gamma", "d_prime")
  )
}

message("\n=== RECOVERY + LEAKAGE (SIMULATION, main seed) ===")
rec_main <- rbind(
  recovery_leakage(seed_main, 200),
  recovery_leakage(seed_main, 100)
)
print(rec_main, digits = 3, row.names = FALSE)

message("\n=== ASYMPTOTIC CONSISTENCY, 5000 trials (SIMULATION) ===")
rec_asymp <- recovery_leakage(seed_main, 5000)
print(rec_asymp, digits = 3, row.names = FALSE)

correlated_pair <- function(n, rho) {
  z <- matrix(rnorm(n * 2), n) %*% chol(matrix(c(1, rho, rho, 1), 2))
  list(a = z[, 1], b = z[, 2])
}

detectability <- function(seed, scenario, rho_level, rho_sens,
                          n_trials_a, n_trials_b) {
  set.seed(seed + n_trials_a * 7 + n_trials_b * 13 +
             round(1000 * (rho_level + 3 * rho_sens)))
  bl <- correlated_pair(n_subjects, rho_level)
  sl <- correlated_pair(n_subjects, rho_sens)
  mk <- function(z_beta, z_lambda) data.frame(
    d_prime = pmax(rnorm(n_subjects, d_prime_mean, 0.3), 0.2),
    lambda = 0.15 + 0.7 * pnorm(z_lambda),
    beta = 0.8 * z_beta
  )
  ta <- mk(bl$a, sl$a)
  tb <- mk(bl$b, sl$b)
  ca <- run_cohort(ta, n_trials_a)
  cb <- run_cohort(tb, n_trials_b)

  sb <- function(odd_a, even_a) {
    r <- cor(odd_a, even_a)
    2 * r / (1 + r)
  }
  measure_row <- function(component, est_col, planted_col) {
    truth <- cor(ta[[planted_col]], tb[[planted_col]])
    bc <- boot_cor(ca[[paste0("all.", est_col)]], cb[[paste0("all.", est_col)]],
                   n_boot)
    rel_a <- sb(ca[[paste0("odd.", est_col)]], ca[[paste0("even.", est_col)]])
    rel_b <- sb(cb[[paste0("odd.", est_col)]], cb[[paste0("even.", est_col)]])
    data.frame(
      data_type = "SIMULATION", scenario = scenario,
      n_trials_a = n_trials_a, n_trials_b = n_trials_b,
      component = component, measure = est_col,
      realized_planted_r = truth,
      estimated_r = bc["r"], ci_lo = bc["ci_lo"], ci_hi = bc["ci_hi"],
      truth_in_ci = bc["ci_lo"] <= truth & truth <= bc["ci_hi"],
      reliability_a = rel_a, reliability_b = rel_b,
      disattenuated_r = unname(bc["r"]) / sqrt(rel_a * rel_b)
    )
  }
  rows <- rbind(
    measure_row("level", "mean_conf", "beta"),
    measure_row("sensitivity", "auroc2", "lambda"),
    measure_row("sensitivity", "gamma", "lambda")
  )

  gap_boot <- replicate(n_boot, {
    idx <- sample.int(n_subjects, replace = TRUE)
    raw <- cor(ca$all.mean_conf[idx], cb$all.mean_conf[idx]) -
      cor(ca$all.auroc2[idx], cb$all.auroc2[idx])
    raw
  })
  raw_gap <- cor(ca$all.mean_conf, cb$all.mean_conf) -
    cor(ca$all.auroc2, cb$all.auroc2)
  gap <- data.frame(
    data_type = "SIMULATION", scenario = scenario,
    n_trials_a = n_trials_a, n_trials_b = n_trials_b,
    realized_planted_gap = rows$realized_planted_r[1] - rows$realized_planted_r[2],
    estimated_raw_gap = raw_gap,
    gap_ci_lo = quantile(gap_boot, 0.025, names = FALSE),
    gap_ci_hi = quantile(gap_boot, 0.975, names = FALSE),
    disattenuated_gap = rows$disattenuated_r[1] - rows$disattenuated_r[2]
  )
  list(rows = rows, gap = gap)
}

configs <- list(c(200, 200), c(100, 100), c(200, 100))
scenarios <- list(
  gap = c(rho_level = 0.6, rho_sens = 0.2),
  null = c(rho_level = 0.4, rho_sens = 0.4)
)

run_detect_suite <- function(seed) {
  out <- list()
  for (sc in names(scenarios)) {
    for (cfg in configs) {
      out[[length(out) + 1]] <- detectability(
        seed, sc, scenarios[[sc]]["rho_level"], scenarios[[sc]]["rho_sens"],
        cfg[1], cfg[2]
      )
    }
  }
  list(
    rows = do.call(rbind, lapply(out, `[[`, "rows")),
    gaps = do.call(rbind, lapply(out, `[[`, "gap"))
  )
}

message("\n=== THESIS-DETECTABILITY (SIMULATION, main seed) ===")
det_main <- run_detect_suite(seed_main)
print(det_main$rows, digits = 3, row.names = FALSE)
print(det_main$gaps, digits = 3, row.names = FALSE)

message("\n=== SECOND SEED SPOT-CHECK (SIMULATION) ===")
rec_check <- rbind(
  recovery_leakage(seed_check, 200),
  recovery_leakage(seed_check, 100)
)
det_check <- run_detect_suite(seed_check)

seed_compare <- merge(
  rec_main[, c("n_trials", "kind", "estimated", "planted", "r")],
  rec_check[, c("n_trials", "kind", "estimated", "planted", "r")],
  by = c("n_trials", "kind", "estimated", "planted"),
  suffixes = c("_seed1", "_seed2")
)
seed_compare$data_type <- "SIMULATION"
print(seed_compare, digits = 3, row.names = FALSE)

gap_compare <- merge(
  det_main$gaps[, c("scenario", "n_trials_a", "n_trials_b",
                    "estimated_raw_gap", "disattenuated_gap")],
  det_check$gaps[, c("scenario", "n_trials_a", "n_trials_b",
                     "estimated_raw_gap", "disattenuated_gap")],
  by = c("scenario", "n_trials_a", "n_trials_b"),
  suffixes = c("_seed1", "_seed2")
)
gap_compare$data_type <- "SIMULATION"
print(gap_compare, digits = 3, row.names = FALSE)

crit <- function(df, n_trials, est, planted) {
  df[df$n_trials == n_trials & df$estimated == est & df$planted == planted, ]
}
lev200 <- crit(rec_main, 200, "all.mean_conf", "beta")
lev100 <- crit(rec_main, 100, "all.mean_conf", "beta")
sen200 <- crit(rec_main, 200, "all.auroc2", "lambda")
gam200 <- crit(rec_main, 200, "all.gamma", "lambda")
leak_l200 <- crit(rec_main, 200, "all.mean_conf", "lambda")
leak_s200 <- crit(rec_main, 200, "all.auroc2", "beta")

leak_pass <- function(x) abs(x$r) < 0.10 && x$ci_lo <= 0 && x$ci_hi >= 0

gap_rows <- det_main$gaps[det_main$gaps$scenario == "gap", ]
null_rows <- det_main$gaps[det_main$gaps$scenario == "null", ]
det_rows_gap <- det_main$rows[det_main$rows$scenario == "gap" &
                                det_main$rows$measure != "gamma", ]

criteria <- list(
  level_recovery_200 = list(r = lev200$r, threshold = 0.9,
                            pass = lev200$r >= 0.9),
  level_recovery_100 = list(r = lev100$r, threshold = 0.9,
                            pass = lev100$r >= 0.9),
  sensitivity_auroc2_recovery_200 = list(r = sen200$r, threshold = 0.7,
                                         pass = sen200$r >= 0.7),
  sensitivity_gamma_recovery_200 = list(r = gam200$r, threshold = 0.7,
                                        pass = gam200$r >= 0.7),
  leakage_level_from_sensitivity_200 = list(
    r = leak_l200$r, ci = c(leak_l200$ci_lo, leak_l200$ci_hi),
    pass = leak_pass(leak_l200)
  ),
  leakage_sensitivity_from_level_200 = list(
    r = leak_s200$r, ci = c(leak_s200$ci_lo, leak_s200$ci_hi),
    pass = leak_pass(leak_s200)
  ),
  detectability_gap_recovered = list(
    level_gt_sensitivity = all(
      tapply(seq_len(nrow(det_rows_gap)), paste(det_rows_gap$n_trials_a,
                                                det_rows_gap$n_trials_b),
             function(i) {
               d <- det_rows_gap[i, ]
               d$estimated_r[d$component == "level"] >
                 d$estimated_r[d$component == "sensitivity"]
             })
    ),
    truth_in_ci_all = all(det_rows_gap$truth_in_ci),
    pass = NA
  ),
  null_no_spurious_raw_gap = list(
    raw_gap_ci_excludes_zero = null_rows$gap_ci_lo > 0,
    disattenuated_gaps = null_rows$disattenuated_gap,
    pass = NA
  )
)
criteria$detectability_gap_recovered$pass <-
  criteria$detectability_gap_recovered$level_gt_sensitivity &&
  criteria$detectability_gap_recovered$truth_in_ci_all
criteria$null_no_spurious_raw_gap$pass <-
  !any(criteria$null_no_spurious_raw_gap$raw_gap_ci_excludes_zero) ||
  all(abs(null_rows$disattenuated_gap) < 0.10)

gate_passed <- criteria$level_recovery_200$pass &&
  criteria$level_recovery_100$pass &&
  criteria$sensitivity_auroc2_recovery_200$pass &&
  criteria$leakage_level_from_sensitivity_200$pass &&
  criteria$leakage_sensitivity_from_level_200$pass &&
  criteria$detectability_gap_recovered$pass &&
  criteria$null_no_spurious_raw_gap$pass

message("\n=== PRE-SPECIFIED CRITERIA (SIMULATION) ===")
str(criteria)
message("GATE (step 0b portion): ", if (gate_passed) "PASSED" else "FAILED")

write.csv(rbind(rec_main, rec_asymp),
          here("results", "step0b_SIMULATION_recovery_leakage.csv"),
          row.names = FALSE)
write.csv(det_main$rows, here("results", "step0b_SIMULATION_detectability.csv"),
          row.names = FALSE)
write.csv(det_main$gaps, here("results", "step0b_SIMULATION_detectability_gaps.csv"),
          row.names = FALSE)
write.csv(seed_compare, here("results", "step0b_SIMULATION_seedcheck_recovery.csv"),
          row.names = FALSE)
write.csv(gap_compare, here("results", "step0b_SIMULATION_seedcheck_gaps.csv"),
          row.names = FALSE)

meta <- list(
  data_type = "SIMULATION",
  script = "analysis/01_recovery_separability.R",
  run_utc = format(Sys.time(), "%Y-%m-%d %H:%M:%S", tz = "UTC"),
  r_version = R.version.string,
  seeds = list(main = seed_main, second = seed_check),
  n_subjects = n_subjects,
  n_boot = n_boot,
  d_prime_mean = d_prime_mean,
  trait_distributions = list(
    d_prime = "truncnorm(mean 1.2, sd 0.3, min 0.2)",
    lambda = "uniform(0.15, 0.85)",
    beta = "normal(0, 0.8)"
  ),
  empirical_calibration_touch = paste(
    "read results/step0a_basic_summary.csv for per-domain accuracy range",
    sprintf("(%.3f-%.3f) to set d' target; no thesis analysis on real data",
            acc_range[1], acc_range[2])
  ),
  selfcheck = selfcheck_flags,
  asymptotic_recovery_5000_trials = list(
    level = crit(rec_asymp, 5000, "all.mean_conf", "beta")$r,
    auroc2 = crit(rec_asymp, 5000, "all.auroc2", "lambda")$r,
    gamma = crit(rec_asymp, 5000, "all.gamma", "lambda")$r
  ),
  criteria = criteria,
  gate_step0b_passed = gate_passed
)
jsonlite::write_json(meta, here("results", "step0b_SIMULATION_meta.json"),
                     auto_unbox = TRUE, pretty = TRUE, digits = 8)

message("\nStep 0b complete. Outputs in results/step0b_SIMULATION_*.")
