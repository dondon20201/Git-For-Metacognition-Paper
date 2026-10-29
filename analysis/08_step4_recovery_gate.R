library(here)
source(here("R", "multilevel.R"))

t_start <- Sys.time()

domains <- data.frame(
  domain = c("memory", "vision", "Calories", "GDP"),
  n_trials = c(200L, 200L, 100L, 100L),
  p_acc = c(0.77, 0.72, 0.71, 0.70)
)
ds <- read.csv(here("results", "step1_REALDATA_domain_summary.csv"))
ds <- ds[ds$trial_set == "primary", ]
domains$beta <- ds$mean_conf_mean[match(domains$domain, ds$domain)]
stopifnot(!anyNA(domains$beta))

n_sub <- 322L
sd_int_total <- 0.8
sd_slope_total <- 0.4
sd_resid <- 1.2
gamma_slope <- 0.8

cells <- data.frame(
  icc_level = c(0.8, 0.8, 0.5, 0.5, 0.2),
  icc_sens = c(0.2, 0.5, 0.2, 0.5, 0.8)
)
n_rep <- 10L
seed_base <- 20260710L

simulate_step4 <- function(icc_level, icc_sens, seed) {
  set.seed(seed)
  u <- rnorm(n_sub, 0, sd_int_total * sqrt(icc_level))
  w <- rnorm(n_sub, 0, sd_slope_total * sqrt(icc_sens))
  out <- vector("list", nrow(domains))
  for (d in seq_len(nrow(domains))) {
    nd <- domains$n_trials[d]
    v <- rnorm(n_sub, 0, sd_int_total * sqrt(1 - icc_level))
    z <- rnorm(n_sub, 0, sd_slope_total * sqrt(1 - icc_sens))
    s <- rep(seq_len(n_sub), each = nd)
    acc <- rbinom(n_sub * nd, 1, domains$p_acc[d])
    acc_c <- acc - domains$p_acc[d]
    conf_star <- domains$beta[d] + u[s] + v[s] +
      (gamma_slope + w[s] + z[s]) * acc_c +
      rnorm(n_sub * nd, 0, sd_resid)
    out[[d]] <- data.frame(
      subject = s, domain = domains$domain[d], accuracy = acc,
      conf = pmin(pmax(round(conf_star), 1L), 7L)
    )
  }
  dat <- do.call(rbind, out)
  dat$domain <- factor(dat$domain, levels = domains$domain)
  dat$subject <- factor(dat$subject)
  dat$subject_domain <- interaction(dat$subject, dat$domain, drop = TRUE)
  add_acc_dc(dat)
}

rows <- list()
ips_rows <- list()
for (ci in seq_len(nrow(cells))) {
  for (r in seq_len(n_rep)) {
    seed <- seed_base + r
    dat <- simulate_step4(cells$icc_level[ci], cells$icc_sens[ci], seed)
    t0 <- Sys.time()
    res <- step4_ladder(dat)
    secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    rows[[length(rows) + 1]] <- data.frame(
      data_type = "SIMULATION",
      cell = sprintf("L%.1f_S%.1f", cells$icc_level[ci], cells$icc_sens[ci]),
      planted_icc_level = cells$icc_level[ci],
      planted_icc_sens = cells$icc_sens[ci],
      replicate = r, seed = seed,
      icc_level_hat = res$icc_level, icc_sens_hat = res$icc_sens,
      delta_icc_hat = res$delta_icc,
      v_int_subject = res$v_int_subject,
      v_int_subject_domain = res$v_int_subject_domain,
      v_slope_subject = res$v_slope_subject,
      v_slope_subject_domain = res$v_slope_subject_domain,
      v_residual = res$v_residual,
      rung = res$rung,
      drop_slope_at = ifelse(is.null(res$drop_slope_at), "",
                             res$drop_slope_at),
      converged = res$converged, singular = res$singular,
      fit_seconds = secs
    )
    message(sprintf(
      "cell (%.1f, %.1f) rep %02d [seed %d]: ICC_level=%.3f ICC_sens=%.3f rung=%s conv=%s sing=%s (%.0fs)",
      cells$icc_level[ci], cells$icc_sens[ci], r, seed,
      res$icc_level, ifelse(is.na(res$icc_sens), NA, res$icc_sens),
      res$rung, res$converged, res$singular, secs
    ))
    if (cells$icc_level[ci] == 0.5 && cells$icc_sens[ci] == 0.5) {
      ips <- ipsative_exhibit(dat)
      ips$data_type <- "SIMULATION"
      ips$replicate <- r
      ips$seed <- seed
      ips$planted_icc_level <- 0.5
      ips$planted_icc_sens <- 0.5
      ips$lmm_icc_level_hat <- res$icc_level
      ips$lmm_icc_sens_hat <- res$icc_sens
      ips$theory_note <- "mechanical expectation ~ -1/3 for 4 exchangeable domains after grand-mean centering"
      ips_rows[[length(ips_rows) + 1]] <- ips
    }
  }
}
rec <- do.call(rbind, rows)
write.csv(rec, here("results", "step8_SIMULATION_step4gate_recovery.csv"),
          row.names = FALSE)
ips_all <- do.call(rbind, ips_rows)
write.csv(ips_all, here("results", "step8_SIMULATION_ipsative_exhibit.csv"),
          row.names = FALSE)

agg <- do.call(rbind, lapply(split(rec, rec$cell), function(d) {
  data.frame(
    cell = d$cell[1],
    planted_icc_level = d$planted_icc_level[1],
    planted_icc_sens = d$planted_icc_sens[1],
    n_rep = nrow(d),
    mean_icc_level = mean(d$icc_level_hat, na.rm = TRUE),
    mean_icc_sens = mean(d$icc_sens_hat, na.rm = TRUE),
    n_icc_sens_estimable = sum(!is.na(d$icc_sens_hat)),
    bias_level = mean(d$icc_level_hat, na.rm = TRUE) - d$planted_icc_level[1],
    bias_sens = mean(d$icc_sens_hat, na.rm = TRUE) - d$planted_icc_sens[1],
    n_sign_correct = sum(!is.na(d$delta_icc_hat) &
                           sign(d$delta_icc_hat) ==
                             sign(d$planted_icc_level[1] -
                                    d$planted_icc_sens[1])),
    n_rung_m1m2 = sum(d$rung %in% c("M1", "M2"))
  )
}))
rownames(agg) <- NULL

cell_of <- function(l, s) agg[agg$planted_icc_level == l &
                                agg$planted_icc_sens == s, ]

agg$pass_recovery <- abs(agg$bias_level) <= 0.10 & abs(agg$bias_sens) <= 0.10
agg$delta_planted <- agg$planted_icc_level - agg$planted_icc_sens
agg$sign_binding <- abs(agg$delta_planted) >= 0.3
agg$pass_sign <- !agg$sign_binding | agg$n_sign_correct >= 9
agg$pass_convergence <- agg$n_rung_m1m2 >= 9

leak <- data.frame(
  comparison = c("sens .2->.5 at level .8 (moves level)",
                 "sens .2->.5 at level .5 (moves level)",
                 "level .8->.5 at sens .2 (moves sens)",
                 "level .8->.5 at sens .5 (moves sens)"),
  delta_other = c(
    cell_of(0.8, 0.5)$mean_icc_level - cell_of(0.8, 0.2)$mean_icc_level,
    cell_of(0.5, 0.5)$mean_icc_level - cell_of(0.5, 0.2)$mean_icc_level,
    cell_of(0.5, 0.2)$mean_icc_sens - cell_of(0.8, 0.2)$mean_icc_sens,
    cell_of(0.5, 0.5)$mean_icc_sens - cell_of(0.8, 0.5)$mean_icc_sens
  )
)
leak$pass <- abs(leak$delta_other) <= 0.10

pass_c1 <- all(agg$pass_recovery)
pass_c2 <- all(leak$pass)
pass_c3 <- all(agg$pass_sign)
pass_c4 <- all(agg$pass_convergence)
overall <- pass_c1 && pass_c2 && pass_c3 && pass_c4

summary_df <- agg
summary_df$data_type <- "SIMULATION"
summary_df$gate_overall <- ifelse(overall, "PASS", "FAIL")
write.csv(summary_df,
          here("results", "step8_SIMULATION_step4gate_summary.csv"),
          row.names = FALSE)

message("\n=== Leakage comparisons ===")
print(leak, row.names = FALSE, digits = 3)
message(sprintf(
  "\n=== STEP 4 GATE: %s ===\nC1 recovery: %s | C2 leakage: %s | C3 sign: %s | C4 convergence: %s",
  ifelse(overall, "PASS", "FAIL"),
  pass_c1, pass_c2, pass_c3, pass_c4
))
message(sprintf(
  "Ipsative exhibit, (0.5,0.5) cell: mean centered pairwise r = %.3f (theory ~ -0.333) while LMM mean ICC_level = %.3f / ICC_sens = %.3f (planted 0.5/0.5)",
  mean(ips_all$r_centered),
  cell_of(0.5, 0.5)$mean_icc_level, cell_of(0.5, 0.5)$mean_icc_sens
))

meta <- list(
  data_type = "SIMULATION",
  design = list(
    n_subjects = n_sub,
    domains = domains,
    sd_int_total = sd_int_total,
    sd_slope_total = sd_slope_total,
    sd_slope_total_note = "IMPLEMENTATION CHOICE — the plan specified mean slope 0.8 and 'slope SDs partitioned by ICC_sens' but not the total slope SD; 0.4 (= mean/2) chosen",
    sd_resid = sd_resid,
    gamma_slope = gamma_slope,
    beta_source = "results/step1_REALDATA_domain_summary.csv (primary, committed)",
    planted_int_slope_correlation = 0,
    discretization = "round(conf*) clamped to 1..7",
    centering_generation = "true p_acc",
    centering_fitting = "empirical domain mean (pre-specified acc_dc)"
  ),
  cells = cells,
  n_replicates_per_cell = n_rep,
  seeds = seed_base + seq_len(n_rep),
  criteria = list(pass_recovery = pass_c1, pass_leakage = pass_c2,
                  pass_sign = pass_c3, pass_convergence = pass_c4,
                  overall = ifelse(overall, "PASS", "FAIL")),
  leakage = leak,
  mean_fit_seconds = mean(rec$fit_seconds),
  total_minutes = as.numeric(difftime(Sys.time(), t_start, units = "mins")),
  r_version = R.version.string,
  lme4_version = as.character(packageVersion("lme4"))
)
jsonlite::write_json(meta, here("results", "step8_SIMULATION_step4gate_meta.json"),
                     auto_unbox = TRUE, digits = 10, pretty = TRUE)
message(sprintf("\nGate outputs written. Total %.1f min.",
                meta$total_minutes))
