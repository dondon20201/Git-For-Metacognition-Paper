library(here)
library(lavaan)
source(here("R", "measures.R"))
source(here("R", "simulate.R"))
source(here("R", "parcels.R"))

seed_main <- 20260710
seed_check <- 42424242
n_subjects <- 322
n_rep_main <- 30
n_rep_check <- 15

scenarios <- list(
  GAP = c(level = 0.6, sens = 0.2, dprime = 0),
  NULLGAP = c(level = 0.4, sens = 0.4, dprime = 0),
  ZEROSENS = c(level = 0.4, sens = 0.0, dprime = 0),
  DCONF = c(level = 0.0, sens = 0.0, dprime = 0.4)
)
configs <- list(c(200, 200), c(200, 100))

correlated_pair <- function(n, rho) {
  z <- matrix(rnorm(n * 2), n) %*% chol(matrix(c(1, rho, rho, 1), 2))
  list(a = z[, 1], b = z[, 2])
}

draw_correlated_traits <- function(n, rho_level, rho_sens, rho_dprime) {
  zb <- correlated_pair(n, rho_level)
  zl <- correlated_pair(n, rho_sens)
  zd <- correlated_pair(n, rho_dprime)
  mk <- function(z_beta, z_lambda, z_d) data.frame(
    d_prime = pmax(1.2 + 0.3 * z_d, 0.2),
    lambda = 0.15 + 0.7 * pnorm(z_lambda),
    beta = 0.8 * z_beta
  )
  list(a = mk(zb$a, zl$a, zd$a), b = mk(zb$b, zl$b, zd$b))
}

sim_domain <- function(traits, n_trials, k, prefix) {
  res <- lapply(seq_len(nrow(traits)), function(i) {
    tr <- simulate_trials(
      n_trials, traits$d_prime[i], traits$lambda[i], traits$beta[i]
    )
    c(parcel_measures(tr, k), dhat = d_prime_hat(tr))
  })
  m <- do.call(rbind, res)
  colnames(m) <- c(
    paste0("s", prefix, seq_len(k)), paste0("l", prefix, seq_len(k)),
    paste0("dhat", prefix)
  )
  as.data.frame(m)
}

partial_out_dprime <- function(dat, k, prefix) {
  d <- dat[[paste0("dhat", prefix)]]
  for (j in seq_len(k)) {
    col <- paste0("s", prefix, j)
    dat[[col]] <- resid(lm(dat[[col]] ~ d, na.action = na.exclude))
  }
  dat
}

run_replicate <- function(scenario, rhos, n_trials_a, n_trials_b, seed, rep) {
  set.seed(seed)
  tr <- draw_correlated_traits(n_subjects, rhos["level"], rhos["sens"],
                               rhos["dprime"])
  k_a <- if (n_trials_a >= 200) 4L else 2L
  k_b <- if (n_trials_b >= 200) 4L else 2L
  dat <- cbind(sim_domain(tr$a, n_trials_a, k_a, "A"),
               sim_domain(tr$b, n_trials_b, k_b, "B"))
  realized <- c(
    level = cor(tr$a$beta, tr$b$beta),
    sens = cor(tr$a$lambda, tr$b$lambda),
    dprime = cor(tr$a$d_prime, tr$b$d_prime)
  )
  base <- data.frame(
    data_type = "SIMULATION", scenario = scenario, rep = rep, seed = seed,
    n_trials_a = n_trials_a, n_trials_b = n_trials_b, k_a = k_a, k_b = k_b,
    realized_level = realized["level"], realized_sens = realized["sens"],
    realized_dprime = realized["dprime"],
    row.names = NULL
  )
  raw <- cbind(base, variant = "raw", fit_latent(dat, k_a, k_b))
  part <- cbind(base, variant = "dprime_partialed",
                fit_latent(partial_out_dprime(partial_out_dprime(dat, k_a, "A"),
                                              k_b, "B"), k_a, k_b))
  rbind(raw, part)
}

run_suite <- function(base_seed, n_rep, label) {
  rows <- list()
  idx <- 0
  for (sc in names(scenarios)) {
    for (cfg in configs) {
      message(sprintf("[%s] scenario %s, trials %d/%d ...",
                      label, sc, cfg[1], cfg[2]))
      for (r in seq_len(n_rep)) {
        idx <- idx + 1
        rows[[idx]] <- run_replicate(
          sc, scenarios[[sc]], cfg[1], cfg[2], base_seed + idx * 1009, r
        )
      }
    }
  }
  out <- do.call(rbind, rows)
  out$suite <- label
  out
}

t0 <- Sys.time()
res_main <- run_suite(seed_main, n_rep_main, "main")
res_check <- run_suite(seed_check, n_rep_check, "seed2")
message(sprintf("Total fit time: %.1f min",
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))

summarize_cell <- function(d) {
  ok <- d$converged & d$admissible
  realized_gap <- d$realized_level - d$realized_sens
  data.frame(
    data_type = "SIMULATION",
    scenario = d$scenario[1], variant = d$variant[1],
    n_trials_a = d$n_trials_a[1], n_trials_b = d$n_trials_b[1],
    n_rep = nrow(d), n_ok = sum(ok),
    pct_bad = 100 * mean(!ok),
    median_rho_sens = median(d$rho_sens[ok]),
    median_sens_bias = median((d$rho_sens - d$realized_sens)[ok]),
    sens_coverage = sum(ok & d$rho_sens_lo <= d$realized_sens &
                          d$realized_sens <= d$rho_sens_hi),
    median_rho_level = median(d$rho_level[ok]),
    median_level_bias = median((d$rho_level - d$realized_level)[ok]),
    median_gap = median(d$gap[ok]),
    median_gap_bias = median((d$gap - realized_gap)[ok]),
    gap_positive = sum(ok & d$gap > 0),
    gap_ci_covers_zero = sum(ok & d$gap_lo <= 0 & 0 <= d$gap_hi)
  )
}

summarize_suite <- function(res) {
  cells <- split(res, list(res$scenario, res$variant, res$n_trials_a,
                           res$n_trials_b), drop = TRUE)
  out <- do.call(rbind, lapply(cells, summarize_cell))
  rownames(out) <- NULL
  out[order(out$scenario, out$variant, -out$n_trials_b), ]
}

sum_main <- summarize_suite(res_main)
sum_check <- summarize_suite(res_check)

criteria_for <- function(smry, variant, n_rep) {
  g <- function(sc, ta, tb) {
    smry[smry$scenario == sc & smry$variant == variant &
           smry$n_trials_a == ta & smry$n_trials_b == tb, ]
  }
  per_config <- lapply(configs, function(cfg) {
    gap <- g("GAP", cfg[1], cfg[2])
    nul <- g("NULLGAP", cfg[1], cfg[2])
    zer <- g("ZEROSENS", cfg[1], cfg[2])
    dcf <- g("DCONF", cfg[1], cfg[2])
    all_cells <- rbind(gap, nul, zer, dcf)
    list(
      config = paste0(cfg[1], "/", cfg[2]),
      a_gap_bias = abs(gap$median_sens_bias) <= 0.07,
      a_gap_coverage = gap$sens_coverage >= ceiling(25 / 30 * n_rep),
      a_gap_positive = gap$gap_positive >= ceiling(27 / 30 * n_rep),
      b_null_gap_bias = abs(nul$median_gap_bias) <= 0.07,
      b_null_gap_ci_zero = nul$gap_ci_covers_zero >= ceiling(25 / 30 * n_rep),
      c_zerosens_median = abs(zer$median_rho_sens) <= 0.07,
      d_dconf_spurious = dcf$median_rho_sens < 0.10,
      e_convergence = all(all_cells$pct_bad <= 20),
      all_pass = NA
    )
  })
  per_config <- lapply(per_config, function(x) {
    x$all_pass <- all(unlist(x[setdiff(names(x), c("config", "all_pass"))]))
    x
  })
  names(per_config) <- vapply(per_config, `[[`, "", "config")
  per_config
}

crit_raw <- criteria_for(sum_main, "raw", n_rep_main)
crit_part <- criteria_for(sum_main, "dprime_partialed", n_rep_main)

dconf_raw_fails <- !all(vapply(crit_raw, `[[`, TRUE, "d_dconf_spurious"))
dconf_part_passes <- all(vapply(crit_part, `[[`, TRUE, "d_dconf_spurious"))
primary_variant <- if (dconf_raw_fails && dconf_part_passes) {
  "dprime_partialed"
} else {
  "raw"
}
crit_primary <- if (primary_variant == "raw") crit_raw else crit_part
gate_passed <- all(vapply(crit_primary, `[[`, TRUE, "all_pass"))

message("\n=== SUMMARY (SIMULATION, main seed) ===")
print(sum_main, digits = 3, row.names = FALSE)
message("\nPrimary variant: ", primary_variant)
message("\n=== PRE-SPECIFIED CRITERIA, primary variant ===")
str(crit_primary)
message("\nSTEP 0b' GATE: ", if (gate_passed) "PASSED" else "FAILED")

message("\n=== SECOND-SEED SPOT-CHECK (medians; ", n_rep_check, " reps) ===")
print(sum_check[sum_check$variant == primary_variant,
                c("scenario", "n_trials_a", "n_trials_b", "median_rho_sens",
                  "median_sens_bias", "median_gap", "median_gap_bias",
                  "pct_bad")],
      digits = 3, row.names = FALSE)

set.seed(seed_main)
trait_check <- do.call(rbind, lapply(c(0.2, 0.4, 0.6), function(rho) {
  r_lambda <- replicate(2000, {
    z <- correlated_pair(n_subjects, rho)
    cor(0.15 + 0.7 * pnorm(z$a), 0.15 + 0.7 * pnorm(z$b))
  })
  r_beta <- replicate(2000, {
    z <- correlated_pair(n_subjects, rho)
    cor(0.8 * z$a, 0.8 * z$b)
  })
  data.frame(
    data_type = "SIMULATION", population_rho = rho,
    lambda_mean_realized = mean(r_lambda), lambda_sd_realized = sd(r_lambda),
    lambda_q025 = quantile(r_lambda, 0.025, names = FALSE),
    lambda_q975 = quantile(r_lambda, 0.975, names = FALSE),
    beta_mean_realized = mean(r_beta), beta_sd_realized = sd(r_beta),
    pctile_of_0.322_when_rho_0.2 = if (rho == 0.2) mean(r_lambda <= 0.322)
                                   else NA_real_
  )
}))
message("\n=== TRAIT-TRANSFORM CHECK (SIMULATION) ===")
print(trait_check, digits = 3, row.names = FALSE)

write.csv(rbind(res_main, res_check),
          here("results", "step0bprime_SIMULATION_replicates.csv"),
          row.names = FALSE)
write.csv(rbind(cbind(sum_main, suite = "main"),
                cbind(sum_check, suite = "seed2")),
          here("results", "step0bprime_SIMULATION_summary.csv"),
          row.names = FALSE)
write.csv(trait_check,
          here("results", "step0bprime_SIMULATION_traitcheck.csv"),
          row.names = FALSE)

meta <- list(
  data_type = "SIMULATION",
  script = "analysis/02_latent_gate.R",
  run_utc = format(Sys.time(), "%Y-%m-%d %H:%M:%S", tz = "UTC"),
  r_version = R.version.string,
  lavaan_version = as.character(packageVersion("lavaan")),
  seeds = list(main = seed_main, second = seed_check),
  n_subjects = n_subjects,
  n_rep = list(main = n_rep_main, second = n_rep_check),
  parcels = "K=4 at 200 trials, K=2 at 100 (equal 50-trial parcels)",
  scenarios = scenarios,
  configs = configs,
  primary_variant = primary_variant,
  criteria_raw = crit_raw,
  criteria_partialed = crit_part,
  gate_step0bprime_passed = gate_passed
)
jsonlite::write_json(meta, here("results", "step0bprime_SIMULATION_meta.json"),
                     auto_unbox = TRUE, pretty = TRUE, digits = 8)

message("\nStep 0b' complete. Outputs in results/step0bprime_SIMULATION_*.")
