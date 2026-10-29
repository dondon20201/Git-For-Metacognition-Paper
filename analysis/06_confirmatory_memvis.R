library(here)
library(lavaan)
source(here("R", "measures.R"))
source(here("R", "parcels.R"))
source(here("R", "gate_harness.R"))

set.seed(20260712)
n_boot <- 2000

stopifnot(identical(
  unname(tools::md5sum(here("data", "raw", "lund2025",
                            "metacognition_trialData_master.csv"))),
  "500b2e6e1eb0e47e717feba5b912ad9f"
))
input_md5 <- c(
  step1_REALDATA_subject_measures.csv = "97bdda9a331fac3baa59c20e96c4b06b",
  step0b2_SIMULATION_calibration_replicates.csv = "151143c2a1b135c4d2ccd147de6c5070",
  step0b2_SIMULATION_calibration_coefs.csv = "6b03f11342321c2ba24419c038a89ad0"
)
lf_md5 <- vapply(names(input_md5), function(f) {
  path <- here("results", f)
  if (!file.exists(path)) return(NA_character_)
  bytes <- readBin(path, "raw", file.size(path))
  crlf <- bytes == as.raw(13) & c(bytes[-1], as.raw(0)) == as.raw(10)
  lf_file <- tempfile()
  on.exit(unlink(lf_file))
  writeBin(bytes[!crlf], lf_file)
  unname(tools::md5sum(lf_file))
}, character(1))
if (!identical(lf_md5, input_md5)) {
  warning("Input files differ from the versions used in the paper; results may differ slightly from the published numbers.")
}

step1 <- read.csv(here("results", "step1_REALDATA_subject_measures.csv"))

cal_reps <- read.csv(here("results",
                          "step0b2_SIMULATION_calibration_replicates.csv"))
cal_coefs <- read.csv(here("results",
                           "step0b2_SIMULATION_calibration_coefs.csv"))
d200 <- cal_reps[cal_reps$config == "200/200" & cal_reps$converged &
                   cal_reps$admissible, ]
cal_lm <- lm(rho_sens ~ realized_sens, data = d200)
frozen <- list(a = unname(coef(cal_lm)[1]), b = unname(coef(cal_lm)[2]),
               V = vcov(cal_lm))
pub <- cal_coefs[cal_coefs$config == "200/200", ]
stopifnot(abs(frozen$a - pub$intercept_a) < 1e-10,
          abs(frozen$b - pub$slope_b) < 1e-10)
message(sprintf("Frozen calibration re-derived: a = %.4f, b = %.4f (matches committed coefs)",
                frozen$a, frozen$b))

correct_sens <- function(rs, rs_se, cl) {
  est <- (rs - cl$a) / cl$b
  grad <- c(-1 / cl$b, -est / cl$b)
  var_cal <- as.numeric(t(grad) %*% cl$V %*% grad)
  se <- sqrt(rs_se^2 / cl$b^2 + var_cal)
  list(est = est, se = se)
}

build_dat <- function(trial_set) {
  m <- step1[step1$trial_set == trial_set & step1$domain == "memory" &
               step1$qualifies, ]
  v <- step1[step1$trial_set == trial_set & step1$domain == "vision" &
               step1$qualifies, ]
  both <- intersect(m$subject, v$subject)
  m <- m[match(both, m$subject), ]
  v <- v[match(both, v$subject), ]
  dat <- data.frame(subject = both)
  for (j in 1:4) {
    dat[[paste0("sA", j)]] <- m[[paste0("parcel_auroc2_", j)]]
    dat[[paste0("lA", j)]] <- m[[paste0("parcel_mean_conf_", j)]]
    dat[[paste0("sB", j)]] <- v[[paste0("parcel_auroc2_", j)]]
    dat[[paste0("lB", j)]] <- v[[paste0("parcel_mean_conf_", j)]]
  }
  dat$dhatA <- m$d_prime_hat
  dat$dhatB <- v$d_prime_hat
  dat$mean_conf_A <- m$mean_conf
  dat$mean_conf_B <- v$mean_conf
  dat$auroc2_A <- m$auroc2
  dat$auroc2_B <- v$auroc2
  dat
}

analyze <- function(dat, label) {
  n_both <- nrow(dat)
  ind <- c(paste0("sA", 1:4), paste0("sB", 1:4),
           paste0("lA", 1:4), paste0("lB", 1:4), "dhatA", "dhatB")
  n_complete <- sum(complete.cases(dat[ind]))
  datp <- partial_out_dprime(partial_out_dprime(dat, 4L, "A"), 4L, "B")
  fit <- fit_latent(datp, 4L, 4L)
  cs <- correct_sens(fit$rho_sens, fit$rho_sens_se, frozen)
  gap_c <- fit$rho_level - cs$est
  gap_c_var <- fit$rho_level_se^2 + cs$se^2 - 2 * fit$cov_rs_rl / frozen$b
  gap_c_se <- sqrt(max(gap_c_var, 0))
  list(
    label = label, n_both = n_both, n_complete = n_complete, fit = fit,
    rho_sens_c = cs$est, rho_sens_c_se = cs$se,
    rho_sens_c_lo = cs$est - 1.96 * cs$se, rho_sens_c_hi = cs$est + 1.96 * cs$se,
    gap_c = gap_c, gap_c_se = gap_c_se,
    gap_c_lo = gap_c - 1.96 * gap_c_se, gap_c_hi = gap_c + 1.96 * gap_c_se
  )
}

dat_primary <- build_dat("primary")
prim <- analyze(dat_primary, "primary")
message(sprintf("PRIMARY sample: n both-qualifying = %d; n complete (listwise) = %d",
                prim$n_both, prim$n_complete))
if (!prim$fit$converged || !prim$fit$admissible) {
  stop("PRIMARY CFA non-converged/inadmissible after fallback — STOPPING per pre-specified rule 5.")
}
stopifnot(
  (prim$rho_sens_c_hi - prim$rho_sens_c_lo) >=
    (prim$fit$rho_sens_hi - prim$fit$rho_sens_lo),
  (prim$gap_c_hi - prim$gap_c_lo) >= (prim$fit$gap_hi - prim$fit$gap_lo)
)
message(sprintf(
  "PRIMARY (CONFIRMATORY): rho_level = %.3f [%.3f, %.3f]; corrected rho_sens = %.3f [%.3f, %.3f]; corrected gap = %.3f [%.3f, %.3f]",
  prim$fit$rho_level, prim$fit$rho_level_lo, prim$fit$rho_level_hi,
  prim$rho_sens_c, prim$rho_sens_c_lo, prim$rho_sens_c_hi,
  prim$gap_c, prim$gap_c_lo, prim$gap_c_hi
))

boot_gap <- rep(NA_real_, n_boot)
boot_rs_c <- rep(NA_real_, n_boot)
boot_bad <- 0L
for (i in seq_len(n_boot)) {
  bd <- dat_primary[sample.int(nrow(dat_primary), replace = TRUE), ]
  bdp <- partial_out_dprime(partial_out_dprime(bd, 4L, "A"), 4L, "B")
  bf <- fit_latent(bdp, 4L, 4L)
  if (!bf$converged || !bf$admissible) {
    boot_bad <- boot_bad + 1L
    next
  }
  bcs <- (bf$rho_sens - frozen$a) / frozen$b
  boot_rs_c[i] <- bcs
  boot_gap[i] <- bf$rho_level - bcs
}
boot_gap_ci <- quantile(boot_gap, c(0.025, 0.975), na.rm = TRUE, names = FALSE)
boot_rs_ci <- quantile(boot_rs_c, c(0.025, 0.975), na.rm = TRUE, names = FALSE)
message(sprintf(
  "SECONDARY bootstrap (2000): corrected gap CI [%.3f, %.3f]; %d bad fits (%.1f%%)",
  boot_gap_ci[1], boot_gap_ci[2], boot_bad, 100 * boot_bad / n_boot
))

sec_all <- analyze(build_dat("all"), "all_trials")

raw_level <- boot_cor(dat_primary$mean_conf_A, dat_primary$mean_conf_B, n_boot)
raw_sens <- boot_cor(dat_primary$auroc2_A, dat_primary$auroc2_B, n_boot)
raw_gap <- raw_level["r"] - raw_sens["r"]
raw_gap_boot <- replicate(n_boot, {
  idx <- sample.int(nrow(dat_primary), replace = TRUE)
  cor(dat_primary$mean_conf_A[idx], dat_primary$mean_conf_B[idx],
      use = "complete.obs") -
    cor(dat_primary$auroc2_A[idx], dat_primary$auroc2_B[idx],
        use = "complete.obs")
})
raw_gap_ci <- quantile(raw_gap_boot, c(0.025, 0.975), names = FALSE)

row_of <- function(block, quantity, est, se, lo, hi, n, note = "") {
  data.frame(data_type = "REALDATA", block = block, quantity = quantity,
             estimate = est, se = se, ci_lo = lo, ci_hi = hi, n = n,
             note = note)
}
out <- rbind(
  row_of("PRIMARY_CONFIRMATORY", "rho_level_latent",
         prim$fit$rho_level, prim$fit$rho_level_se,
         prim$fit$rho_level_lo, prim$fit$rho_level_hi, prim$n_complete,
         "raw latent (level needs no correction; bias <= 0.006 in gates)"),
  row_of("PRIMARY_CONFIRMATORY", "rho_sens_latent_corrected",
         prim$rho_sens_c, prim$rho_sens_c_se,
         prim$rho_sens_c_lo, prim$rho_sens_c_hi, prim$n_complete,
         "d'-partialed, frozen 200/200 correction, delta-method CI"),
  row_of("PRIMARY_CONFIRMATORY", "gap_corrected",
         prim$gap_c, prim$gap_c_se, prim$gap_c_lo, prim$gap_c_hi,
         prim$n_complete, "THE pre-specified primary outcome"),
  row_of("SECONDARY_a_bootstrap", "gap_corrected_bootCI",
         prim$gap_c, NA, boot_gap_ci[1], boot_gap_ci[2], prim$n_complete,
         sprintf("subject bootstrap, calibration fixed; %d/%d bad fits",
                 boot_bad, n_boot)),
  row_of("SECONDARY_a_bootstrap", "rho_sens_corrected_bootCI",
         prim$rho_sens_c, NA, boot_rs_ci[1], boot_rs_ci[2], prim$n_complete, ""),
  row_of("SECONDARY_b_alltrials", "rho_level_latent",
         sec_all$fit$rho_level, sec_all$fit$rho_level_se,
         sec_all$fit$rho_level_lo, sec_all$fit$rho_level_hi,
         sec_all$n_complete, ""),
  row_of("SECONDARY_b_alltrials", "rho_sens_latent_corrected",
         sec_all$rho_sens_c, sec_all$rho_sens_c_se,
         sec_all$rho_sens_c_lo, sec_all$rho_sens_c_hi, sec_all$n_complete, ""),
  row_of("SECONDARY_b_alltrials", "gap_corrected",
         sec_all$gap_c, sec_all$gap_c_se, sec_all$gap_c_lo, sec_all$gap_c_hi,
         sec_all$n_complete, ""),
  row_of("SECONDARY_c_uncorrected", "rho_sens_latent_uncorrected",
         prim$fit$rho_sens, prim$fit$rho_sens_se,
         prim$fit$rho_sens_lo, prim$fit$rho_sens_hi, prim$n_complete,
         "uncorrected-latent exhibit"),
  row_of("SECONDARY_c_uncorrected", "gap_uncorrected",
         prim$fit$gap, prim$fit$gap_se, prim$fit$gap_lo, prim$fit$gap_hi,
         prim$n_complete, ""),
  row_of("SECONDARY_d_manifest_EXHIBIT", "raw_r_mean_conf",
         raw_level["r"], NA, raw_level["ci_lo"], raw_level["ci_hi"],
         prim$n_both, "METHODS EXHIBIT ONLY — raw contrasts invent the thesis under a null"),
  row_of("SECONDARY_d_manifest_EXHIBIT", "raw_r_auroc2",
         raw_sens["r"], NA, raw_sens["ci_lo"], raw_sens["ci_hi"], prim$n_both,
         "unpartialed subject-level AUROC2"),
  row_of("SECONDARY_d_manifest_EXHIBIT", "raw_gap_naive",
         raw_gap, NA, raw_gap_ci[1], raw_gap_ci[2], prim$n_both,
         "METHODS EXHIBIT ONLY")
)
out$converged <- prim$fit$converged
out$admissible <- prim$fit$admissible
write.csv(out, here("results", "step6_REALDATA_confirmatory.csv"),
          row.names = FALSE)

gap <- prim$gap_c
lo <- prim$gap_c_lo
hi <- prim$gap_c_hi
fmt <- function(x) sprintf("%.2f", x)
primary_sentence <- if (lo > 0) {
  sprintf("On the confirmatory memory×vision pair, confidence level was more domain-general than confidence informativeness (corrected gap = %s, 95%% CI [%s, %s]).",
          fmt(gap), fmt(lo), fmt(hi))
} else if (hi < 0) {
  sprintf("The confirmatory test showed the REVERSE pattern (corrected gap = %s, 95%% CI [%s, %s]).",
          fmt(gap), fmt(lo), fmt(hi))
} else {
  sprintf("The confirmatory test does not exclude a null gap (corrected gap = %s, 95%% CI [%s, %s]). Given ~60%% simulated power at a true gap of 0.4 (plausibly somewhat higher given memory's measured reliability), a null gap is **not excluded**; this is not evidence of absence.",
          fmt(gap), fmt(lo), fmt(hi))
}

md <- c(
  "# CONFIRMATORY RESULT — memory × vision",
  "",
  sprintf("**Date:** %s · **Code:** `analysis/06_confirmatory_memvis.R` · **Data:** REALDATA, primary trial set · n = %d both-qualifying, %d complete-case (listwise, pre-specified)",
          format(Sys.Date()), prim$n_both, prim$n_complete),
  "",
  "## PRIMARY (CONFIRMATORY — pre-specified wording verbatim)",
  "",
  primary_sentence,
  "",
  "| Quantity | Estimate | 95% CI (delta) | Status |",
  "|---|---|---|---|",
  sprintf("| Latent ρ_level (raw) | %.3f | [%.3f, %.3f] | VERIFIED, CONFIRMATORY |",
          prim$fit$rho_level, prim$fit$rho_level_lo, prim$fit$rho_level_hi),
  sprintf("| Latent ρ_sens (d′-partialed, corrected) | %.3f | [%.3f, %.3f] | VERIFIED, CONFIRMATORY |",
          prim$rho_sens_c, prim$rho_sens_c_lo, prim$rho_sens_c_hi),
  sprintf("| **Corrected gap (primary outcome)** | **%.3f** | **[%.3f, %.3f]** | VERIFIED, CONFIRMATORY |",
          prim$gap_c, prim$gap_c_lo, prim$gap_c_hi),
  "",
  sprintf("CFA converged and admissible: %s/%s. Corrected CI widths ≥ uncorrected (asserted in code). Frozen calibration a = %.4f, b = %.4f.",
          prim$fit$converged, prim$fit$admissible, frozen$a, frozen$b),
  "",
  "## SECONDARY",
  "",
  sprintf("- (a) Subject bootstrap (2000, calibration fixed): corrected gap CI [%.3f, %.3f]; corrected ρ_sens CI [%.3f, %.3f]; bad fits %d/2000.",
          boot_gap_ci[1], boot_gap_ci[2], boot_rs_ci[1], boot_rs_ci[2], boot_bad),
  sprintf("- (b) All-trials variant: gap = %.3f [%.3f, %.3f] (n = %d).",
          sec_all$gap_c, sec_all$gap_c_lo, sec_all$gap_c_hi, sec_all$n_complete),
  sprintf("- (c) Uncorrected latent: ρ_sens = %.3f [%.3f, %.3f]; gap = %.3f [%.3f, %.3f].",
          prim$fit$rho_sens, prim$fit$rho_sens_lo, prim$fit$rho_sens_hi,
          prim$fit$gap, prim$fit$gap_lo, prim$fit$gap_hi),
  sprintf("- (d) METHODS EXHIBIT ONLY (never evidence): raw manifest r(mean conf) = %.3f [%.3f, %.3f]; raw r(AUROC2) = %.3f [%.3f, %.3f]; naive raw gap = %.3f [%.3f, %.3f].",
          raw_level["r"], raw_level["ci_lo"], raw_level["ci_hi"],
          raw_sens["r"], raw_sens["ci_lo"], raw_sens["ci_hi"],
          raw_gap, raw_gap_ci[1], raw_gap_ci[2]),
  ""
)
writeLines(md, here("results", "step6_REALDATA_confirmatory.md"))
message("\nConfirmatory outputs written. THE PRIMARY NUMBERS ARE NOW FINAL.")
