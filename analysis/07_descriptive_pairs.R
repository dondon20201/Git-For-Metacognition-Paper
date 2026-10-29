library(here)
library(lavaan)
source(here("R", "measures.R"))
source(here("R", "parcels.R"))
source(here("R", "gate_harness.R"))

stopifnot(identical(
  unname(tools::md5sum(here("data", "raw", "lund2025",
                            "metacognition_trialData_master.csv"))),
  "500b2e6e1eb0e47e717feba5b912ad9f"
))
input_md5 <- c(
  step1_REALDATA_subject_measures.csv = "97bdda9a331fac3baa59c20e96c4b06b",
  step6_REALDATA_confirmatory.csv = "5472f1b8c8e76edb16ec0b0e90dc7336"
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
step6 <- read.csv(here("results", "step6_REALDATA_confirmatory.csv"))

k_domain <- c(memory = 4L, vision = 4L, Calories = 2L, GDP = 2L)
pairs <- data.frame(
  dom_a = c("memory", "memory", "vision", "vision", "Calories"),
  dom_b = c("Calories", "GDP", "Calories", "GDP", "GDP"),
  band  = c(0.09, 0.09, 0.09, 0.09, 0.10)
)

build_pair <- function(dom_a, dom_b) {
  a <- step1[step1$trial_set == "primary" & step1$domain == dom_a &
               step1$qualifies, ]
  b <- step1[step1$trial_set == "primary" & step1$domain == dom_b &
               step1$qualifies, ]
  both <- intersect(a$subject, b$subject)
  a <- a[match(both, a$subject), ]
  b <- b[match(both, b$subject), ]
  dat <- data.frame(subject = both)
  for (j in seq_len(k_domain[[dom_a]])) {
    dat[[paste0("sA", j)]] <- a[[paste0("parcel_auroc2_", j)]]
    dat[[paste0("lA", j)]] <- a[[paste0("parcel_mean_conf_", j)]]
  }
  for (j in seq_len(k_domain[[dom_b]])) {
    dat[[paste0("sB", j)]] <- b[[paste0("parcel_auroc2_", j)]]
    dat[[paste0("lB", j)]] <- b[[paste0("parcel_mean_conf_", j)]]
  }
  dat$dhatA <- a$d_prime_hat
  dat$dhatB <- b$d_prime_hat
  dat
}

rows <- lapply(seq_len(nrow(pairs)), function(i) {
  dom_a <- pairs$dom_a[i]
  dom_b <- pairs$dom_b[i]
  k_a <- k_domain[[dom_a]]
  k_b <- k_domain[[dom_b]]
  dat <- build_pair(dom_a, dom_b)
  ind <- c(paste0("sA", seq_len(k_a)), paste0("sB", seq_len(k_b)),
           paste0("lA", seq_len(k_a)), paste0("lB", seq_len(k_b)),
           "dhatA", "dhatB")
  n_both <- nrow(dat)
  n_complete <- sum(complete.cases(dat[ind]))
  datp <- partial_out_dprime(partial_out_dprime(dat, k_a, "A"), k_b, "B")
  fit <- fit_latent(datp, k_a, k_b)
  estimable <- fit$converged && fit$admissible
  message(sprintf(
    "%s x %s (K=%d/%d): n_both=%d n_complete=%d converged=%s admissible=%s",
    dom_a, dom_b, k_a, k_b, n_both, n_complete, fit$converged, fit$admissible
  ))
  if (estimable) {
    message(sprintf(
      "  rho_level=%.3f [%.3f,%.3f]; rho_sens(uncorr)=%.3f [%.3f,%.3f]; gap=%.3f [%.3f,%.3f] +/- %.2f",
      fit$rho_level, fit$rho_level_lo, fit$rho_level_hi,
      fit$rho_sens, fit$rho_sens_lo, fit$rho_sens_hi,
      fit$gap, fit$gap_lo, fit$gap_hi, pairs$band[i]
    ))
  } else {
    message("  NOT ESTIMABLE — reported as such per pre-specified rule 5.")
  }
  data.frame(
    data_type = "REALDATA", tier = "NON-CONFIRMATORY_DESCRIPTIVE",
    pair = paste0(dom_a, "_x_", dom_b),
    domain_a = dom_a, domain_b = dom_b, k_a = k_a, k_b = k_b,
    n_both = n_both, n_complete = n_complete,
    rho_level = fit$rho_level, rho_level_se = fit$rho_level_se,
    rho_level_lo = fit$rho_level_lo, rho_level_hi = fit$rho_level_hi,
    rho_sens_uncorrected = fit$rho_sens,
    rho_sens_se = fit$rho_sens_se,
    rho_sens_lo = fit$rho_sens_lo, rho_sens_hi = fit$rho_sens_hi,
    gap_uncorrected = fit$gap, gap_se = fit$gap_se,
    gap_lo = fit$gap_lo, gap_hi = fit$gap_hi,
    bias_band = pairs$band[i],
    converged = fit$converged, admissible = fit$admissible,
    estimable = estimable,
    note = if (estimable) {
      "uncorrected latent; gap carries frozen systematic band (0b'' calibration FAILED at this config)"
    } else {
      "NOT ESTIMABLE — CFA non-converged/inadmissible after bounded fallback (pre-specified rule 5)"
    }
  )
})
out <- do.call(rbind, rows)
write.csv(out, here("results", "step7_REALDATA_descriptive.csv"),
          row.names = FALSE)

g6 <- function(q) step6[step6$block == "PRIMARY_CONFIRMATORY" &
                          step6$quantity == q, ]
c_lev <- g6("rho_level_latent")
c_sen <- g6("rho_sens_latent_corrected")
c_gap <- g6("gap_corrected")
c_n_both <- step6[step6$quantity == "raw_r_mean_conf", "n"]

fmt3 <- function(x) sprintf("%.3f", x)
ci <- function(lo, hi) sprintf("[%s, %s]", fmt3(lo), fmt3(hi))

desc_rows <- vapply(seq_len(nrow(out)), function(i) {
  r <- out[i, ]
  if (r$estimable) {
    sprintf(
      "| %s × %s | %d/%d | %d / %d | %s %s | %s %s | %s %s ± %.2f (systematic) | NON-CONFIRMATORY, DESCRIPTIVE |",
      r$domain_a, r$domain_b, r$k_a, r$k_b, r$n_both, r$n_complete,
      fmt3(r$rho_level), ci(r$rho_level_lo, r$rho_level_hi),
      fmt3(r$rho_sens_uncorrected), ci(r$rho_sens_lo, r$rho_sens_hi),
      fmt3(r$gap_uncorrected), ci(r$gap_lo, r$gap_hi), r$bias_band
    )
  } else {
    sprintf(
      "| %s × %s | %d/%d | %d / %d | NOT ESTIMABLE | NOT ESTIMABLE | NOT ESTIMABLE (converged=%s, admissible=%s) | NON-CONFIRMATORY, DESCRIPTIVE |",
      r$domain_a, r$domain_b, r$k_a, r$k_b, r$n_both, r$n_complete,
      r$converged, r$admissible
    )
  }
}, character(1))

est <- out[out$estimable, ]
range_sentence <- if (nrow(est) > 0) {
  sprintf(
    paste0("Across the %d estimable descriptive pairs, uncorrected latent ",
           "ρ_level ranged %s–%s, uncorrected ρ_sens ranged %s–%s, and the ",
           "uncorrected gap ranged %s–%s."),
    nrow(est),
    fmt3(min(est$rho_level)), fmt3(max(est$rho_level)),
    fmt3(min(est$rho_sens_uncorrected)), fmt3(max(est$rho_sens_uncorrected)),
    fmt3(min(est$gap_uncorrected)), fmt3(max(est$gap_uncorrected))
  )
} else {
  "No descriptive pair was estimable."
}
not_est <- out[!out$estimable, ]
not_est_sentence <- if (nrow(not_est) > 0) {
  sprintf(" %d pair(s) were NOT ESTIMABLE (%s) — a reportable fragility finding, not a bug.",
          nrow(not_est), paste(not_est$pair, collapse = ", "))
} else {
  ""
}

md <- c(
  "# DESCRIPTIVE TIER — five remaining domain pairs",
  "",
  sprintf("**Date:** %s · **Code:** `analysis/07_descriptive_pairs.R` · **Data:** REALDATA, primary trial set, committed step-1 parcel columns AS-IS",
          format(Sys.Date())),
  "",
  "**Every row below except the marked confirmatory one is NON-CONFIRMATORY / DESCRIPTIVE.**",
  "These four 200/100 configurations and the 100/100 configuration FAILED the blind",
  "calibration gate in step 0b'' — so their sensitivity correlations and gaps are",
  "reported UNCORRECTED, with the frozen systematic bias band attached to each gap",
  "(±0.09 for 200/100 pairs, ±0.10 for Calories × GDP) instead of a calibration",
  "correction. The memory × vision row is the confirmatory result,",
  "copied verbatim from `results/step6_REALDATA_confirmatory.csv` — NOT recomputed.",
  "",
  "| Pair | K | n both / complete | ρ_level [95% CI] | ρ_sens [95% CI] | Gap [95% CI] ± band | Status |",
  "|---|---|---|---|---|---|---|",
  sprintf("| **memory × vision** | 4/4 | %d / %d | %s %s | %s %s (corrected) | %s %s (corrected, no band) | **CONFIRMATORY** |",
          c_n_both, c_gap$n,
          fmt3(c_lev$estimate), ci(c_lev$ci_lo, c_lev$ci_hi),
          fmt3(c_sen$estimate), ci(c_sen$ci_lo, c_sen$ci_hi),
          fmt3(c_gap$estimate), ci(c_gap$ci_lo, c_gap$ci_hi)),
  desc_rows,
  "",
  "## Honest reading (what is and is not claimable from this tier)",
  "",
  paste0(
    range_sentence, not_est_sentence,
    " The confidence intervals above capture SAMPLING uncertainty only; the ±band on ",
    "each descriptive gap is a separate, additive SYSTEMATIC allowance — the frozen ",
    "worst-case estimator bias observed in the 0b'' simulations at these trial counts, ",
    "where the calibration correction could not be validated. A descriptive gap whose ",
    "CI excludes zero is therefore still NOT a positive finding if it sits within its ",
    "band of zero (and even one that clears both is descriptive, not a test — nothing ",
    "in this tier was pre-specified as confirmatory). What this table legitimately ",
    "supports: (1) the cross-domain correlation LANDSCAPE for the paper's summary ",
    "table, with honest uncertainty in both channels; (2) a consistency check — do the ",
    "descriptive pairs point the same direction as the confirmatory memory × vision ",
    "result (gap = 0.113 [−0.138, 0.364], null not excluded)? What it does NOT ",
    "support: any claim that a specific descriptive gap is nonzero, any cross-pair ",
    "meta-analytic combination (the pairs share subjects and the bands are correlated ",
    "across pairs), or any confirmatory-style wording. The sensitivity correlations ",
    "here are UNCORRECTED and known from simulation to be attenuated at these trial ",
    "counts — they are lower bounds in expectation, not estimates of the latent truth."
  ),
  ""
)
writeLines(md, here("results", "step7_REALDATA_descriptive.md"))
message("\nDescriptive-tier outputs written (step7_REALDATA_descriptive.csv/.md).")
