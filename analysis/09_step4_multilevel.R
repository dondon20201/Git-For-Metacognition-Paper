library(here)
source(here("R", "multilevel.R"))

t_start <- Sys.time()

gate <- read.csv(here("results", "step8_SIMULATION_step4gate_summary.csv"))
if (!all(gate$gate_overall == "PASS")) {
  stop("STEP 4 GATE DID NOT PASS — the multilevel decomposition must not run (pre-specified).")
}
stopifnot(identical(
  unname(tools::md5sum(here("data", "raw", "lund2025",
                            "metacognition_trialData_master.csv"))),
  "500b2e6e1eb0e47e717feba5b912ad9f"
))
input_md5 <- c(
  step1_REALDATA_subject_measures.csv = "97bdda9a331fac3baa59c20e96c4b06b",
  step8_SIMULATION_step4gate_summary.csv = "41c98f5703e7e178a1fc1ecf166e478a"
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

trials <- read.csv(here("data", "raw", "lund2025",
                        "metacognition_trialData_master.csv"))
trials <- trials[!is.na(trials$confidence) & trials$repeat_exclude == 0, ]
step1 <- read.csv(here("results", "step1_REALDATA_subject_measures.csv"))
sd_ok <- step1[step1$trial_set == "primary" & step1$qualifies, ]
keep <- paste(trials$subject, trials$modality) %in%
  paste(sd_ok$subject, sd_ok$domain)
message(sprintf("Trials: %d primary; %d in qualifying subject-domains (%d dropped with excluded subject-domains)",
                nrow(trials), sum(keep), sum(!keep)))
trials <- trials[keep, ]

dat <- data.frame(
  subject = factor(trials$subject),
  domain = factor(trials$modality,
                  levels = c("memory", "vision", "Calories", "GDP")),
  accuracy = trials$accuracy,
  conf = trials$confidence
)
dat$subject_domain <- interaction(dat$subject, dat$domain, drop = TRUE)
dat <- add_acc_dc(dat)
n_sub <- nlevels(dat$subject)
n_sd <- nlevels(dat$subject_domain)
message(sprintf("%d subjects, %d subject-domains, %d trials",
                n_sub, n_sd, nrow(dat)))

t0 <- Sys.time()
point <- step4_ladder(dat)
fit_secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
message(sprintf(
  "POINT (rung %s, conv %s, sing %s, %.0fs): ICC_level = %.3f, ICC_sens = %s, Delta = %s",
  point$rung, point$converged, point$singular, fit_secs,
  point$icc_level,
  ifelse(is.na(point$icc_sens), "NOT ESTIMABLE", sprintf("%.3f", point$icc_sens)),
  ifelse(is.na(point$delta_icc), "NOT ESTIMABLE", sprintf("%.3f", point$delta_icc))
))
if (!point$converged) {
  stop("Point-estimate fit did not converge even at M3 — NOT ESTIMABLE; report and stop (pre-specified).")
}

b_target <- 500L
projected_h <- function(B) B * fit_secs / 3600
B <- b_target
ci_resolution <- "FULL"
if (projected_h(B) > 2.5) B <- 200L
if (projected_h(B) > 2.5) {
  B <- 100L
  ci_resolution <- "LOW-RESOLUTION"
}
message(sprintf("Bootstrap: B = %d (projected %.1f h at %.0f s/fit) [%s]",
                B, projected_h(B), fit_secs, ci_resolution))

set.seed(20260713)
subs <- levels(dat$subject)
by_sub <- split(seq_len(nrow(dat)), dat$subject)
boot_rows <- vector("list", B)
n_escalated <- 0L
n_failed <- 0L
for (b in seq_len(B)) {
  draw <- sample(subs, length(subs), replace = TRUE)
  idx <- unlist(by_sub[draw], use.names = FALSE)
  bd <- dat[idx, ]
  bd$subject <- factor(rep(paste0("b", seq_along(draw)),
                           times = lengths(by_sub[draw])))
  bd$subject_domain <- interaction(bd$subject, bd$domain, drop = TRUE)
  bd <- add_acc_dc(bd)
  res <- step4_ladder(bd, start_rung = point$rung)
  ok <- res$converged && !res$singular
  if (!ok) n_failed <- n_failed + 1L
  if (res$rung != point$rung) n_escalated <- n_escalated + 1L
  boot_rows[[b]] <- data.frame(
    b = b, ok = ok, rung = res$rung,
    icc_level = res$icc_level, icc_sens = res$icc_sens,
    delta_icc = res$delta_icc,
    v_int_subject = res$v_int_subject,
    v_int_subject_domain = res$v_int_subject_domain,
    v_slope_subject = res$v_slope_subject,
    v_slope_subject_domain = res$v_slope_subject_domain,
    v_residual = res$v_residual
  )
  if (b %% 25 == 0) {
    message(sprintf("  bootstrap %d/%d (%.1f min elapsed; %d failed, %d escalated)",
                    b, B, as.numeric(difftime(Sys.time(), t_start, units = "mins")),
                    n_failed, n_escalated))
  }
}
boot <- do.call(rbind, boot_rows)
bok <- boot[boot$ok, ]
ci_stability <- ifelse(n_failed > 0.10 * B, "UNSTABLE", "STABLE")
pci <- function(x) quantile(x, c(0.025, 0.975), na.rm = TRUE, names = FALSE)

qty <- function(name, est, cis, note = "") {
  data.frame(
    data_type = "REALDATA",
    tier = "STEP4_MULTILEVEL (pre-specified, NON-CONFIRMATORY)",
    quantity = name, estimate = est, ci_lo = cis[1], ci_hi = cis[2],
    ci_method = sprintf("cluster bootstrap B=%d, percentile, %s, %s",
                        B, ci_resolution, ci_stability),
    rung = point$rung, n_subjects = n_sub, n_subject_domains = n_sd,
    n_trials = nrow(dat), n_boot_ok = nrow(bok),
    n_boot_failed = n_failed, n_boot_escalated = n_escalated,
    note = note
  )
}
out <- rbind(
  qty("v_int_subject", point$v_int_subject, pci(bok$v_int_subject),
      "person-general set-point variance"),
  qty("v_int_subject_domain", point$v_int_subject_domain,
      pci(bok$v_int_subject_domain), "domain-specific set-point variance"),
  qty("v_slope_subject", point$v_slope_subject, pci(bok$v_slope_subject),
      "person-general resolution-slope variance"),
  qty("v_slope_subject_domain", point$v_slope_subject_domain,
      pci(bok$v_slope_subject_domain),
      "domain-specific resolution-slope variance"),
  qty("v_residual", point$v_residual, pci(bok$v_residual), ""),
  qty("ICC_level", point$icc_level, pci(bok$icc_level),
      "person-general share of set-point variance"),
  qty("ICC_sens", point$icc_sens, pci(bok$icc_sens),
      "person-general share of resolution variance; NOT d'-partialed"),
  qty("Delta_ICC", point$delta_icc, pci(bok$delta_icc),
      "ICC_level - ICC_sens; variance-share estimand, NOT the Step-3 gap")
)
write.csv(out, here("results", "step9_REALDATA_step4_multilevel.csv"),
          row.names = FALSE)

dat_z <- dat
dm <- tapply(dat$conf, dat$domain, mean)
dsd <- tapply(dat$conf, dat$domain, sd)
dat_z$conf <- (dat$conf - dm[as.character(dat$domain)]) /
  dsd[as.character(dat$domain)]
r1 <- step4_ladder(dat_z)
message(sprintf("R1 (z within domain, rung %s): ICC_level = %.3f, ICC_sens = %s",
                r1$rung, r1$icc_level,
                ifelse(is.na(r1$icc_sens), "NOT ESTIMABLE",
                       sprintf("%.3f", r1$icc_sens))))

r2_status <- "SKIPPED (runtime)"
r2_icc_level <- NA_real_
r2_icc_sens <- NA_real_
r2 <- tryCatch({
  setTimeLimit(elapsed = 3600, transient = TRUE)
  on.exit(setTimeLimit(elapsed = Inf), add = TRUE)
  dato <- dat
  dato$conf_ord <- factor(dato$conf, levels = 1:7, ordered = TRUE)
  fit <- ordinal::clmm(
    conf_ord ~ domain + domain:acc_dc + (1 + acc_dc | subject) +
      (1 + acc_dc | subject_domain),
    data = dato, Hess = FALSE
  )
  vc <- lapply(ordinal::VarCorr(fit), function(m) diag(as.matrix(m)))
  list(fit = fit, vc = vc)
}, error = function(e) {
  message("R2 clmm attempt ended: ", conditionMessage(e))
  NULL
})
setTimeLimit(elapsed = Inf)
if (!is.null(r2)) {
  vs <- r2$vc$subject
  vsd <- r2$vc$subject_domain
  r2_icc_level <- vs["(Intercept)"] / (vs["(Intercept)"] + vsd["(Intercept)"])
  r2_icc_sens <- vs["acc_dc"] / (vs["acc_dc"] + vsd["acc_dc"])
  r2_status <- "COMPLETED"
  message(sprintf("R2 (clmm): ICC_level = %.3f, ICC_sens = %.3f",
                  r2_icc_level, r2_icc_sens))
}

blups <- function(fit, group, term) {
  re <- lme4::ranef(fit)
  hit <- NULL
  for (i in seq_along(re)) {
    if (names(re)[i] == group && term %in% colnames(re[[i]])) {
      hit <- setNames(re[[i]][[term]], rownames(re[[i]]))
    }
  }
  hit
}
slope_s <- blups(point$fit, "subject", "acc_dc")
slope_sd <- blups(point$fit, "subject_domain", "acc_dc")
diag_rows <- list()
if (!is.null(slope_s)) {
  dsub <- tapply(sd_ok$d_prime_hat, sd_ok$subject, mean)
  common <- intersect(names(slope_s), names(dsub))
  ct <- cor.test(slope_s[common], dsub[common])
  diag_rows[[1]] <- data.frame(
    block = "EXPLORATORY_DIAGNOSTIC", quantity = "r_subject_slopeBLUP_vs_mean_dprime",
    estimate = unname(ct$estimate), ci_lo = ct$conf.int[1],
    ci_hi = ct$conf.int[2], n = length(common),
    note = "LMM slope is NOT d'-partialed; confound quantified, not removed"
  )
}
if (!is.null(slope_sd)) {
  key <- paste(sd_ok$subject, sd_ok$domain, sep = ".")
  dmap <- setNames(sd_ok$d_prime_hat, key)
  common <- intersect(names(slope_sd), names(dmap))
  ct <- cor.test(slope_sd[common], dmap[common])
  diag_rows[[2]] <- data.frame(
    block = "EXPLORATORY_DIAGNOSTIC", quantity = "r_subjdom_slopeBLUP_vs_dprime",
    estimate = unname(ct$estimate), ci_lo = ct$conf.int[1],
    ci_hi = ct$conf.int[2], n = length(common),
    note = "matched step-1 d'-hat per subject-domain"
  )
}
diag_df <- rbind(
  data.frame(block = "SECONDARY_R1_zscored", quantity = "ICC_level",
             estimate = r1$icc_level, ci_lo = NA, ci_hi = NA, n = n_sub,
             note = sprintf("rung %s; point estimate only", r1$rung)),
  data.frame(block = "SECONDARY_R1_zscored", quantity = "ICC_sens",
             estimate = r1$icc_sens, ci_lo = NA, ci_hi = NA, n = n_sub,
             note = sprintf("rung %s; point estimate only", r1$rung)),
  data.frame(block = "SECONDARY_R2_ordinal", quantity = "ICC_level",
             estimate = r2_icc_level, ci_lo = NA, ci_hi = NA, n = n_sub,
             note = r2_status),
  data.frame(block = "SECONDARY_R2_ordinal", quantity = "ICC_sens",
             estimate = r2_icc_sens, ci_lo = NA, ci_hi = NA, n = n_sub,
             note = r2_status),
  do.call(rbind, diag_rows)
)
diag_df$data_type <- "REALDATA"
write.csv(diag_df, here("results", "step9_REALDATA_step4_diagnostics.csv"),
          row.names = FALSE)

ips <- ipsative_exhibit(dat)
ips$data_type <- "REALDATA"
ips$block <- "IPSATIVE_EXHIBIT"
ips$note <- "METHODS EXHIBIT ONLY — what the artifact-prone naive mean-subtraction analysis would have said; NOT evidence"
write.csv(ips, here("results", "step9_REALDATA_ipsative_exhibit.csv"),
          row.names = FALSE)
message(sprintf("Ipsative exhibit (real data): centered pairwise r range [%.3f, %.3f], mean %.3f",
                min(ips$r_centered), max(ips$r_centered), mean(ips$r_centered)))

fmt <- function(x) ifelse(is.na(x), "NOT ESTIMABLE", sprintf("%.2f", x))
d_est <- point$delta_icc
d_lo <- pci(bok$delta_icc)[1]
d_hi <- pci(bok$delta_icc)[2]
delta_sentence <- if (is.na(d_est)) {
  "Delta_ICC is NOT ESTIMABLE at the landed ladder rung."
} else if (d_lo > 0 || d_hi < 0) {
  sprintf("As a pre-specified secondary finding, stated with hedges: in this dataset and under this model, the person-general variance share was %s for the set-point and %s for resolution (Delta_ICC = %s, 95%% bootstrap CI [%s, %s], CI excluding zero). This is a variance-share result, not a re-test of the confirmatory latent-correlation gap.",
          fmt(point$icc_level), fmt(point$icc_sens), fmt(d_est), fmt(d_lo), fmt(d_hi))
} else {
  sprintf("A Delta_ICC of zero is not excluded (Delta_ICC = %s, 95%% bootstrap CI [%s, %s]).",
          fmt(d_est), fmt(d_lo), fmt(d_hi))
}

md <- c(
  "# STEP 4 — Non-ipsative multilevel decomposition",
  "",
  "**Tier: STEP4_MULTILEVEL (pre-specified, NON-CONFIRMATORY).**",
  "",
  sprintf("**Date:** %s · **Code:** `analysis/09_step4_multilevel.R` · **Data:** REALDATA, primary trial set, step-1 qualifying subject-domains (n = %d subjects, %d subject-domains, %d trials)",
          format(Sys.Date()), n_sub, n_sd, nrow(dat)),
  "",
  "**Estimand note (binding):** the ICCs below are VARIANCE-SHARE",
  "generalizability indices from a trial-level mixed model — a different",
  "estimand from Step 3's latent cross-domain correlations. The random",
  "accuracy slope measures confidence RESOLUTION (Δconf correct−incorrect),",
  "not AUROC2, and is NOT d′-partialed (the d′ confound is quantified in the",
  "diagnostics, not removed). Nothing here confirms, disconfirms, or re-tests",
  "the 006 confirmatory result; the 'null not excluded' wording of 006 stands.",
  "",
  "## Decomposition",
  "",
  sprintf("Ladder rung: **%s** (converged %s, singular %s). Bootstrap: B = %d, %d ok / %d failed / %d escalated (%s, %s).",
          point$rung, point$converged, point$singular, B, nrow(bok),
          n_failed, n_escalated, ci_resolution, ci_stability),
  "",
  "| Quantity | Estimate | 95% CI (cluster bootstrap) | Status |",
  "|---|---|---|---|",
  sprintf("| ICC_level (person-general share of set-point variance) | %s | [%s, %s] | VERIFIED, NON-CONFIRMATORY |",
          fmt(point$icc_level), fmt(pci(bok$icc_level)[1]), fmt(pci(bok$icc_level)[2])),
  sprintf("| ICC_sens (person-general share of resolution variance) | %s | [%s, %s] | VERIFIED, NON-CONFIRMATORY |",
          fmt(point$icc_sens), fmt(pci(bok$icc_sens)[1]), fmt(pci(bok$icc_sens)[2])),
  sprintf("| Δ_ICC | %s | [%s, %s] | VERIFIED, NON-CONFIRMATORY |",
          fmt(d_est), fmt(d_lo), fmt(d_hi)),
  "",
  delta_sentence,
  "",
  "## Secondary / robustness",
  "",
  sprintf("- R1 (confidence z-scored within domain, rung %s): ICC_level = %s, ICC_sens = %s — %s the raw-scale decomposition.",
          r1$rung, fmt(r1$icc_level), fmt(r1$icc_sens),
          ifelse(!is.na(r1$icc_level) && !is.na(point$icc_level) &&
                   abs(r1$icc_level - point$icc_level) < 0.1,
                 "consistent with", "compare against")),
  sprintf("- R2 (ordinal cumulative-link mixed model): %s%s", r2_status,
          ifelse(r2_status == "COMPLETED",
                 sprintf(" — ICC_level = %s, ICC_sens = %s",
                         fmt(r2_icc_level), fmt(r2_icc_sens)), "")),
  sprintf("- Exploratory diagnostics: see `step9_REALDATA_step4_diagnostics.csv` — the resolution-slope BLUPs' correlation with d′-hat quantifies the standing d′ confound (the LMM slope is not d′-partialed)."),
  "",
  "## Ipsative exhibit (METHODS EXHIBIT ONLY — never evidence)",
  "",
  sprintf("The naive mean-subtraction analysis on the same trials yields centered domain-mean correlations from %.3f to %.3f (mean %.3f) across the six pairs — the mechanical negativity the non-ipsative model exists to avoid (cf. the simulation exhibit, where the same computation manufactures ≈ −1/3 under a planted null structure while the LMM recovers the planted shares).",
          min(ips$r_centered), max(ips$r_centered), mean(ips$r_centered)),
  "",
  "## Honest reading",
  "",
  "What this tier supports: a decomposition of confidence variance into",
  "person-general vs person×domain shares for set-point and resolution,",
  "with cluster-bootstrap uncertainty, free of the ipsative artifact.",
  "What it does NOT support: any claim about the Step-3 latent gap, any",
  "confirmatory-style wording, or treating the resolution slope as",
  "metacognitive sensitivity in the AUROC2 sense. The resolution slope",
  "carries a known, quantified d′ confound and its domain-specific share is",
  "partly measurement noise at 100-trial domains.",
  ""
)
writeLines(md, here("results", "step9_REALDATA_step4_multilevel.md"))

meta <- list(
  data_type = "REALDATA",
  tier = "STEP4_MULTILEVEL (pre-specified, NON-CONFIRMATORY)",
  gate = "PASSED (results/step8_SIMULATION_step4gate_summary.csv)",
  n_subjects = n_sub, n_subject_domains = n_sd, n_trials = nrow(dat),
  point_rung = point$rung, point_fit_seconds = fit_secs,
  bootstrap = list(B = B, seed = 20260713, n_ok = nrow(bok),
                   n_failed = n_failed, n_escalated = n_escalated,
                   resolution = ci_resolution, stability = ci_stability),
  r2_status = r2_status,
  total_minutes = as.numeric(difftime(Sys.time(), t_start, units = "mins")),
  r_version = R.version.string,
  lme4_version = as.character(packageVersion("lme4"))
)
jsonlite::write_json(meta, here("results", "step9_REALDATA_step4_meta.json"),
                     auto_unbox = TRUE, digits = 10, pretty = TRUE)
message(sprintf("\nStep 4 outputs written. Total %.1f min.", meta$total_minutes))
