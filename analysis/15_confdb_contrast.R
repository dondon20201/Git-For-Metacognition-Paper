library(here)
source(here("R", "confdb.R"))
source(here("R", "measures.R"))

n_boot <- 2000L
boot_seed <- 20260711L
min_trials <- 60L
raw_dir <- here("data", "raw", "confidence_db")

md5_expected <- c(
  data_Bang_2019_Exp2.csv     = "e697b4615d6f3da5fa7c10342d16dfe2",
  data_Haddara_2022_Expt1.csv = "fbe82cb5185b41d67894e429c8be1699"
)
for (fn in names(md5_expected)) {
  got <- unname(tools::md5sum(file.path(raw_dir, fn)))
  if (!identical(got, unname(md5_expected[fn]))) {
    stop("md5 mismatch for ", fn, ": got ", got, " expected ",
         md5_expected[fn], " — refusing to proceed on an unverified file")
  }
}
message("md5 guards passed: ", paste(names(md5_expected), collapse = ", "))

measures <- read.csv(here("results", "step11_CONFDB_measures.csv"),
                     stringsAsFactors = FALSE,
                     colClasses = c(task = "character",
                                    subj_idx = "character"))
reliab <- read.csv(here("results", "step11_CONFDB_reliability.csv"),
                   stringsAsFactors = FALSE,
                   colClasses = c(task = "character"))
stopifnot(
  nrow(measures) == 1525,
  all(c("study", "task", "subj_idx", "level", "auroc2", "dprime_hat",
        "n_trials_clean") %in% names(measures)),
  all(c("rel_sb", "ci_lo", "ci_hi") %in% names(reliab))
)

nullref <- read.csv(here("results", "step14_SIMULATION_confdb_nullref.csv"),
                    stringsAsFactors = FALSE)
nullref_meta <- jsonlite::read_json(
  here("results", "step14_SIMULATION_confdb_nullref_meta.json"))
nullref_pass <- identical(nullref_meta$gate, "PASS")
message("Null-reference gate: ", nullref_meta$gate)

dprime_confdb <- function(stim, resp, signal_value) {
  is_sig <- stim == signal_value
  n_s <- sum(is_sig)
  n_n <- sum(!is_sig)
  if (n_s == 0 || n_n == 0) return(NA_real_)
  h <- (sum(resp == signal_value & is_sig) + 0.5) / (n_s + 1)
  f <- (sum(resp == signal_value & !is_sig) + 0.5) / (n_n + 1)
  qnorm(h) - qnorm(f)
}

bang_parsed <- read_confdb_csv(file.path(raw_dir, "data_Bang_2019_Exp2.csv"))
stopifnot(bang_parsed$status %in% c("OK", "OK_WITH_FIX"))
bang_ct <- clean_trials(bang_parsed$data)
bang_cd <- bang_ct$data
catch_col <- find_col(bang_cd, "Catch_trial")
stopifnot(!is.na(catch_col))
catch_v <- suppressWarnings(as.numeric(bang_cd[[catch_col]]))
stopifnot(!anyNA(catch_v), all(catch_v %in% c(0, 1)))
bang_n_clean_before <- nrow(bang_cd)
bang_n_catch <- sum(catch_v == 1)
bang_cd <- bang_cd[catch_v != 1, , drop = FALSE]
bang_row1 <- measures[measures$study == "Bang_2019_Exp2", ][1, ]
stopifnot(min(bang_cd$.conf) == bang_row1$conf_scale_min,
          max(bang_cd$.conf) == bang_row1$conf_scale_max)
bang_cmin <- bang_row1$conf_scale_min
bang_cmax <- bang_row1$conf_scale_max
message(sprintf(
  "Bang recompute: %d clean trials, %d catch trials dropped (%.1f%%) -> %d",
  bang_n_clean_before, bang_n_catch,
  100 * bang_n_catch / bang_n_clean_before, nrow(bang_cd)))

s_col <- find_col(bang_cd, "Stimulus")
r_col <- find_col(bang_cd, "Response")
bang_rows <- list()
bi <- 0L
for (task in sort(unique(bang_cd$.task))) {
  td <- bang_cd[bang_cd$.task == task, , drop = FALSE]
  stim_vals <- sort(unique(trimws(as.character(td[[s_col]]))))
  stopifnot(length(stim_vals) == 2)
  cnt <- table(td$.subj)
  for (sj in names(cnt)[cnt >= min_trials]) {
    sd_ <- td[td$.subj == sj, , drop = FALSE]
    bi <- bi + 1L
    bang_rows[[bi]] <- data.frame(
      study = "Bang_2019_Exp2", task = task, subj_idx = sj,
      n_trials_clean = nrow(sd_),
      level = (mean(sd_$.conf) - bang_cmin) / (bang_cmax - bang_cmin),
      auroc2 = auroc2(sd_$.conf, sd_$.acc),
      dprime_hat = dprime_confdb(trimws(as.character(sd_[[s_col]])),
                                 trimws(as.character(sd_[[r_col]])),
                                 stim_vals[2]),
      stringsAsFactors = FALSE
    )
  }
}
bang_recomp <- do.call(rbind, bang_rows)
message("Bang recomputed subject-task rows: ", nrow(bang_recomp))

had_parsed <- read_confdb_csv(file.path(raw_dir, "data_Haddara_2022_Expt1.csv"))
stopifnot(had_parsed$status %in% c("OK", "OK_WITH_FIX"))
had_raw <- had_parsed$data
fb_col <- find_col(had_raw, "Feedback")
task_col <- find_task_col(had_raw)
subj_col <- find_subj_col(had_raw)
stopifnot(!is.na(fb_col), !is.na(task_col), !is.na(subj_col))
h1 <- had_raw[as.character(had_raw[[task_col]]) == "1", ]
fb_by_subj <- tapply(as.numeric(h1[[fb_col]]),
                     as.character(h1[[subj_col]]),
                     function(v) sort(unique(v[!is.na(v)])))
stopifnot(all(lengths(fb_by_subj) == 1))
h2 <- had_raw[as.character(had_raw[[task_col]]) == "2", ]
stopifnot(all(as.numeric(h2[[fb_col]]) == 0, na.rm = TRUE))
had_group <- ifelse(vapply(fb_by_subj, `[`, numeric(1), 1) == 1,
                    "feedback", "no_feedback")
message("Haddara groups: ", paste(names(table(had_group)),
                                  table(had_group), collapse = ", "))

make_wide <- function(m, study, task_a, task_b, subj_keep = NULL) {
  cols <- c("subj_idx", "level", "auroc2", "dprime_hat")
  a <- m[m$study == study & m$task == task_a, cols]
  b <- m[m$study == study & m$task == task_b, cols]
  w <- merge(a, b, by = "subj_idx", suffixes = c("_a", "_b"))
  if (!is.null(subj_keep)) w <- w[w$subj_idx %in% subj_keep, ]
  w[complete.cases(w), , drop = FALSE]
}

contrast_stats <- function(w) {
  res_a <- resid(lm(auroc2_a ~ dprime_hat_a, data = w))
  res_b <- resid(lm(auroc2_b ~ dprime_hat_b, data = w))
  r_level <- cor(w$level_a, w$level_b)
  r_sens <- cor(res_a, res_b)
  r_sens_raw <- cor(w$auroc2_a, w$auroc2_b)
  c(r_level = r_level, r_sens = r_sens, r_sens_raw = r_sens_raw,
    gap = r_level - r_sens, gap_raw = r_level - r_sens_raw,
    r_level_spearman = cor(w$level_a, w$level_b, method = "spearman"),
    r_sens_spearman = cor(res_a, res_b, method = "spearman"))
}

boot_contrast <- function(w) {
  set.seed(boot_seed)
  n <- nrow(w)
  point <- contrast_stats(w)
  boots <- replicate(n_boot, {
    idx <- sample.int(n, replace = TRUE)
    tryCatch(contrast_stats(w[idx, , drop = FALSE]),
             error = function(e) rep(NA_real_, length(point)))
  })
  n_bad <- sum(colSums(is.na(boots)) > 0)
  list(point = point,
       ci_lo = apply(boots, 1, quantile, 0.025, na.rm = TRUE, names = FALSE),
       ci_hi = apply(boots, 1, quantile, 0.975, na.rm = TRUE, names = FALSE),
       n = n, n_bad_boot = n_bad)
}

disatt_range <- function(r, rel_a, rel_b) {
  vals <- r / sqrt(rel_a * rel_b)
  range(vals)
}

rel_vec <- function(study, task, measure) {
  r <- reliab[reliab$study == study & reliab$task == task &
                reliab$measure == measure, ]
  stopifnot(nrow(r) == 1)
  c(r$rel_sb, r$ci_lo, r$ci_hi)
}

primary_defs <- list(
  list(study = "Arbuzova_2021_expt2", task_a = "1", task_b = "2",
       pair_label = "task 1 x task 2 (visuomotor, condition x condition)",
       data = measures),
  list(study = "Arbuzova_unpub_3", task_a = "memory", task_b = "visual",
       pair_label = "memory x visual (CROSS-DOMAIN)",
       data = measures),
  list(study = "Bang_2019_Exp2", task_a = "1", task_b = "2",
       pair_label = "task 1 x task 2 (coarse x fine Gabor, within vision)",
       data = bang_recomp),
  list(study = "Haddara_2022_Expt1", task_a = "1", task_b = "2",
       pair_label = "task 1 x task 2 (letter x color grids, within vision)",
       data = measures)
)

rows <- list()
ri <- 0L
add_row <- function(study, pair, subset, tier, quantity, est, lo, hi, n,
                    disatt_lo = NA_real_, disatt_hi = NA_real_,
                    nr30_lo = NA_real_, nr30_hi = NA_real_,
                    nr45_lo = NA_real_, nr45_hi = NA_real_, note = "") {
  ri <<- ri + 1L
  rows[[ri]] <<- data.frame(
    data_type = "REALDATA-EXTERNAL (CONFDB)", tier_label = "NON-CONFIRMATORY (robustness)",
    study = study, pair = pair, subset = subset, tier = tier,
    quantity = quantity, estimate = est, ci_lo = lo, ci_hi = hi,
    n_subjects = n, n_boot = n_boot, boot_seed = boot_seed,
    disatt_lo = disatt_lo, disatt_hi = disatt_hi,
    nullref_rho30_lo = nr30_lo, nullref_rho30_hi = nr30_hi,
    nullref_rho45_lo = nr45_lo, nullref_rho45_hi = nr45_hi,
    note = note, stringsAsFactors = FALSE)
}

primary_results <- list()
for (pd in primary_defs) {
  w <- make_wide(pd$data, pd$study, pd$task_a, pd$task_b)
  bc <- boot_contrast(w)
  message(sprintf("%s (%s): n=%d, r_level=%.3f, r_sens=%.3f, gap=%.3f%s",
                  pd$study, pd$pair_label, bc$n, bc$point["r_level"],
                  bc$point["r_sens"], bc$point["gap"],
                  if (bc$n_bad_boot > 0)
                    sprintf(" [%d bad boot reps]", bc$n_bad_boot) else ""))

  d_lev <- disatt_range(bc$point["r_level"],
                        rel_vec(pd$study, pd$task_a, "level"),
                        rel_vec(pd$study, pd$task_b, "level"))
  d_sen <- disatt_range(bc$point["r_sens"],
                        rel_vec(pd$study, pd$task_a, "auroc2"),
                        rel_vec(pd$study, pd$task_b, "auroc2"))
  nr <- nullref[nullref$study == pd$study, ]
  nr30 <- nr[abs(nr$rho_null - 0.30) < 1e-9, ]
  nr45 <- nr[abs(nr$rho_null - 0.45) < 1e-9, ]
  bang_note <- if (pd$study == "Bang_2019_Exp2")
    "recomputed from raw (md5-guarded), catch trials dropped, author-excluded subjects KEPT; disattenuation uses committed (catch-included) reliabilities" else ""

  q <- bc$point
  add_row(pd$study, pd$pair_label, "primary", "PRIMARY", "r_level",
          q["r_level"], bc$ci_lo[1], bc$ci_hi[1], bc$n,
          disatt_lo = d_lev[1], disatt_hi = d_lev[2], note = bang_note)
  add_row(pd$study, pd$pair_label, "primary", "PRIMARY", "r_sens",
          q["r_sens"], bc$ci_lo[2], bc$ci_hi[2], bc$n,
          disatt_lo = d_sen[1], disatt_hi = d_sen[2],
          note = paste("d'-partialed;", bang_note))
  add_row(pd$study, pd$pair_label, "primary", "SECONDARY", "r_sens_raw",
          q["r_sens_raw"], bc$ci_lo[3], bc$ci_hi[3], bc$n,
          note = "un-partialed (d' confound exhibit)")
  add_row(pd$study, pd$pair_label, "primary", "PRIMARY", "gap",
          q["gap"], bc$ci_lo[4], bc$ci_hi[4], bc$n,
          note = "r_level - r_sens (d'-partialed)")
  if (nullref_pass) {
    add_row(pd$study, pd$pair_label, "primary", "SECONDARY", "gap_raw",
            q["gap_raw"], bc$ci_lo[5], bc$ci_hi[5], bc$n,
            nr30_lo = nr30$gap_raw_null_lo, nr30_hi = nr30$gap_raw_null_hi,
            nr45_lo = nr45$gap_raw_null_lo, nr45_hi = nr45$gap_raw_null_hi,
            note = "raw gap; compare to SIMULATION null-reference interval before any interpretation")
  } else {
    add_row(pd$study, pd$pair_label, "primary", "SECONDARY", "gap_raw",
            q["gap_raw"], bc$ci_lo[5], bc$ci_hi[5], bc$n,
            note = "raw gap; null-reference gate FAILED - no null-reference columns; disattenuation ranges are the only yardstick")
  }
  add_row(pd$study, pd$pair_label, "primary", "SECONDARY",
          "r_level_spearman", q["r_level_spearman"], bc$ci_lo[6],
          bc$ci_hi[6], bc$n, note = "scale-robustness secondary")
  add_row(pd$study, pd$pair_label, "primary", "SECONDARY",
          "r_sens_spearman", q["r_sens_spearman"], bc$ci_lo[7],
          bc$ci_hi[7], bc$n, note = "scale-robustness secondary")

  primary_results[[pd$study]] <- list(
    bc = bc, d_lev = d_lev, d_sen = d_sen,
    nr30 = nr30, nr45 = nr45, pair_label = pd$pair_label)
}

run_secondary <- function(data, study, task_a, task_b, subset, note,
                          pair_label, subj_keep = NULL) {
  w <- make_wide(data, study, task_a, task_b, subj_keep)
  bc <- boot_contrast(w)
  q <- bc$point
  for (spec in list(c("r_level", 1), c("r_sens", 2), c("gap", 4))) {
    add_row(study, pair_label, subset, "SECONDARY", spec[1],
            q[spec[1]], bc$ci_lo[as.integer(spec[2])],
            bc$ci_hi[as.integer(spec[2])], bc$n, note = note)
  }
  message(sprintf("  [%s] %s: n=%d, r_level=%.3f, r_sens=%.3f, gap=%.3f",
                  subset, study, bc$n, q["r_level"], q["r_sens"], q["gap"]))
  bc
}

message("\nSecondaries:")
bs1 <- run_secondary(measures, "Bang_2019_Exp2", "1", "2", "B-s1",
                     "committed step11 measures AS-IS (catch trials included)",
                     primary_results[["Bang_2019_Exp2"]]$pair_label)
bang_keep <- unique(bang_recomp$subj_idx[
  !(suppressWarnings(as.numeric(bang_recomp$subj_idx)) %in% 179:201)])
bs2 <- run_secondary(bang_recomp, "Bang_2019_Exp2", "1", "2", "B-s2",
                     "recomputed (catch-dropped) WITHOUT author-excluded subjects 179-201",
                     primary_results[["Bang_2019_Exp2"]]$pair_label,
                     subj_keep = bang_keep)
hs1 <- list()
for (g in c("feedback", "no_feedback")) {
  keep <- names(had_group)[had_group == g]
  hs1[[g]] <- run_secondary(
    measures, "Haddara_2022_Expt1", "1", "2", paste0("H-s1_", g),
    sprintf("Task-1 %s group only (between-subject manipulation); heterogeneity check, no test", g),
    primary_results[["Haddara_2022_Expt1"]]$pair_label, subj_keep = keep)
}

message("\nExploratory Arbuzova_unpub_3 pairs:")
arb_tasks <- sort(unique(measures$task[measures$study == "Arbuzova_unpub_3"]))
stopifnot(setequal(arb_tasks, c("memory", "motor", "visual", "visuomotor")))
arb_expl <- list()
for (cmb in combn(arb_tasks, 2, simplify = FALSE)) {
  w <- make_wide(measures, "Arbuzova_unpub_3", cmb[1], cmb[2])
  bc <- boot_contrast(w)
  q <- bc$point
  pl <- paste(cmb[1], "x", cmb[2])
  for (spec in list(c("r_level", 1), c("r_sens", 2), c("gap", 4))) {
    add_row("Arbuzova_unpub_3", pl, "exploratory_pairs", "EXPLORATORY",
            spec[1], q[spec[1]], bc$ci_lo[as.integer(spec[2])],
            bc$ci_hi[as.integer(spec[2])], bc$n,
            note = "EXPLORATORY 6-pair table; motor/visuomotor are two conditions of one skittles paradigm")
  }
  arb_expl[[pl]] <- bc
  message(sprintf("  %s: n=%d, r_level=%.3f, r_sens=%.3f, gap=%.3f",
                  pl, bc$n, q["r_level"], q["r_sens"], q["gap"]))
}

out <- do.call(rbind, rows)
write.csv(out, here("results", "step15_CONFDB_contrast.csv"),
          row.names = FALSE)

verdicts <- lapply(primary_results, function(pr) {
  q <- pr$bc$point
  cond_a <- q["r_level"] > pr$d_sen[2]
  cond_b <- if (nullref_pass) {
    q["gap_raw"] > max(pr$nr30$gap_raw_null_hi, pr$nr45$gap_raw_null_hi)
  } else NA
  list(cond_a = unname(cond_a), cond_b = unname(cond_b),
       verdict = if (isTRUE(cond_a) && isTRUE(cond_b)) "holds"
                 else if (isTRUE(cond_a) || isTRUE(cond_b)) "partial"
                 else "fails")
})

fmt <- function(x) sprintf("%.3f", x)
fmt_ci <- function(e, lo, hi) sprintf("%.3f [%.3f, %.3f]", e, lo, hi)
fmt_rng <- function(lo, hi) sprintf("[%.3f, %.3f]", lo, hi)

md <- c(
  "# Step 5b — Confidence Database level-vs-sensitivity contrast",
  "",
  "**Data:** REALDATA-EXTERNAL (CONFDB), md5-verified at fetch and re-guarded at every raw read this run. **Status: VERIFIED** (every number below computed this run by `analysis/15_confdb_contrast.R`; saved in `results/step15_CONFDB_contrast.csv`).",
  "",
  "**TIER: CONFDB ROBUSTNESS, NON-CONFIRMATORY.** No calibration correction exists for these study configurations — the frozen Lund 200/200 calibration does not transfer — so nothing here is a confirmatory test. The yardsticks are (a) per-study SIMULATION null-reference intervals for the raw gap (`analysis/14_confdb_nullref.R`) and (b) disattenuation RANGES (never a corrected point value).",
  "",
  "## Taxonomy (binding for all write-ups)",
  "",
  "This tier tests **TASK-generalization**, with exactly one external cross-domain pair. Only **Arbuzova_unpub_3 memory x visual** is a cross-DOMAIN pair. Arbuzova_2021_expt2 is a visuomotor **condition x condition** pair; Bang_2019_Exp2 (coarse x fine Gabor) and Haddara_2022_Expt1 (letter x color grids) are cross-TASK pairs **within vision**. These four studies must never be presented as four cross-domain replications.",
  "",
  "## Primary results (pre-specified pairs; subject bootstrap, B = 2000, seed 20260711)",
  "",
  "| Study | Pair | n | r_level [95% CI] | r_sens (d'-partialed) [95% CI] | gap [95% CI] |",
  "|---|---|---|---|---|---|"
)
for (st in names(primary_results)) {
  pr <- primary_results[[st]]
  q <- pr$bc$point
  md <- c(md, sprintf("| %s | %s | %d | %s | %s | %s |",
                      st, pr$pair_label, pr$bc$n,
                      fmt_ci(q["r_level"], pr$bc$ci_lo[1], pr$bc$ci_hi[1]),
                      fmt_ci(q["r_sens"], pr$bc$ci_lo[2], pr$bc$ci_hi[2]),
                      fmt_ci(q["gap"], pr$bc$ci_lo[4], pr$bc$ci_hi[4])))
}
md <- c(md, "",
  "## Yardsticks per study (disattenuation ranges + SIMULATION null reference)",
  "",
  "| Study | r_level disatt. range | r_sens disatt. range | raw gap [95% CI] | null-ref interval (rho=0.30) | null-ref interval (rho=0.45) |",
  "|---|---|---|---|---|---|")
for (st in names(primary_results)) {
  pr <- primary_results[[st]]
  q <- pr$bc$point
  nr30_s <- if (nullref_pass) fmt_rng(pr$nr30$gap_raw_null_lo, pr$nr30$gap_raw_null_hi) else "SKIPPED (gate failed)"
  nr45_s <- if (nullref_pass) fmt_rng(pr$nr45$gap_raw_null_lo, pr$nr45$gap_raw_null_hi) else "SKIPPED (gate failed)"
  md <- c(md, sprintf("| %s | %s | %s | %s | %s | %s |",
                      st, fmt_rng(pr$d_lev[1], pr$d_lev[2]),
                      fmt_rng(pr$d_sen[1], pr$d_sen[2]),
                      fmt_ci(q["gap_raw"], pr$bc$ci_lo[5], pr$bc$ci_hi[5]),
                      nr30_s, nr45_s))
}
md <- c(md, "",
  "Disattenuation ranges are r / sqrt(rel_a x rel_b) evaluated at the committed step11 reliability point estimates and at both reliability CI bounds — a sensitivity RANGE, not a corrected estimate. The null-reference intervals are the 2.5–97.5% span of the RAW gap under a SIMULATED true-null world (true r_level = true r_sens = rho) matched to each study's n / trials / scale / reliabilities: they quantify how much raw-gap inflation differential attenuation alone produces. They are a simulation reference, not a correction and not a test.",
  "",
  "## Reading each study (pre-specified criteria)",
  "",
  "The Lund pattern **holds** in a study only if r_level exceeds the TOP of r_sens's disattenuated range AND the raw gap exceeds the null-reference interval under both null scenarios; one condition met = **partial**; neither = **fails**.",
  "")
for (st in names(primary_results)) {
  pr <- primary_results[[st]]
  v <- verdicts[[st]]
  q <- pr$bc$point
  cond_b_txt <- if (nullref_pass) {
    sprintf("raw gap %s the null-reference interval (%.3f vs upper bounds %.3f / %.3f)",
            if (isTRUE(v$cond_b)) "EXCEEDS" else "does NOT exceed",
            q["gap_raw"], pr$nr30$gap_raw_null_hi, pr$nr45$gap_raw_null_hi)
  } else "null-reference comparison SKIPPED (null-reference gate failed)"
  md <- c(md, sprintf(
    "- **%s: %s.** r_level = %s %s the top of r_sens's disattenuated range (%.3f); %s.",
    st, v$verdict, fmt(q["r_level"]),
    if (isTRUE(v$cond_a)) "EXCEEDS" else "does NOT exceed",
    pr$d_sen[2], cond_b_txt))
}
md <- c(md, "",
  "## Secondaries",
  "",
  "| Unit | Definition | n | r_level [95% CI] | r_sens [95% CI] | gap [95% CI] |",
  "|---|---|---|---|---|---|")
sec_units <- list(
  list(lab = "B-s1", def = "Bang, committed catch-included measures", bc = bs1),
  list(lab = "B-s2", def = "Bang, recomputed, without author-excluded 179-201", bc = bs2),
  list(lab = "H-s1 feedback", def = "Haddara, Task-1 feedback group", bc = hs1[["feedback"]]),
  list(lab = "H-s1 no-feedback", def = "Haddara, Task-1 no-feedback group", bc = hs1[["no_feedback"]])
)
for (u in sec_units) {
  q <- u$bc$point
  md <- c(md, sprintf("| %s | %s | %d | %s | %s | %s |",
                      u$lab, u$def, u$bc$n,
                      fmt_ci(q["r_level"], u$bc$ci_lo[1], u$bc$ci_hi[1]),
                      fmt_ci(q["r_sens"], u$bc$ci_lo[2], u$bc$ci_hi[2]),
                      fmt_ci(q["gap"], u$bc$ci_lo[4], u$bc$ci_hi[4])))
}
md <- c(md, "",
  "H-s1 is a heterogeneity check only — two correlations reported side by side, no test.",
  "",
  "### Scale-robustness (Spearman) and the d' confound exhibit",
  "",
  "| Study | Spearman r_level [95% CI] | Spearman r_sens [95% CI] | r_sens_raw (un-partialed) [95% CI] |",
  "|---|---|---|---|")
for (st in names(primary_results)) {
  pr <- primary_results[[st]]
  q <- pr$bc$point
  md <- c(md, sprintf("| %s | %s | %s | %s |",
                      st,
                      fmt_ci(q["r_level_spearman"], pr$bc$ci_lo[6], pr$bc$ci_hi[6]),
                      fmt_ci(q["r_sens_spearman"], pr$bc$ci_lo[7], pr$bc$ci_hi[7]),
                      fmt_ci(q["r_sens_raw"], pr$bc$ci_lo[3], pr$bc$ci_hi[3])))
}
md <- c(md, "",
  "## Exploratory: Arbuzova_unpub_3, all six pairs",
  "",
  "EXPLORATORY. motor and visuomotor are two conditions of one skittles paradigm (study readme), so only memory x visual is a clean cross-domain pair.",
  "",
  "| Pair | n | r_level [95% CI] | r_sens [95% CI] | gap [95% CI] |",
  "|---|---|---|---|---|")
for (pl in names(arb_expl)) {
  bc <- arb_expl[[pl]]
  q <- bc$point
  md <- c(md, sprintf("| %s | %d | %s | %s | %s |",
                      pl, bc$n,
                      fmt_ci(q["r_level"], bc$ci_lo[1], bc$ci_hi[1]),
                      fmt_ci(q["r_sens"], bc$ci_lo[2], bc$ci_hi[2]),
                      fmt_ci(q["gap"], bc$ci_lo[4], bc$ci_hi[4])))
}

n_holds <- sum(vapply(verdicts, function(v) v$verdict == "holds", logical(1)))
n_partial <- sum(vapply(verdicts, function(v) v$verdict == "partial", logical(1)))
n_fails <- sum(vapply(verdicts, function(v) v$verdict == "fails", logical(1)))
md <- c(md, "",
  "## Honest reading",
  "",
  paste0(
    "(i) This entire tier is **NON-CONFIRMATORY / ROBUSTNESS**: the four studies were selected mechanically, but no calibration correction exists for their configurations, so none of these estimates carries the confirmatory status of the Lund memory x vision result. ",
    "(ii) **Raw gaps are inflated by differential attenuation**: LEVEL is measured near-perfectly (committed split-half reliabilities 0.97–1.00) while AUROC2 is not (0.62–0.90), so a positive raw gap is expected even in a true-null world — each observed raw gap must be compared to its SIMULATION null-reference interval before any interpretation. ",
    sprintf("(iii) Per the pre-specified criteria, the Lund pattern **holds** in %d of 4 studies, is **partial** in %d, and **fails** in %d (study-by-study statements above). ",
            n_holds, n_partial, n_fails),
    "Where the pattern does not hold, a null is **not excluded** (it is never 'absent' — sensitivity here is weakly measured, not proven undifferentiated). ",
    "No meta-analytic pooling is performed across the four studies (heterogeneous tasks, scales, and populations; forest display only). ",
    "Because Bang and Haddara are within-vision task pairs and Arbuzova_2021 is a condition pair, **no claim of cross-domain generality is made from them**; the only external cross-domain evidence in this tier is the single Arbuzova_unpub_3 memory x visual pair."),
  "",
  "## Provenance / run info",
  "",
  sprintf("- Inputs: committed `results/step11_CONFDB_measures.csv` (asserted 1,525 rows) and `results/step11_CONFDB_reliability.csv`; raw files re-read this run with md5 guards: %s.",
          paste(sprintf("`%s` (%s)", names(md5_expected), md5_expected), collapse = ", ")),
  sprintf("- Bang recompute: %d clean trials, %d catch trials dropped (%.1f%%); %d subject-task rows; author-excluded subjects 179-201 KEPT in the primary (mechanical criteria only), removed in B-s2.",
          bang_n_clean_before, bang_n_catch,
          100 * bang_n_catch / bang_n_clean_before, nrow(bang_recomp)),
  sprintf("- Haddara feedback groups (from raw, between-subject, Task 1 only): %s.",
          paste(names(table(had_group)), as.integer(table(had_group)),
                sep = " = ", collapse = ", ")),
  sprintf("- Estimator: r_sens = Pearson correlation of AUROC2 residualized (OLS, within study-task, on the analysis sample) on the same-task d'-hat; residualization recomputed inside every bootstrap resample. Subjects: >= %d clean trials in BOTH tasks, listwise complete.",
          min_trials),
  sprintf("- Null-reference gate: %s.", nullref_meta$gate),
  sprintf("- %s; bootstrap B = %d, percentile, seed %d per analysis unit.",
          R.version.string, n_boot, boot_seed)
)
writeLines(md, here("results", "step15_CONFDB_contrast.md"))
message("\nWrote results/step15_CONFDB_contrast.csv and .md")
message("Verdicts: ", paste(names(verdicts),
        vapply(verdicts, `[[`, "", "verdict"), sep = "=", collapse = ", "))
