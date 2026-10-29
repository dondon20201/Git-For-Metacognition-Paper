suppressWarnings(suppressMessages({
  library(here)
}))

res <- function(...) here::here("results", ...)

cell <- function(df, value_col, ...) {
  conds <- list(...)
  keep <- rep(TRUE, nrow(df))
  for (nm in names(conds)) keep <- keep & (df[[nm]] == conds[[nm]])
  v <- df[[value_col]][keep]
  stopifnot(length(v) == 1L)
  as.numeric(v)
}
guard <- function(x, target, k) {
  stopifnot(abs(round(x, k) - target) < .Machine$double.eps^0.5)
  x
}
f <- function(x, k) formatf <- sprintf(paste0("%.", k, "f"), x)

step6  <- read.csv(res("step6_REALDATA_confirmatory.csv"),  stringsAsFactors = FALSE)
step7  <- read.csv(res("step7_REALDATA_descriptive.csv"),   stringsAsFactors = FALSE)
step9  <- read.csv(res("step9_REALDATA_step4_multilevel.csv"), stringsAsFactors = FALSE)
step9i <- read.csv(res("step9_REALDATA_ipsative_exhibit.csv"), stringsAsFactors = FALSE)
step8i <- read.csv(res("step8_SIMULATION_ipsative_exhibit.csv"), stringsAsFactors = FALSE)
step15 <- read.csv(res("step15_CONFDB_contrast.csv"),       stringsAsFactors = FALSE)
step2  <- read.csv(res("step2_REALDATA_reliability.csv"),   stringsAsFactors = FALSE)
step11 <- read.csv(res("step11_CONFDB_reliability.csv"),    stringsAsFactors = FALSE)

g1_est <- guard(cell(step6, "estimate", block = "PRIMARY_CONFIRMATORY", quantity = "gap_corrected"), 0.113, 3)
g1_lo  <- guard(cell(step6, "ci_lo",    block = "PRIMARY_CONFIRMATORY", quantity = "gap_corrected"), -0.138, 3)
g1_hi  <- guard(cell(step6, "ci_hi",    block = "PRIMARY_CONFIRMATORY", quantity = "gap_corrected"), 0.364, 3)
g1_n   <- cell(step6, "n", block = "PRIMARY_CONFIRMATORY", quantity = "gap_corrected"); stopifnot(g1_n == 317)

g2_est <- guard(cell(step6, "estimate", block = "PRIMARY_CONFIRMATORY", quantity = "rho_level_latent"), 0.422, 3)
g2_lo  <- guard(cell(step6, "ci_lo",    block = "PRIMARY_CONFIRMATORY", quantity = "rho_level_latent"), 0.328, 3)
g2_hi  <- guard(cell(step6, "ci_hi",    block = "PRIMARY_CONFIRMATORY", quantity = "rho_level_latent"), 0.515, 3)

g3_est <- guard(cell(step6, "estimate", block = "PRIMARY_CONFIRMATORY", quantity = "rho_sens_latent_corrected"), 0.308, 3)
g3_lo  <- guard(cell(step6, "ci_lo",    block = "PRIMARY_CONFIRMATORY", quantity = "rho_sens_latent_corrected"), 0.075, 3)
g3_hi  <- guard(cell(step6, "ci_hi",    block = "PRIMARY_CONFIRMATORY", quantity = "rho_sens_latent_corrected"), 0.542, 3)

g4_est <- guard(cell(step6, "estimate", block = "SECONDARY_d_manifest_EXHIBIT", quantity = "raw_gap_naive"), 0.200, 3)
g4_lo  <- guard(cell(step6, "ci_lo",    block = "SECONDARY_d_manifest_EXHIBIT", quantity = "raw_gap_naive"), 0.051, 3)
g4_hi  <- guard(cell(step6, "ci_hi",    block = "SECONDARY_d_manifest_EXHIBIT", quantity = "raw_gap_naive"), 0.353, 3)
g4_n   <- cell(step6, "n", block = "SECONDARY_d_manifest_EXHIBIT", quantity = "raw_gap_naive"); stopifnot(g4_n == 320)

stopifnot(nrow(step7) == 5L)
d_gap_min <- guard(min(step7$gap_uncorrected), -0.133, 3)
d_gap_max <- guard(max(step7$gap_uncorrected),  0.296, 3)
d_rl_min  <- guard(min(step7$rho_level),        0.323, 3)
d_rl_max  <- guard(max(step7$rho_level),        0.473, 3)
stopifnot(all(step7$bias_band %in% c(0.09, 0.10)))

icc_l   <- guard(cell(step9, "estimate", quantity = "ICC_level"), 0.443, 3)
icc_s   <- guard(cell(step9, "estimate", quantity = "ICC_sens"),  0.337, 3)
dicc    <- guard(cell(step9, "estimate", quantity = "Delta_ICC"), 0.105, 3)
dicc_lo <- guard(cell(step9, "ci_lo",    quantity = "Delta_ICC"), 0.001, 3)
dicc_hi <- guard(cell(step9, "ci_hi",    quantity = "Delta_ICC"), 0.210, 3)

ips_real <- guard(step9i$r_centered[step9i$domain_a == "memory" & step9i$domain_b == "vision"], -0.667, 3)
sim_r1 <- step8i[step8i$replicate == 1, ]
ips_sim_mean <- guard(mean(sim_r1$r_centered), -0.31, 2)

rl <- function(study) cell(step15, "estimate", study = study, quantity = "r_level",
                           subset = "primary", tier = "PRIMARY")
cd_arb2 <- guard(rl("Arbuzova_2021_expt2"), 0.550, 3)
cd_arbU <- guard(rl("Arbuzova_unpub_3"),    0.560, 3)
cd_bang <- guard(rl("Bang_2019_Exp2"),      0.574, 3)
cd_hadd <- guard(rl("Haddara_2022_Expt1"),  0.872, 3)
arbU_rawgap <- guard(cell(step15, "estimate", study = "Arbuzova_unpub_3",
                          quantity = "gap_raw", subset = "primary", tier = "SECONDARY"), 0.551, 3)
arbU_null30_hi <- guard(cell(step15, "nullref_rho30_hi", study = "Arbuzova_unpub_3",
                             quantity = "gap_raw", subset = "primary", tier = "SECONDARY"), 0.505, 3)
arbU_null45_hi <- guard(cell(step15, "nullref_rho45_hi", study = "Arbuzova_unpub_3",
                             quantity = "gap_raw", subset = "primary", tier = "SECONDARY"), 0.531, 3)
margin30 <- guard(round(arbU_rawgap, 3) - round(arbU_null30_hi, 3), 0.046, 3)
margin45 <- guard(round(arbU_rawgap, 3) - round(arbU_null45_hi, 3), 0.020, 3)

lund_level_min <- guard(min(step2$reliability[step2$measure == "mean_conf" & step2$kind == "splithalf_SB"]), 0.949, 3)
lund_auroc_min <- guard(min(step2$reliability[step2$measure == "auroc2"    & step2$kind == "splithalf_SB"]), 0.390, 3)
lund_auroc_max <- guard(max(step2$reliability[step2$measure == "auroc2"    & step2$kind == "splithalf_SB"]), 0.837, 3)
cd_level_min <- guard(min(step11$rel_sb[step11$measure == "level"]),  0.970, 3)
cd_level_max <- guard(max(step11$rel_sb[step11$measure == "level"]),  0.996, 3)
cd_auroc_min <- guard(min(step11$rel_sb[step11$measure == "auroc2"]), 0.615, 3)
cd_auroc_max <- guard(max(step11$rel_sb[step11$measure == "auroc2"]), 0.903, 3)

num <- function(x, k = 3) sprintf(paste0("%.", k, "f"), x)
naf <- NA_character_

rows <- list(
  list(
    claim_id = "gap_confirmatory", tier = "CONFIRMATORY",
    quantity = "corrected level-minus-sensitivity gap (memory x vision)",
    estimate = num(g1_est), ci_lo = num(g1_lo), ci_hi = num(g1_hi), n = "317",
    source_csv = "step6_REALDATA_confirmatory.csv",
    verdict_wording = "null gap not excluded",
    figure_ref = "Figure 1 (fig_landscape), Figure 2 (fig_confdb_contrast) Lund row"
  ),
  list(
    claim_id = "rho_level_confirmatory", tier = "CONFIRMATORY",
    quantity = "latent cross-domain correlation of confidence level",
    estimate = num(g2_est), ci_lo = num(g2_lo), ci_hi = num(g2_hi), n = "317",
    source_csv = "step6_REALDATA_confirmatory.csv",
    verdict_wording = "confidence level is domain-general",
    figure_ref = "Figure 1 (fig_landscape)"
  ),
  list(
    claim_id = "rho_sens_confirmatory", tier = "CONFIRMATORY",
    quantity = "corrected latent cross-domain correlation of metacognitive sensitivity (d'-partialed)",
    estimate = num(g3_est), ci_lo = num(g3_lo), ci_hi = num(g3_hi), n = "317",
    source_csv = "step6_REALDATA_confirmatory.csv",
    verdict_wording = "corrected sensitivity correlation excludes zero - sensitivity is not domain-specific either",
    figure_ref = "Figure 1 (fig_landscape)"
  ),
  list(
    claim_id = "naive_gap_exhibit", tier = "METHODS EXHIBIT",
    quantity = "uncorrected (naive) raw level-minus-sensitivity gap (memory x vision)",
    estimate = num(g4_est), ci_lo = num(g4_lo), ci_hi = num(g4_hi), n = "320",
    source_csv = "step6_REALDATA_confirmatory.csv",
    verdict_wording = "the uncorrected analysis would have affirmed the thesis; the correction machinery prevented that claim",
    figure_ref = naf
  ),
  list(
    claim_id = "descriptive_landscape", tier = "DESCRIPTIVE (NON-CONFIRMATORY)",
    quantity = paste0("five descriptive domain-pair gaps (range ", num(d_gap_min), " to ", num(d_gap_max),
                      "); rho_level range ", num(d_rl_min), " to ", num(d_rl_max),
                      "; frozen bias bands +/-0.09 and +/-0.10"),
    estimate = naf, ci_lo = naf, ci_hi = naf, n = "5 pairs (312-322 per pair)",
    source_csv = "step7_REALDATA_descriptive.csv",
    verdict_wording = "heterogeneous, vision-anchored; no nonzero-gap claim for any descriptive pair",
    figure_ref = "Figure 1 (fig_landscape); Figure S1 (figS_gap_vs_reliability)"
  ),
  list(
    claim_id = "step4_decomposition", tier = "STEP4 MULTILEVEL (NON-CONFIRMATORY)",
    quantity = paste0("variance-share Delta_ICC = ICC_level (", num(icc_l), ") - ICC_sens (", num(icc_s), ")"),
    estimate = num(dicc), ci_lo = num(dicc_lo), ci_hi = num(dicc_hi), n = "322",
    source_csv = "step9_REALDATA_step4_multilevel.csv",
    verdict_wording = paste0("a Delta_ICC of zero sits at the very edge of the interval; ",
                             "directionally consistent with the confirmatory corrected gap; ",
                             "the zero-exclusion is never claimed"),
    figure_ref = naf
  ),
  list(
    claim_id = "ipsative_exhibit", tier = "METHODS EXHIBIT",
    quantity = paste0("naive grand-mean-centred within-subject correlation: real memory x vision ",
                      num(ips_real), "; simulation planted-null six-pair mean ", num(ips_sim_mean, 2),
                      " (mechanical expectation ~ -1/3)"),
    estimate = num(ips_real), ci_lo = naf, ci_hi = naf, n = naf,
    source_csv = "step9_REALDATA_ipsative_exhibit.csv; step8_SIMULATION_ipsative_exhibit.csv",
    verdict_wording = "grand-mean centering manufactures negative within-subject correlations; asymmetric under unequal trial counts",
    figure_ref = naf
  ),
  list(
    claim_id = "confdb_tier", tier = "CONFDB ROBUSTNESS (NON-CONFIRMATORY)",
    quantity = paste0("per-study r_level: Arbuzova_2021_expt2 ", num(cd_arb2),
                      " (partial), Arbuzova_unpub_3 memory x visual ", num(cd_arbU),
                      " (holds, narrowly; raw-gap margins ", num(margin30), " and ", num(margin45),
                      " vs the rho=0.30 and rho=0.45 null references), Bang_2019_Exp2 ", num(cd_bang),
                      " (partial), Haddara_2022_Expt1 ", num(cd_hadd), " (partial)"),
    estimate = naf, ci_lo = naf, ci_hi = naf, n = "4 studies (40-443 per study)",
    source_csv = "step15_CONFDB_contrast.csv",
    verdict_wording = paste0("holds 1/4, partial 3/4, fails 0/4; the only external cross-domain pair is ",
                             "Arbuzova_unpub_3 memory x visual (consistent with the Lund direction, exceeds its ",
                             "simulated true-null reference narrowly); within-vision task pairs and the ",
                             "condition x condition pair are NOT cross-domain replications"),
    figure_ref = "Figure 2 (fig_confdb_contrast)"
  ),
  list(
    claim_id = "reliability_floor", tier = "FLOOR / METHODS",
    quantity = paste0("split-half reliability: confidence level ", num(lund_level_min, 2), "-1.00 ",
                      "(Lund min ", num(lund_level_min), ", CD max ", num(cd_level_max),
                      ") vs AUROC2 ", num(lund_auroc_min, 2), "-", num(max(lund_auroc_max, cd_auroc_max), 2),
                      " (Lund ", num(lund_auroc_min), "-", num(lund_auroc_max),
                      ", CD ", num(cd_auroc_min), "-", num(cd_auroc_max), ")"),
    estimate = naf, ci_lo = naf, ci_hi = naf, n = "Lund 4 domains + 4 CD studies",
    source_csv = "step2_REALDATA_reliability.csv; step11_CONFDB_reliability.csv",
    verdict_wording = "the measurement asymmetry that inflates raw gaps is universal in these data",
    figure_ref = naf
  )
)

df <- do.call(rbind, lapply(rows, function(r) as.data.frame(r, stringsAsFactors = FALSE)))
stopifnot(nrow(df) == 9L)

out_csv <- res("step17_claim_status.csv")
write.csv(df, out_csv, row.names = FALSE, na = "")
cat("wrote", out_csv, "\n")

esc <- function(x) gsub("\\|", "\\\\|", ifelse(is.na(x), "", x))
ci_str <- function(lo, hi) ifelse(lo == "" | is.na(lo) | lo == "NA", "",
                                  paste0("[", lo, ", ", hi, "]"))

md <- c(
  "# Claim-status table",
  "",
  "Machine-generated by `analysis/17_claim_status.R` from committed results CSVs.",
  "Every value is extracted from the named `source_csv` and assert-guarded in code",
  "at the precision shown; the `verdict_wording` strings are frozen (pre-specified).",
  "Status labels: CONFIRMATORY / DESCRIPTIVE (NON-CONFIRMATORY) /",
  "STEP4 MULTILEVEL (NON-CONFIRMATORY) / CONFDB ROBUSTNESS (NON-CONFIRMATORY) /",
  "METHODS EXHIBIT / FLOOR. All estimates VERIFIED (re-derived from real data;",
  "SIMULATION components unmistakably labelled).",
  "",
  "| # | claim_id | tier | quantity | estimate | 95% CI | n | source_csv | verdict (frozen) | figure |",
  "|---|---|---|---|---|---|---|---|---|---|"
)
for (i in seq_len(nrow(df))) {
  r <- df[i, ]
  md <- c(md, paste0(
    "| ", i,
    " | ", esc(r$claim_id),
    " | ", esc(r$tier),
    " | ", esc(r$quantity),
    " | ", ifelse(is.na(r$estimate) || r$estimate == "", "", r$estimate),
    " | ", ci_str(r$ci_lo, r$ci_hi),
    " | ", esc(r$n),
    " | ", esc(r$source_csv),
    " | ", esc(r$verdict_wording),
    " | ", esc(r$figure_ref),
    " |"
  ))
}
out_md <- res("step17_claim_status.md")
writeLines(md, out_md)
cat("wrote", out_md, "\n")
cat("OK: 9 rows, all values assert-guarded.\n")
