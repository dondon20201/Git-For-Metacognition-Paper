library(here)

con <- read.csv(here("results", "step15_CONFDB_contrast.csv"),
                stringsAsFactors = FALSE)
nul <- read.csv(here("results", "step14_SIMULATION_confdb_nullref.csv"),
                stringsAsFactors = FALSE)
lund <- read.csv(here("results", "step6_REALDATA_confirmatory.csv"),
                 stringsAsFactors = FALSE)

lund_get <- function(qty) {
  r <- lund[lund$block == "PRIMARY_CONFIRMATORY" & lund$quantity == qty, ]
  stopifnot(nrow(r) == 1)
  c(est = r$estimate, lo = r$ci_lo, hi = r$ci_hi)
}
l_lev <- lund_get("rho_level_latent")
l_sen <- lund_get("rho_sens_latent_corrected")
l_gap <- lund_get("gap_corrected")
stopifnot(round(l_lev["est"], 3) == 0.422,
          round(l_sen["est"], 3) == 0.308,
          round(l_gap["est"], 3) == 0.113)

cd_studies <- c("Arbuzova_2021_expt2", "Arbuzova_unpub_3", "Bang_2019_Exp2",
                "Haddara_2022_Expt1")
cd_get <- function(study, qty) {
  r <- con[con$study == study & con$subset == "primary" &
             con$quantity == qty, ]
  stopifnot(nrow(r) == 1)
  r
}
nul_get <- function(study, rho) {
  r <- nul[nul$study == study & abs(nul$rho_null - rho) < 1e-9, ]
  stopifnot(nrow(r) == 1)
  r
}

row_labels <- c(
  Lund = "Lund memory x vision\nCONFIRMATORY (corrected)",
  Arbuzova_2021_expt2 = "Arbuzova 2021 E2, task 1 x 2\n(visuomotor cond x cond) NON-CONF.",
  Arbuzova_unpub_3 = "Arbuzova unpub 3, memory x visual\n(CROSS-DOMAIN) NON-CONF.",
  Bang_2019_Exp2 = "Bang 2019 E2, coarse x fine Gabor\n(within vision) NON-CONF.",
  Haddara_2022_Expt1 = "Haddara 2022 E1, letter x color\n(within vision) NON-CONF."
)
ys <- c(Lund = 5, Arbuzova_2021_expt2 = 4, Arbuzova_unpub_3 = 3,
        Bang_2019_Exp2 = 2, Haddara_2022_Expt1 = 1)

pd <- list()
pi_ <- 0L
add_pd <- function(row_id, panel, quantity, est, lo, hi, style, note = "") {
  pi_ <<- pi_ + 1L
  pd[[pi_]] <<- data.frame(row_id = row_id, y = ys[[row_id]], panel = panel,
                           quantity = quantity, estimate = est, lo = lo,
                           hi = hi, style = style, note = note,
                           stringsAsFactors = FALSE)
}
add_pd("Lund", "correlations", "r_level", l_lev["est"], l_lev["lo"],
       l_lev["hi"], "point_filled", "REALDATA latent, delta-method CI")
add_pd("Lund", "correlations", "r_sens_corrected", l_sen["est"], l_sen["lo"],
       l_sen["hi"], "point_open", "calibration-corrected, CONFIRMATORY")
add_pd("Lund", "gap", "gap_corrected", l_gap["est"], l_gap["lo"], l_gap["hi"],
       "diamond_filled", "corrected gap; null not excluded")
for (st in cd_studies) {
  rl <- cd_get(st, "r_level")
  rs <- cd_get(st, "r_sens")
  gr <- cd_get(st, "gap_raw")
  add_pd(st, "correlations", "r_level", rl$estimate, rl$ci_lo, rl$ci_hi,
         "point_filled", "subject bootstrap CI")
  add_pd(st, "correlations", "r_sens", rs$estimate, rs$ci_lo, rs$ci_hi,
         "point_open", "d'-partialed; subject bootstrap CI")
  add_pd(st, "correlations", "r_sens_disatt_range", NA, rs$disatt_lo,
         rs$disatt_hi, "bracket", "disattenuation RANGE, not a corrected point")
  add_pd(st, "gap", "gap_raw", gr$estimate, gr$ci_lo, gr$ci_hi,
         "diamond_open", "raw gap (r_level - r_sens_raw)")
  n30 <- nul_get(st, 0.30)
  n45 <- nul_get(st, 0.45)
  add_pd(st, "gap", "nullref_rho30", NA, n30$gap_raw_null_lo,
         n30$gap_raw_null_hi, "nullref_dotted",
         "SIMULATION null reference (true r_level = r_sens = 0.30)")
  add_pd(st, "gap", "nullref_rho45", NA, n45$gap_raw_null_lo,
         n45$gap_raw_null_hi, "nullref_dashed",
         "SIMULATION null reference (true r_level = r_sens = 0.45)")
}
pdat <- do.call(rbind, pd)
write.csv(pdat, here("results", "fig_confdb_contrast_data.csv"),
          row.names = FALSE)

draw_fig <- function() {
  layout(matrix(1:2, nrow = 1), widths = c(1.45, 1))
  xlim_l <- c(-0.45, 1.0)
  xlim_r <- c(-0.45, 1.0)
  gray_ref <- "gray45"

  par(mar = c(4.2, 12.5, 3, 0.8))
  plot(NA, xlim = xlim_l, ylim = c(-0.35, 5.6), yaxt = "n", xlab = "",
       ylab = "", main = "Cross-task correlation", cex.main = 1)
  abline(v = 0, col = "gray70", lty = 3)
  axis(2, at = ys, labels = row_labels[names(ys)], las = 1, cex.axis = 0.62,
       tick = FALSE)
  mtext("Pearson r", side = 1, line = 2.4, cex = 0.8)
  for (rid in names(ys)) {
    y <- ys[[rid]]
    sub <- pdat[pdat$row_id == rid & pdat$panel == "correlations", ]
    lv <- sub[sub$quantity == "r_level", ]
    sn <- sub[grepl("^r_sens($|_corrected)", sub$quantity), ]
    br <- sub[sub$quantity == "r_sens_disatt_range", ]
    segments(lv$lo, y + 0.18, lv$hi, y + 0.18)
    points(lv$estimate, y + 0.18, pch = 16, cex = 1.1)
    segments(sn$lo, y - 0.05, sn$hi, y - 0.05)
    points(sn$estimate, y - 0.05, pch = 21, bg = "white", cex = 1.1)
    if (nrow(br) == 1) {
      yb <- y - 0.30
      segments(br$lo, yb, br$hi, yb, col = gray_ref, lwd = 1.6)
      segments(c(br$lo, br$hi), yb - 0.06, c(br$lo, br$hi), yb + 0.06,
               col = gray_ref, lwd = 1.6)
    }
  }
  legend("bottomright", cex = 0.62, bty = "n",
         legend = c("confidence LEVEL (r_level)",
                    "sensitivity (r_sens; CD: d'-partialed / Lund: corrected)",
                    "r_sens disattenuation range (not a corrected point)"),
         pch = c(16, 21, NA), lty = c(1, 1, 1),
         col = c("black", "black", gray_ref),
         pt.bg = c(NA, "white", NA), lwd = c(1, 1, 1.6))

  par(mar = c(4.2, 0.8, 3, 0.8))
  plot(NA, xlim = xlim_r, ylim = c(-0.35, 5.6), yaxt = "n", xlab = "",
       ylab = "", main = "Gap vs simulated null reference", cex.main = 1)
  abline(v = 0, col = "gray70", lty = 3)
  mtext("r_level - r_sens_raw (Lund row: corrected gap)", side = 1,
        line = 2.4, cex = 0.8)
  for (rid in names(ys)) {
    y <- ys[[rid]]
    sub <- pdat[pdat$row_id == rid & pdat$panel == "gap", ]
    gp <- sub[grepl("^gap", sub$quantity), ]
    segments(gp$lo, y + 0.15, gp$hi, y + 0.15)
    points(gp$estimate, y + 0.15,
           pch = if (gp$style == "diamond_filled") 18 else 5, cex = 1.2)
    n30 <- sub[sub$quantity == "nullref_rho30", ]
    n45 <- sub[sub$quantity == "nullref_rho45", ]
    if (nrow(n30) == 1) {
      segments(n30$lo, y - 0.10, n30$hi, y - 0.10, col = gray_ref, lty = 3,
               lwd = 1.8)
      segments(n45$lo, y - 0.30, n45$hi, y - 0.30, col = gray_ref, lty = 2,
               lwd = 1.4)
    }
  }
  legend("bottomright", cex = 0.62, bty = "n",
         legend = c("raw gap (Lund row: corrected gap, CONFIRMATORY)",
                    "null reference, true rho = 0.30 (SIMULATION)",
                    "null reference, true rho = 0.45 (SIMULATION)"),
         pch = c(5, NA, NA), lty = c(1, 3, 2),
         col = c("black", gray_ref, gray_ref), lwd = c(1, 1.8, 1.4))
}

pdf(here("figures", "fig_confdb_contrast.pdf"), width = 9, height = 5.5)
draw_fig()
dev.off()
png(here("figures", "fig_confdb_contrast.png"), width = 9 * 300,
    height = 5.5 * 300, res = 300)
draw_fig()
dev.off()
message("Wrote figures/fig_confdb_contrast.pdf/.png and results/fig_confdb_contrast_data.csv")

cap_path <- here("figures", "fig_captions.md")
existing <- readLines(cap_path)
if (!any(grepl("fig_confdb_contrast", existing))) {
  gr_arb <- cd_get("Arbuzova_unpub_3", "gap_raw")
  cap <- c("",
    "## Figure 2 — `figures/fig_confdb_contrast.pdf/.png` (data: `results/fig_confdb_contrast_data.csv`)",
    "",
    paste0(
      "**Confidence Database level-vs-sensitivity contrast (NON-CONFIRMATORY, robustness tier) with the Lund confirmatory pair for comparison.** ",
      "Left: cross-task correlations of confidence LEVEL (filled) and metacognitive sensitivity (open; CD rows d'-partialed AUROC2, Lund row calibration-corrected) with 95% CIs; ",
      "gray brackets are the r_sens disattenuation RANGES (r / sqrt(rel_a x rel_b) at the committed reliability point estimates and both CI bounds) — sensitivity ranges, never corrected point values. ",
      "Right: the RAW gap (r_level - r_sens_raw) per CD study with 95% CI, against SIMULATION null-reference intervals — the 2.5-97.5% span of the raw gap in a true-null world (true r_level = true r_sens = 0.30 dotted / 0.45 dashed) matched to each study's n, trials, scale, and reliabilities. ",
      "**Raw gaps are inflated by differential attenuation** (level split-half reliabilities 0.97-1.00 vs AUROC2 0.62-0.90), so a positive raw gap is expected even under a true null; only ",
      sprintf("Arbuzova_unpub_3 memory x visual exceeds its null-reference intervals (raw gap %.3f vs upper bounds %.3f / %.3f). ",
              gr_arb$estimate, gr_arb$nullref_rho30_hi, gr_arb$nullref_rho45_hi),
      "This tier tests TASK-generalization with exactly one external cross-domain pair: only Arbuzova_unpub_3 memory x visual is cross-DOMAIN; Arbuzova_2021_expt2 is a visuomotor condition x condition pair, and Bang_2019_Exp2 and Haddara_2022_Expt1 are within-vision task pairs — these four studies are not four cross-domain replications. ",
      "The Lund row (top) is the pre-specified CONFIRMATORY (corrected) result: its gap (0.113, 95% CI [-0.138, 0.364]) does not exclude a null. ",
      "Null-reference intervals are a simulation reference, not a correction and not a test."))
  writeLines(c(existing, cap), cap_path)
  message("Caption appended to figures/fig_captions.md")
} else {
  message("Caption already present in figures/fig_captions.md — not duplicated")
}
