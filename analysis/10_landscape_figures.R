library(here)

stopifnot(dir.exists(here("figures")))

step6 <- read.csv(here("results", "step6_REALDATA_confirmatory.csv"))
step7 <- read.csv(here("results", "step7_REALDATA_descriptive.csv"))
step2 <- read.csv(here("results", "step2_REALDATA_reliability.csv"))

r3 <- function(x) round(x, 3)
assert_eq <- function(what, got, want) {
  if (!isTRUE(all.equal(r3(got), want, tolerance = 1e-8))) {
    stop(sprintf("ASSERTION FAILED for %s: read %.3f, expected %.3f — committed file does not match the audited value; stopping.",
                 what, got, want))
  }
}

c_lev <- step6[step6$block == "PRIMARY_CONFIRMATORY" &
                 step6$quantity == "rho_level_latent", ]
c_gap <- step6[step6$block == "PRIMARY_CONFIRMATORY" &
                 step6$quantity == "gap_corrected", ]
assert_eq("step6 rho_level", c_lev$estimate, 0.422)
assert_eq("step6 gap", c_gap$estimate, 0.113)
assert_eq("step6 gap ci_lo", c_gap$ci_lo, -0.138)
assert_eq("step6 gap ci_hi", c_gap$ci_hi, 0.364)

want <- data.frame(
  pair = c("memory_x_Calories", "memory_x_GDP", "vision_x_Calories",
           "vision_x_GDP", "Calories_x_GDP"),
  rho_level = c(0.375, 0.323, 0.450, 0.464, 0.473),
  gap = c(-0.020, -0.133, 0.194, 0.296, 0.237),
  band = c(0.09, 0.09, 0.09, 0.09, 0.10)
)
step7 <- step7[match(want$pair, step7$pair), ]
for (i in seq_len(nrow(want))) {
  assert_eq(paste("step7 rho_level", want$pair[i]),
            step7$rho_level[i], want$rho_level[i])
  assert_eq(paste("step7 gap", want$pair[i]),
            step7$gap_uncorrected[i], want$gap[i])
  stopifnot(step7$bias_band[i] == want$band[i])
}

rel <- step2[step2$measure == "auroc2" & step2$kind == "parcel_reliability", ]
rel_want <- c(memory = 0.546, vision = 0.273, Calories = 0.306, GDP = 0.242)
for (d in names(rel_want)) {
  assert_eq(paste("step2 parcel reliability", d),
            rel$reliability[rel$domain == d], rel_want[[d]])
}
rel_of <- function(d) rel$reliability[rel$domain == d]
message("All assertions passed — committed files match audited values.")

fig <- data.frame(
  pair = c("memory_x_vision", want$pair),
  label = c("mem × vis", "mem × Cal", "mem × GDP",
            "vis × Cal", "vis × GDP", "Cal × GDP"),
  domain_a = c("memory", "memory", "memory", "vision", "vision", "Calories"),
  domain_b = c("vision", "Calories", "GDP", "Calories", "GDP", "GDP"),
  tier = c("CONFIRMATORY (corrected)",
           rep("NON-CONFIRMATORY (uncorrected, descriptive)", 5)),
  rho_level = c(c_lev$estimate, step7$rho_level),
  rho_level_lo = c(c_lev$ci_lo, step7$rho_level_lo),
  rho_level_hi = c(c_lev$ci_hi, step7$rho_level_hi),
  gap = c(c_gap$estimate, step7$gap_uncorrected),
  gap_lo = c(c_gap$ci_lo, step7$gap_lo),
  gap_hi = c(c_gap$ci_hi, step7$gap_hi),
  band = c(NA, step7$bias_band)
)
fig$x <- seq_len(nrow(fig))
write.csv(fig, here("results", "fig_landscape_data.csv"), row.names = FALSE)

draw_landscape <- function() {
  op <- par(mfrow = c(2, 1), mar = c(1.5, 4.5, 2.5, 1), oma = c(6, 0, 0, 0),
            mgp = c(2.8, 0.8, 0))
  on.exit(par(op))
  conf <- fig$tier == "CONFIRMATORY (corrected)"
  pch <- ifelse(conf, 16, 21)
  plot(NA, xlim = c(0.5, 6.5), ylim = c(0, 0.65), xaxt = "n",
       xlab = "", ylab = expression(rho["level"]~"(latent, 95% CI)"),
       main = "Cross-domain correlation landscape")
  arrows(fig$x, fig$rho_level_lo, fig$x, fig$rho_level_hi,
         angle = 90, code = 3, length = 0.05)
  points(fig$x, fig$rho_level, pch = pch, bg = "white", cex = 1.3)
  axis(1, at = fig$x, labels = FALSE)
  ylim <- range(c(fig$gap_lo, fig$gap_hi,
                  fig$gap + fig$band, fig$gap - fig$band), na.rm = TRUE)
  ylim <- ylim + c(-0.05, 0.05)
  plot(NA, xlim = c(0.5, 6.5), ylim = ylim, xaxt = "n", xlab = "",
       ylab = expression("gap ="~rho["level"] - rho["sens"]~"(95% CI ± band)"))
  hasband <- !is.na(fig$band)
  rect(fig$x[hasband] - 0.28, (fig$gap - fig$band)[hasband],
       fig$x[hasband] + 0.28, (fig$gap + fig$band)[hasband],
       col = grDevices::adjustcolor("steelblue", alpha.f = 0.25), border = NA)
  abline(h = 0, lty = 2, col = "gray40")
  arrows(fig$x, fig$gap_lo, fig$x, fig$gap_hi,
         angle = 90, code = 3, length = 0.05)
  points(fig$x, fig$gap, pch = pch, bg = "white", cex = 1.3)
  axis(1, at = fig$x, labels = fig$label, las = 2, cex.axis = 0.95)
  legend("topleft", bty = "n", cex = 0.85,
         pch = c(16, 21, 22),
         pt.bg = c(NA, "white", grDevices::adjustcolor("steelblue", 0.25)),
         pt.cex = c(1.2, 1.2, 2),
         col = c("black", "black", "gray70"),
         legend = c("CONFIRMATORY (corrected)",
                    "NON-CONFIRMATORY (uncorrected, descriptive)",
                    "frozen systematic band (±0.09 / ±0.10)"))
}
pdf(here("figures", "fig_landscape.pdf"), width = 7, height = 8)
draw_landscape()
dev.off()
png(here("figures", "fig_landscape.png"), width = 7 * 300, height = 8 * 300,
    res = 300)
draw_landscape()
dev.off()

figs <- fig[, c("pair", "label", "domain_a", "domain_b", "tier", "gap")]
figs$rel_a <- vapply(figs$domain_a, rel_of, numeric(1))
figs$rel_b <- vapply(figs$domain_b, rel_of, numeric(1))
figs$geomean_rel <- sqrt(figs$rel_a * figs$rel_b)
write.csv(figs, here("results", "figS_gap_vs_reliability_data.csv"),
          row.names = FALSE)

draw_scatter <- function() {
  op <- par(mar = c(4.5, 4.5, 2.5, 1))
  on.exit(par(op))
  conf <- figs$tier == "CONFIRMATORY (corrected)"
  plot(figs$geomean_rel, figs$gap,
       pch = ifelse(conf, 16, 21), bg = "white", cex = 1.4,
       xlim = range(figs$geomean_rel) + c(-0.03, 0.05),
       ylim = range(figs$gap) + c(-0.08, 0.08),
       xlab = "geometric mean of pair's raw AUROC2 parcel reliabilities",
       ylab = expression("gap ="~rho["level"] - rho["sens"]),
       main = "EXPLORATORY: gap vs sensitivity parcel reliability")
  abline(h = 0, lty = 2, col = "gray40")
  text(figs$geomean_rel, figs$gap, figs$label, pos = 4, offset = 0.4,
       cex = 0.85)
  legend("bottomleft", bty = "n", cex = 0.85, pch = c(16, 21),
         pt.bg = "white",
         legend = c("CONFIRMATORY (corrected)",
                    "NON-CONFIRMATORY (uncorrected)"))
}
pdf(here("figures", "figS_gap_vs_reliability.pdf"), width = 6.5, height = 5.5)
draw_scatter()
dev.off()
png(here("figures", "figS_gap_vs_reliability.png"), width = 6.5 * 300,
    height = 5.5 * 300, res = 300)
draw_scatter()
dev.off()

writeLines(c(
  "# Figure captions",
  "",
  "## Figure 1 — `figures/fig_landscape.pdf/.png` (data: `results/fig_landscape_data.csv`)",
  "",
  "**Cross-domain correlation landscape.** Top: latent correlation of confidence",
  "level (set-point) for each of the six domain pairs, with 95% delta-method CIs.",
  "Bottom: the gap (ρ_level − ρ_sens) per pair. memory × vision (filled point) is",
  "the pre-specified CONFIRMATORY pair: its ρ_sens is",
  "calibration-corrected and its gap carries no systematic band. The five open",
  "points are NON-CONFIRMATORY / DESCRIPTIVE: their sensitivity",
  "correlations are uncorrected (these trial-count configurations failed the",
  "blind calibration gate), and each gap therefore carries a shaded rectangle —",
  "the frozen systematic-bias allowance (±0.09 for 200/100-trial pairs, ±0.10 for",
  "Calories × GDP). **The shaded bands are systematic-bias allowances, NOT",
  "confidence intervals**; whiskers are sampling CIs only. Note vision × GDP: its",
  "CI excludes 0, but its lower bound (0.055) lies inside its ±0.09 band, so no",
  "nonzero-gap claim is made for any descriptive pair. The confirmatory result",
  "(gap = 0.113, 95% CI [−0.138, 0.364]) does not exclude a null gap.",
  "",
  "## Figure S1 — `figures/figS_gap_vs_reliability.pdf/.png` (data: `results/figS_gap_vs_reliability_data.csv`)",
  "",
  "**EXPLORATORY: gap vs sensitivity measurement quality.** Each pair's gap",
  "(as in Figure 1) against the geometric mean of its two domains' raw,",
  "PRE-PARTIALING AUROC2 parcel reliabilities (memory 0.546, vision 0.273,",
  "Calories 0.306, GDP 0.242; `results/step2_REALDATA_reliability.csv`). No line",
  "is fitted, no correlation computed, no test performed: with six non-independent",
  "points this display is EXPLORATORY. The measurement-lens reading — the",
  "level-advantage appears where sensitivity is least reliably measured and",
  "reverses where it is best measured (memory pairs) — is offered as",
  "interpretation of residual attenuation, not as a tested result."
), here("figures", "fig_captions.md"))
message("Figures + captions + underlying-data CSVs written.")
