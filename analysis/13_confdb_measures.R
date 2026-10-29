library(here)
source(here("R", "confdb.R"))
source(here("R", "measures.R"))

set.seed(20260710)
n_boot <- 500L
min_trials <- 60L

raw_dir <- here("data", "raw", "confidence_db")
selection <- read.csv(here("results", "step11_CONFDB_selection.csv"),
                      stringsAsFactors = FALSE)
selected <- selection$study[selection$selected]
message("Selected studies: ", length(selected))
stopifnot(length(selected) > 0)

dprime_confdb <- function(stim, resp, signal_value) {
  is_sig <- stim == signal_value
  n_s <- sum(is_sig)
  n_n <- sum(!is_sig)
  if (n_s == 0 || n_n == 0) return(NA_real_)
  h <- (sum(resp == signal_value & is_sig) + 0.5) / (n_s + 1)
  f <- (sum(resp == signal_value & !is_sig) + 0.5) / (n_n + 1)
  qnorm(h) - qnorm(f)
}

sb2 <- function(r) 2 * r / (1 + r)

measures_rows <- list()
rel_rows <- list()
mi <- 0L
ri <- 0L

for (study in selected) {
  fn <- paste0("data_", study, ".csv")
  message("Measuring ", study, " ...")
  p <- read_confdb_csv(file.path(raw_dir, fn))
  stopifnot(p$status %in% c("OK", "OK_WITH_FIX"))
  ct <- clean_trials(p$data)
  cd <- ct$data

  conf_min <- min(cd$.conf)
  conf_max <- max(cd$.conf)
  stopifnot(conf_max > conf_min)

  s_col <- find_col(cd, "Stimulus")
  r_col <- find_col(cd, "Response")

  for (task in sort(unique(cd$.task))) {
    td <- cd[cd$.task == task, , drop = FALSE]

    stim_chr <- if (!is.na(s_col)) trimws(as.character(td[[s_col]])) else NULL
    stim_vals <- if (!is.null(stim_chr)) {
      sort(unique(stim_chr[!is.na(stim_chr) & stim_chr != ""]))
    } else character(0)
    dprime_ok <- length(stim_vals) == 2
    dprime_note <- if (dprime_ok) "" else
      sprintf("d'-hat undefined: %d distinct stimulus values",
              length(stim_vals))

    cnt <- table(td$.subj)
    subjects <- names(cnt)[cnt >= min_trials]
    if (length(subjects) == 0) next

    half_lv_o <- half_lv_e <- half_au_o <- half_au_e <-
      rep(NA_real_, length(subjects))

    for (k in seq_along(subjects)) {
      sd_ <- td[td$.subj == subjects[k], , drop = FALSE]
      lev_raw <- mean(sd_$.conf)
      dh <- NA_real_
      if (dprime_ok) {
        rs <- trimws(as.character(sd_[[r_col]]))
        dh <- dprime_confdb(trimws(as.character(sd_[[s_col]])), rs,
                            stim_vals[2])
      }
      mi <- mi + 1L
      measures_rows[[mi]] <- data.frame(
        data_type = "REALDATA-EXTERNAL (CONFDB)",
        study = study, task = task, subj_idx = subjects[k],
        n_trials_clean = nrow(sd_),
        mean_accuracy = mean(sd_$.acc),
        mean_conf_raw = lev_raw,
        conf_scale_min = conf_min, conf_scale_max = conf_max,
        level = (lev_raw - conf_min) / (conf_max - conf_min),
        auroc2 = auroc2(sd_$.conf, sd_$.acc),
        dprime_hat = dh,
        dprime_note = dprime_note,
        stringsAsFactors = FALSE
      )
      odd <- sd_[c(TRUE, FALSE), , drop = FALSE]
      even <- sd_[c(FALSE, TRUE), , drop = FALSE]
      half_lv_o[k] <- (mean(odd$.conf) - conf_min) / (conf_max - conf_min)
      half_lv_e[k] <- (mean(even$.conf) - conf_min) / (conf_max - conf_min)
      half_au_o[k] <- auroc2(odd$.conf, odd$.acc)
      half_au_e[k] <- auroc2(even$.conf, even$.acc)
    }

    for (m in c("level", "auroc2")) {
      x <- if (m == "level") half_lv_o else half_au_o
      y <- if (m == "level") half_lv_e else half_au_e
      ok <- !is.na(x) & !is.na(y)
      xs <- x[ok]; ys <- y[ok]
      n <- length(xs)
      if (n < 3 || sd(xs) == 0 || sd(ys) == 0) {
        rel <- ci <- c(NA_real_, NA_real_, NA_real_)
      } else {
        rel <- sb2(cor(xs, ys))
        boots <- replicate(n_boot, {
          idx <- sample.int(n, replace = TRUE)
          if (sd(xs[idx]) == 0 || sd(ys[idx]) == 0) NA_real_
          else sb2(cor(xs[idx], ys[idx]))
        })
        ci <- quantile(boots, c(0.025, 0.975), names = FALSE, na.rm = TRUE)
        rel <- c(rel, ci)
      }
      ri <- ri + 1L
      rel_rows[[ri]] <- data.frame(
        data_type = "REALDATA-EXTERNAL (CONFDB)",
        study = study, task = task, measure = m,
        n_subjects = n, split = "odd_even_row_order",
        rel_sb = rel[1], ci_lo = rel[2], ci_hi = rel[3],
        n_boot = n_boot,
        stringsAsFactors = FALSE
      )
    }
  }
}

measures <- do.call(rbind, measures_rows)
reliability <- do.call(rbind, rel_rows)

stopifnot(!anyDuplicated(measures[, c("study", "task", "subj_idx")]))

write.csv(measures, here("results", "step11_CONFDB_measures.csv"),
          row.names = FALSE)
write.csv(reliability, here("results", "step11_CONFDB_reliability.csv"),
          row.names = FALSE)

message("\n=== Summary ===")
message("subject-task rows: ", nrow(measures),
        " across ", length(unique(measures$study)), " studies / ",
        nrow(unique(measures[, c("study", "task")])), " study-task cells")
message("\nSplit-half reliabilities (Spearman-Brown, 95% subject-bootstrap CI):")
print(reliability[, c("study", "task", "measure", "n_subjects",
                      "rel_sb", "ci_lo", "ci_hi")], digits = 2, row.names = FALSE)
