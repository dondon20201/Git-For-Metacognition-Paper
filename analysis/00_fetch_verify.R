library(here)

pinned_commit <- "0c38aba4eb720a143d126b7c7d4c313a900bc536"
raw_base <- paste0(
  "https://raw.githubusercontent.com/embodied-computation-group/",
  "dg-metacognition/", pinned_commit, "/data_summary/"
)

files <- list(
  trial = list(
    name = "metacognition_trialData_master.csv",
    expected_md5 = "500b2e6e1eb0e47e717feba5b912ad9f",
    expected_bytes = 11837935L
  ),
  mle = list(
    name = "mle_fit_data_wide.csv",
    expected_md5 = NA_character_,
    expected_bytes = NA_integer_
  )
)
wrongcase_md5 <- "653437023bcee998481bc33a77fd2729"

raw_dir <- here("data", "raw", "lund2025")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(here("results"), recursive = TRUE, showWarnings = FALSE)

fetch <- function(url, dest, tries = 4) {
  for (i in seq_len(tries)) {
    ok <- tryCatch({
      download.file(url, dest, mode = "wb", quiet = TRUE, method = "libcurl")
      TRUE
    }, error = function(e) {
      message(sprintf("  attempt %d failed: %s", i, conditionMessage(e)))
      FALSE
    })
    if (ok) return(invisible(TRUE))
    Sys.sleep(2^i)
  }
  stop("Download failed after ", tries, " attempts: ", url)
}

for (f in files) {
  dest <- file.path(raw_dir, f$name)
  if (file.exists(dest)) {
    message("Already present, not re-downloading: ", f$name)
  } else {
    message("Fetching ", f$name, " ...")
    fetch(paste0(raw_base, f$name), dest)
  }
}

trial_path <- file.path(raw_dir, files$trial$name)
mle_path <- file.path(raw_dir, files$mle$name)

trial_md5 <- unname(tools::md5sum(trial_path))
trial_bytes <- file.size(trial_path)
mle_md5 <- unname(tools::md5sum(mle_path))
mle_bytes <- file.size(mle_path)

message(sprintf("trial csv: md5 %s, %d bytes", trial_md5, trial_bytes))
message(sprintf("mle csv:   md5 %s, %d bytes", mle_md5, mle_bytes))

if (identical(trial_md5, wrongcase_md5)) {
  stop(
    "FATAL: fetched the wrong-case file (metacognition_TrialData_master.csv, ",
    "md5 653437...). Fix the URL casing before proceeding."
  )
}
md5_match <- identical(trial_md5, files$trial$expected_md5)
size_match <- identical(as.integer(trial_bytes), files$trial$expected_bytes)
if (!md5_match) {
  message(
    "WARNING: trial-file md5 does not match expected ",
    files$trial$expected_md5, " — structure checks below must hold, and the ",
    "discrepancy must be reported prominently."
  )
}

trials <- read.csv(trial_path, stringsAsFactors = FALSE)

expected_cols <- c(
  "subject", "modality", "trial", "signal", "response", "accuracy",
  "confidence", "rt", "rt_conf", "contrast", "repeat_exclude"
)

n_rows <- nrow(trials)
col_names <- names(trials)
n_subjects <- length(unique(trials$subject))
modalities <- sort(unique(trials$modality))
conf_values <- sort(unique(trials$confidence[!is.na(trials$confidence)]))

trials_per_subj_mod <- aggregate(
  trial ~ subject + modality, data = trials, FUN = length
)
median_trials <- aggregate(
  trial ~ modality, data = trials_per_subj_mod, FUN = median
)
names(median_trials) <- c("modality", "median_trials_per_subject")

structure_checks <- list(
  n_rows = n_rows,
  n_rows_expected = 193200L,
  columns = col_names,
  columns_match = identical(sort(col_names), sort(expected_cols)),
  n_subjects = n_subjects,
  n_subjects_expected = 322L,
  modalities = modalities,
  modalities_expected = sort(c("memory", "Calories", "GDP", "vision")),
  confidence_values = conf_values,
  confidence_expected = 1:7,
  median_trials_per_subject = median_trials
)

structure_ok <-
  n_rows == 193200L &&
  structure_checks$columns_match &&
  n_subjects == 322L &&
  identical(modalities, structure_checks$modalities_expected) &&
  identical(as.integer(conf_values), 1:7)

message("structure_ok: ", structure_ok)
print(median_trials)

if (!md5_match && !structure_ok) {
  stop("FATAL: md5 AND structure both deviate from expectations. Stopping.")
}
if (!structure_ok) {
  message("WARNING: structure deviates from expectations — see provenance JSON.")
}

per_subj <- function(dat) {
  ps <- aggregate(
    cbind(confidence, accuracy) ~ subject + modality, data = dat, FUN = mean
  )
  names(ps) <- c("subject", "modality", "mean_confidence", "mean_accuracy")
  ps
}
ps_all <- per_subj(trials)
ps_noexcl <- per_subj(trials[trials$repeat_exclude == 0, ])

domain_summary <- function(ps, label) {
  out <- aggregate(
    cbind(mean_confidence, mean_accuracy) ~ modality, data = ps, FUN = mean
  )
  out$n_subjects <- aggregate(
    subject ~ modality, data = ps, FUN = function(x) length(unique(x))
  )$subject
  out$trial_filter <- label
  out
}
own <- rbind(
  domain_summary(ps_all, "all_trials"),
  domain_summary(ps_noexcl, "repeat_exclude_dropped")
)

mle <- read.csv(mle_path, stringsAsFactors = FALSE)
avg_conf_cols <- grep("^avg_conf", names(mle), value = TRUE)
message("shipped avg_conf columns: ", paste(avg_conf_cols, collapse = ", "))
message("shipped mle file: ", nrow(mle), " rows, ", ncol(mle), " cols")

shipped <- data.frame(
  shipped_column = avg_conf_cols,
  shipped_mean_confidence = vapply(
    avg_conf_cols, function(cl) mean(mle[[cl]], na.rm = TRUE), numeric(1)
  ),
  shipped_n = vapply(
    avg_conf_cols, function(cl) sum(!is.na(mle[[cl]])), integer(1)
  )
)

suffix <- tolower(sub("^avg_conf_?", "", shipped$shipped_column))
mod_lookup <- setNames(modalities, tolower(modalities))
shipped$modality <- unname(mod_lookup[suffix])

comparison <- merge(own, shipped[, c("modality", "shipped_column",
                                     "shipped_mean_confidence", "shipped_n")],
                    by = "modality", all.x = TRUE)
comparison$abs_diff <- abs(
  comparison$mean_confidence - comparison$shipped_mean_confidence
)
comparison <- comparison[order(comparison$trial_filter, comparison$modality), ]

stopifnot("subject" %in% names(mle))
subj_check <- do.call(rbind, lapply(seq_along(avg_conf_cols), function(i) {
  cl <- avg_conf_cols[i]
  mod <- unname(mod_lookup[tolower(sub("^avg_conf_?", "", cl))])
  m <- merge(
    ps_all[ps_all$modality == mod, c("subject", "mean_confidence")],
    mle[!is.na(mle[[cl]]), c("subject", cl)],
    by = "subject"
  )
  data.frame(
    modality = mod,
    n_matched_subjects = nrow(m),
    max_abs_per_subject_diff = max(abs(m$mean_confidence - m[[cl]])),
    cor_own_vs_shipped = cor(m$mean_confidence, m[[cl]])
  )
}))
message("\n=== Per-subject agreement on matched subjects (all trials) ===")
print(subj_check, digits = 4)

message("\n=== Basic-summary cross-check (mean confidence, own vs shipped) ===")
print(comparison, digits = 4)

write.csv(comparison, here("results", "step0a_basic_summary.csv"),
          row.names = FALSE)

provenance <- list(
  fetched_utc = format(Sys.time(), "%Y-%m-%d %H:%M:%S", tz = "UTC"),
  source_repo = "https://github.com/embodied-computation-group/dg-metacognition",
  pinned_commit = pinned_commit,
  files = list(
    trial = list(
      name = files$trial$name, url = paste0(raw_base, files$trial$name),
      md5 = trial_md5, bytes = trial_bytes,
      expected_md5 = files$trial$expected_md5,
      expected_bytes = files$trial$expected_bytes,
      md5_match = md5_match, size_match = size_match
    ),
    mle = list(
      name = files$mle$name, url = paste0(raw_base, files$mle$name),
      md5 = mle_md5, bytes = mle_bytes
    )
  ),
  structure = structure_checks,
  structure_ok = structure_ok,
  max_abs_diff_mean_confidence = max(comparison$abs_diff, na.rm = TRUE),
  per_subject_agreement = subj_check
)
jsonlite::write_json(
  provenance, here("results", "step0a_provenance.json"),
  auto_unbox = TRUE, pretty = TRUE, digits = 8
)

message("\nStep 0a complete. Provenance written to results/step0a_provenance.json")
