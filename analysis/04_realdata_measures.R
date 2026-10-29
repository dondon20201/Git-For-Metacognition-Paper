library(here)
source(here("R", "measures.R"))
source(here("R", "parcels.R"))

raw_path <- here("data", "raw", "lund2025", "metacognition_trialData_master.csv")
md5 <- unname(tools::md5sum(raw_path))
if (!identical(md5, "500b2e6e1eb0e47e717feba5b912ad9f")) {
  stop("Raw trial file md5 mismatch (", md5, ") — provenance broken, stopping.")
}
trials <- read.csv(raw_path, stringsAsFactors = FALSE)

k_domain <- c(memory = 4L, vision = 4L, Calories = 2L, GDP = 2L)
min_trials <- 60L

n_conf_na <- sum(is.na(trials$confidence))
message(sprintf("Trials with missing confidence dropped: %d of %d (%.2f%%)",
                n_conf_na, nrow(trials), 100 * n_conf_na / nrow(trials)))
trials <- trials[!is.na(trials$confidence), ]

subject_rows <- function(dat, trial_set_label) {
  out <- list()
  i <- 0
  for (dom in names(k_domain)) {
    k <- k_domain[[dom]]
    dd <- dat[dat$modality == dom, ]
    for (subj in sort(unique(dd$subject))) {
      sd_ <- dd[dd$subject == subj, ]
      sd_ <- sd_[order(sd_$trial), ]
      n <- nrow(sd_)
      qualifies <- n >= min_trials &&
        length(unique(sd_$accuracy)) == 2L
      i <- i + 1
      row <- data.frame(
        data_type = "REALDATA", trial_set = trial_set_label,
        subject = subj, domain = dom, n_trials = n, qualifies = qualifies
      )
      if (qualifies) {
        pm <- parcel_measures(sd_, k)
        row$mean_conf <- mean(sd_$confidence)
        row$mean_accuracy <- mean(sd_$accuracy)
        row$overconfidence <- (row$mean_conf - 1) / 6 - row$mean_accuracy
        row$auroc2 <- auroc2(sd_$confidence, sd_$accuracy)
        row$d_prime_hat <- d_prime_hat(sd_)
        for (j in 1:4) {
          row[[paste0("parcel_auroc2_", j)]] <-
            if (j <= k) unname(pm[paste0("s", j)]) else NA_real_
          row[[paste0("parcel_mean_conf_", j)]] <-
            if (j <= k) unname(pm[paste0("l", j)]) else NA_real_
        }
      } else {
        for (cl in c("mean_conf", "mean_accuracy", "overconfidence", "auroc2",
                     "d_prime_hat", paste0("parcel_auroc2_", 1:4),
                     paste0("parcel_mean_conf_", 1:4))) {
          row[[cl]] <- NA_real_
        }
      }
      out[[i]] <- row
    }
  }
  do.call(rbind, out)
}

measures <- rbind(
  subject_rows(trials[trials$repeat_exclude == 0, ], "primary"),
  subject_rows(trials, "all")
)

domain_summary <- do.call(rbind, lapply(split(
  measures, list(measures$trial_set, measures$domain), drop = TRUE
), function(d) {
  q <- d[d$qualifies, ]
  data.frame(
    data_type = "REALDATA", trial_set = d$trial_set[1], domain = d$domain[1],
    n_subject_domains = nrow(d),
    n_qualifying = nrow(q),
    n_excluded = nrow(d) - nrow(q),
    pct_excluded = 100 * (nrow(d) - nrow(q)) / nrow(d),
    median_trials = median(q$n_trials),
    mean_conf_mean = mean(q$mean_conf), mean_conf_sd = sd(q$mean_conf),
    accuracy_mean = mean(q$mean_accuracy), accuracy_sd = sd(q$mean_accuracy),
    overconf_mean = mean(q$overconfidence), overconf_sd = sd(q$overconfidence),
    auroc2_mean = mean(q$auroc2), auroc2_sd = sd(q$auroc2),
    d_prime_mean = mean(q$d_prime_hat), d_prime_sd = sd(q$d_prime_hat)
  )
}))
rownames(domain_summary) <- NULL

message("\n=== REALDATA domain summary ===")
print(domain_summary, digits = 3, row.names = FALSE)

write.csv(measures, here("results", "step1_REALDATA_subject_measures.csv"),
          row.names = FALSE)
write.csv(domain_summary, here("results", "step1_REALDATA_domain_summary.csv"),
          row.names = FALSE)
message("\nStep 1 complete. NO cross-domain quantity computed (firewall).")
