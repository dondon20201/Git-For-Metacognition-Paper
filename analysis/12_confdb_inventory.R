library(here)
source(here("R", "confdb.R"))

raw_dir <- here("data", "raw", "confidence_db")
manifest <- read.csv(here("results", "step11_CONFDB_manifest.csv"),
                     stringsAsFactors = FALSE)
data_files <- manifest[manifest$role == "data", ]
message("Manifest: ", nrow(data_files), " study data files, ",
        sum(data_files$status == "VERIFIED"), " VERIFIED")

std_cols <- c("Subj_idx", "Stimulus", "Response", "Confidence",
              "RT_dec", "RT_conf", "Task", "Training", "Accuracy")

inv_rows <- list()
sel_rows <- list()

for (i in seq_len(nrow(data_files))) {
  fn <- data_files$filename[i]
  study <- sub("^data_", "", sub("\\.csv$", "", fn))
  path <- file.path(raw_dir, fn)

  inv <- data.frame(filename = fn, study = study,
                    manifest_status = data_files$status[i],
                    stringsAsFactors = FALSE)
  sel <- data.frame(study = study, stringsAsFactors = FALSE)

  if (data_files$status[i] != "VERIFIED" || !file.exists(path)) {
    inv$parse_status <- "NOT_DOWNLOADED"
    inv$parse_note <- paste("manifest status:", data_files$status[i])
    inv_rows[[i]] <- inv
    sel$s1_pass <- FALSE; sel$s2_pass <- NA; sel$s3_pass <- NA
    sel$s4_pass <- NA; sel$selected <- FALSE
    sel$note <- "file not available (see manifest)"
    sel_rows[[i]] <- sel
    next
  }

  p <- read_confdb_csv(path)
  inv$parse_status <- p$status
  inv$parse_note <- p$note
  if (p$status == "FAILED") {
    inv_rows[[i]] <- inv
    sel$s1_pass <- FALSE; sel$s2_pass <- NA; sel$s3_pass <- NA
    sel$s4_pass <- NA; sel$selected <- FALSE
    sel$note <- "parse FAILED (see inventory)"
    sel_rows[[i]] <- sel
    next
  }
  df <- p$data

  inv$n_rows <- nrow(df)
  inv$n_cols <- ncol(df)
  inv$columns <- paste(names(df), collapse = ";")
  for (sc in std_cols) {
    inv[[paste0("has_", sc)]] <- !is.na(find_col(df, sc))
  }

  subj_col <- find_subj_col(df)
  inv$n_subjects <- if (is.na(subj_col)) NA_integer_
                    else length(unique(df[[subj_col]]))

  task_col <- find_task_col(df)
  inv$task_col <- ifelse(is.na(task_col), "", task_col)
  if (!is.na(task_col) && !is.na(subj_col)) {
    tv <- as.character(df[[task_col]])
    inv$n_task_values <- length(unique(tv[!is.na(tv)]))
    per_subj_tasks <- tapply(tv, as.character(df[[subj_col]]),
                             function(x) length(unique(x[!is.na(x)])))
    inv$max_tasks_within_subject <- max(per_subj_tasks)
    inv$n_subjects_ge2_tasks <- sum(per_subj_tasks >= 2)
  } else {
    inv$n_task_values <- if (is.na(task_col)) 0L else NA_integer_
    inv$max_tasks_within_subject <- NA_integer_
    inv$n_subjects_ge2_tasks <- 0L
  }
  inv$multi_task_within_subject <-
    isTRUE(inv$max_tasks_within_subject >= 2)

  if (!is.na(subj_col)) {
    cell <- if (!is.na(task_col)) {
      paste(as.character(df[[subj_col]]), as.character(df[[task_col]]),
            sep = "\r")
    } else {
      as.character(df[[subj_col]])
    }
    cnt <- as.numeric(table(cell))
    inv$trials_per_cell_min <- min(cnt)
    inv$trials_per_cell_median <- median(cnt)
    inv$trials_per_cell_max <- max(cnt)
  } else {
    inv$trials_per_cell_min <- NA_real_
    inv$trials_per_cell_median <- NA_real_
    inv$trials_per_cell_max <- NA_real_
  }

  conf_col <- find_col(df, "Confidence")
  if (!is.na(conf_col)) {
    raw_conf <- df[[conf_col]]
    n_missing_raw <- sum(is.na(raw_conf) |
                           trimws(as.character(raw_conf)) %in%
                             c("", "NA", "NaN"))
    num_conf <- suppressWarnings(as.numeric(raw_conf))
    nn <- num_conf[!is.na(num_conf)]
    inv$conf_numeric <- length(nn) > 0 &&
      sum(is.na(num_conf)) == n_missing_raw
    inv$conf_n_distinct <- length(unique(nn))
    inv$conf_min <- if (length(nn)) min(nn) else NA_real_
    inv$conf_max <- if (length(nn)) max(nn) else NA_real_
    inv$conf_integer <- if (length(nn)) all(abs(nn - round(nn)) < 1e-9)
                        else NA
    inv$pct_missing_conf <- round(100 * n_missing_raw / nrow(df), 2)
  } else {
    inv$conf_numeric <- FALSE
    inv$conf_n_distinct <- 0L
    inv$conf_min <- NA_real_
    inv$conf_max <- NA_real_
    inv$conf_integer <- NA
    inv$pct_missing_conf <- NA_real_
  }

  acc_info <- derive_accuracy(df)
  inv$acc_mechanism <- acc_info$mechanism
  inv$acc_col <- ifelse(is.na(acc_info$acc_col), "", acc_info$acc_col)
  train_col <- find_col(df, "Training")
  inv$training_col_present <- !is.na(train_col)
  inv$training_values <- if (!is.na(train_col)) {
    paste(sort(unique(as.character(df[[train_col]]))), collapse = ";")
  } else ""

  ct <- clean_trials(df)
  inv$n_clean_trials <- ct$n_clean
  inv$n_training_dropped <- ct$n_training_dropped
  cd <- ct$data

  sel$task_col <- inv$task_col
  sel$acc_mechanism <- ct$mechanism
  sel$n_subjects_total <- inv$n_subjects
  sel$n_tasks_observed <- inv$n_task_values

  if (ct$n_clean > 0 && !is.na(task_col)) {
    counts <- aggregate(rep(1L, nrow(cd)),
                        by = list(subj = cd$.subj, task = cd$.task), FUN = sum)
    names(counts)[3] <- "n"
    q <- counts[counts$n >= 60, ]
    per_subj_qtasks <- tapply(q$task, q$subj, function(x) length(unique(x)))
    qualified_subj <- names(per_subj_qtasks)[per_subj_qtasks >= 2]
    sel$s1_n_qualified_subjects <- length(qualified_subj)

    sel$best_pair <- ""
    sel$best_pair_n_subjects <- 0L
    if (length(qualified_subj) > 0) {
      qq <- q[q$subj %in% qualified_subj, ]
      tasks <- sort(unique(qq$task))
      if (length(tasks) >= 2) {
        prs <- t(combn(tasks, 2))
        pair_n <- apply(prs, 1, function(pr) {
          length(intersect(qq$subj[qq$task == pr[1]],
                           qq$subj[qq$task == pr[2]]))
        })
        b <- which.max(pair_n)
        sel$best_pair <- paste(prs[b, ], collapse = " x ")
        sel$best_pair_n_subjects <- pair_n[b]
      }
    }

    sel$s1_pass <- sel$s1_n_qualified_subjects >= 40 &&
      isTRUE(inv$n_task_values >= 2)

    qual_tasks <- sort(unique(q$task[q$subj %in% qualified_subj]))
    sel$qualifying_tasks <- paste(qual_tasks, collapse = ";")
    if (length(qual_tasks) > 0) {
      qc <- cd[cd$.subj %in% qualified_subj & cd$.task %in% qual_tasks, ]
      conf_distinct <- tapply(qc$.conf, qc$.task,
                              function(x) length(unique(x)))
      sel$s2_min_conf_distinct <- min(conf_distinct[qual_tasks])
      sel$s2_pass <- all(conf_distinct[qual_tasks] >= 3)
      acc_classes <- tapply(qc$.acc, qc$.task,
                            function(x) length(unique(x)))
      sel$s3_all_tasks_both_classes <- all(acc_classes[qual_tasks] == 2)
      sel$s3_pass <- ct$mechanism != "none" &&
        isTRUE(sel$s3_all_tasks_both_classes)
    } else {
      sel$s2_min_conf_distinct <- NA_real_
      sel$s2_pass <- NA
      sel$s3_all_tasks_both_classes <- NA
      sel$s3_pass <- NA
    }
  } else {
    sel$s1_n_qualified_subjects <- 0L
    sel$best_pair <- ""
    sel$best_pair_n_subjects <- 0L
    sel$s1_pass <- FALSE
    sel$qualifying_tasks <- ""
    sel$s2_min_conf_distinct <- NA_real_
    sel$s2_pass <- NA
    sel$s3_all_tasks_both_classes <- NA
    sel$s3_pass <- NA
  }

  inv_rows[[i]] <- inv
  sel_rows[[i]] <- sel
  if (i %% 30 == 0) message("  inventoried ", i, "/", nrow(data_files))
}

fill_union <- function(rows) {
  all_names <- unique(unlist(lapply(rows, names)))
  do.call(rbind, lapply(rows, function(r) {
    for (nm in setdiff(all_names, names(r))) r[[nm]] <- NA
    r[all_names]
  }))
}
inventory <- fill_union(inv_rows)
selection <- fill_union(sel_rows)

meta_path <- file.path(raw_dir, "Database_Information.xlsx")
stopifnot(file.exists(meta_path))
meta <- as.data.frame(readxl::read_excel(meta_path))
message("Metadata sheet: ", nrow(meta), " rows, columns: ",
        paste(names(meta), collapse = ", "))
flag_cols <- grep("flag|exclu|problem|warn|caution|issue",
                  names(meta), ignore.case = TRUE, value = TRUE)
name_col <- grep("^name", names(meta), ignore.case = TRUE, value = TRUE)[1]
notes_col <- grep("^notes$", names(meta), ignore.case = TRUE, value = TRUE)
problem_rx <- "exclud|problem|warning|caution|do not use|corrupt|erroneous|unreliable"
meta_key <- trimws(as.character(meta[[name_col]]))
message("S4: dedicated flag columns: ",
        if (length(flag_cols)) paste(flag_cols, collapse = ", ") else "NONE",
        "; Notes column scanned with keyword regex: ", problem_rx)
selection$s4_flag <- vapply(selection$study, function(st) {
  hit <- which(meta_key == st)
  if (length(hit) == 0) return("")
  out <- character(0)
  for (fc in flag_cols) {
    v <- meta[hit[1], fc]
    if (!is.na(v) && trimws(as.character(v)) != "") {
      out <- c(out, paste0(fc, ": ", v))
    }
  }
  for (nc in notes_col) {
    v <- as.character(meta[hit[1], nc])
    if (!is.na(v) && grepl(problem_rx, v, ignore.case = TRUE)) {
      out <- c(out, paste0("Notes: ", v))
    }
  }
  paste(out, collapse = " | ")
}, character(1))
selection$s4_pass <- selection$s4_flag == ""
selection$in_metadata_sheet <- if (!is.na(name_col)) {
  selection$study %in% trimws(as.character(meta[[name_col]]))
} else NA

isTRUE_v <- function(x) !is.na(x) & x
selection$selected <- isTRUE_v(selection$s1_pass) &
  isTRUE_v(selection$s2_pass) & isTRUE_v(selection$s3_pass) &
  isTRUE_v(selection$s4_pass)

write.csv(inventory, here("results", "step11_CONFDB_inventory.csv"),
          row.names = FALSE)
write.csv(selection, here("results", "step11_CONFDB_selection.csv"),
          row.names = FALSE)

message("\n=== Selection funnel (mechanical, pre-specified) ===")
message("data files in manifest:        ", nrow(data_files))
message("parsed OK:                     ",
        sum(inventory$parse_status %in% c("OK", "OK_WITH_FIX")))
message("multi-task within subject:     ",
        sum(inventory$multi_task_within_subject, na.rm = TRUE))
message("S1 (>=40 subj x >=60 x 2 tasks): ",
        sum(isTRUE_v(selection$s1_pass)))
message("S1+S2:                         ",
        sum(isTRUE_v(selection$s1_pass) & isTRUE_v(selection$s2_pass)))
message("S1+S2+S3:                      ",
        sum(isTRUE_v(selection$s1_pass) & isTRUE_v(selection$s2_pass) &
              isTRUE_v(selection$s3_pass)))
message("SELECTED (all of S1-S4):       ", sum(selection$selected))
message("\nSelected studies:")
print(selection$study[selection$selected])
