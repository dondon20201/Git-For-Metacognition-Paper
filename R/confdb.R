read_confdb_csv <- function(path) {
  attempt <- function(...) {
    tryCatch(
      list(ok = TRUE, data = read.csv(path, stringsAsFactors = FALSE, ...)),
      error = function(e) list(ok = FALSE, err = conditionMessage(e)),
      warning = function(w) list(ok = FALSE, err = conditionMessage(w))
    )
  }
  a <- attempt()
  if (a$ok && ncol(a$data) >= 2) {
    return(list(data = a$data, status = "OK", note = "read.csv default"))
  }
  b <- attempt(fileEncoding = "UTF-8-BOM")
  if (b$ok && ncol(b$data) >= 2) {
    return(list(data = b$data, status = "OK_WITH_FIX",
                note = "required fileEncoding=UTF-8-BOM"))
  }
  c_ <- attempt(fileEncoding = "latin1")
  if (c_$ok && ncol(c_$data) >= 2) {
    return(list(data = c_$data, status = "OK_WITH_FIX",
                note = "required fileEncoding=latin1"))
  }
  reason <- if (a$ok) "parsed to a single column (delimiter not comma?)"
            else a$err
  list(data = NULL, status = "FAILED", note = reason)
}

find_col <- function(df, name) {
  hit <- which(tolower(names(df)) == tolower(name))
  if (length(hit) == 0) NA_character_ else names(df)[hit[1]]
}

find_task_col <- function(df) {
  exact <- find_col(df, "Task")
  if (!is.na(exact)) return(exact)
  hit <- grep("task", names(df), ignore.case = TRUE)
  if (length(hit) == 0) NA_character_ else names(df)[hit[1]]
}

derive_accuracy <- function(df) {
  acc_col <- find_col(df, "Accuracy")
  if (!is.na(acc_col)) {
    v <- suppressWarnings(as.numeric(df[[acc_col]]))
    vals <- unique(v[!is.na(v)])
    if (length(vals) > 0 && all(vals %in% c(0, 1))) {
      return(list(mechanism = "column", acc_col = acc_col, acc = v))
    }
  }
  s_col <- find_col(df, "Stimulus")
  r_col <- find_col(df, "Response")
  if (!is.na(s_col) && !is.na(r_col)) {
    s <- trimws(as.character(df[[s_col]]))
    r <- trimws(as.character(df[[r_col]]))
    s[s %in% c("", "NA", "NaN")] <- NA
    r[r %in% c("", "NA", "NaN")] <- NA
    acc <- as.numeric(s == r)
    got <- acc[!is.na(acc)]
    if (length(got) > 0 && any(got == 1) && any(got == 0)) {
      return(list(mechanism = "stim_resp", acc_col = NA_character_, acc = acc))
    }
  }
  list(mechanism = "none", acc_col = NA_character_,
       acc = rep(NA_real_, nrow(df)))
}

find_subj_col <- function(df) {
  sc <- find_col(df, "Subj_idx")
  if (is.na(sc)) sc <- find_col(df, "subject_id")
  sc
}

clean_trials <- function(df) {
  n0 <- nrow(df)
  train_col <- find_col(df, "Training")
  n_training <- 0L
  if (!is.na(train_col)) {
    tv <- tolower(trimws(as.character(df[[train_col]])))
    is_training <- !is.na(df[[train_col]]) & grepl("^1$|^training", tv)
    n_training <- sum(is_training)
    df <- df[!is_training, , drop = FALSE]
  }
  conf_col <- find_col(df, "Confidence")
  conf <- if (is.na(conf_col)) rep(NA_real_, nrow(df))
          else suppressWarnings(as.numeric(df[[conf_col]]))
  acc_info <- derive_accuracy(df)
  subj_col <- find_subj_col(df)
  task_col <- find_task_col(df)
  df$.conf <- conf
  df$.acc <- acc_info$acc
  df$.subj <- if (is.na(subj_col)) rep(NA_character_, nrow(df))
              else as.character(df[[subj_col]])
  df$.task <- if (is.na(task_col)) rep("single_task", nrow(df))
              else as.character(df[[task_col]])
  keep <- !is.na(df$.conf) & !is.na(df$.acc) & !is.na(df$.subj) &
    !is.na(df$.task)
  list(
    data = df[keep, , drop = FALSE],
    mechanism = acc_info$mechanism,
    acc_col = acc_info$acc_col,
    task_col = task_col,
    subj_col = subj_col,
    conf_col = conf_col,
    train_col = train_col,
    n_rows_in = n0,
    n_training_dropped = n_training,
    n_clean = sum(keep)
  )
}
