parcel_measures <- function(trials, k) {
  p <- (seq_len(nrow(trials)) - 1L) %% k + 1L
  s <- vapply(seq_len(k), function(j) {
    auroc2(trials$confidence[p == j], trials$accuracy[p == j])
  }, numeric(1))
  l <- vapply(seq_len(k), function(j) {
    mean(trials$confidence[p == j])
  }, numeric(1))
  c(
    setNames(s, paste0("s", seq_len(k))),
    setNames(l, paste0("l", seq_len(k)))
  )
}

d_prime_hat <- function(trials) {
  trials <- trials[!is.na(trials$signal) & !is.na(trials$response), ]
  n_s <- sum(trials$signal == 1)
  n_n <- sum(trials$signal == 0)
  h <- (sum(trials$response == 1 & trials$signal == 1) + 0.5) / (n_s + 1)
  f <- (sum(trials$response == 1 & trials$signal == 0) + 0.5) / (n_n + 1)
  qnorm(h) - qnorm(f)
}

latent_model <- function(k_a, k_b) {
  fac <- function(name, lab, ind) {
    c(
      paste0(name, " =~ ", paste(paste0("load_", lab, "*", ind),
                                 collapse = " + ")),
      paste0(ind, " ~~ resid_", lab, "*", ind)
    )
  }
  paste(c(
    fac("sensA", "sA", paste0("sA", seq_len(k_a))),
    fac("sensB", "sB", paste0("sB", seq_len(k_b))),
    fac("levA", "lA", paste0("lA", seq_len(k_a))),
    fac("levB", "lB", paste0("lB", seq_len(k_b))),
    "sensA ~~ rs*sensB",
    "levA ~~ rl*levB",
    "gap := rl - rs"
  ), collapse = "\n")
}

fit_latent <- function(dat, k_a, k_b) {
  families <- list(
    paste0("sA", seq_len(k_a)), paste0("sB", seq_len(k_b)),
    paste0("lA", seq_len(k_a)), paste0("lB", seq_len(k_b))
  )
  for (fam in families) {
    pooled_sd <- sqrt(mean(vapply(dat[fam], var, numeric(1), na.rm = TRUE)))
    dat[fam] <- lapply(dat[fam], function(x) {
      (x - mean(x, na.rm = TRUE)) / pooled_sd
    })
  }
  model <- latent_model(k_a, k_b)
  try_fit <- function(...) {
    tryCatch(suppressWarnings(lavaan::cfa(model, data = dat, std.lv = TRUE,
                                          ...)),
             error = function(e) NULL)
  }
  is_ok <- function(f) {
    !is.null(f) && lavaan::lavInspect(f, "converged") &&
      isTRUE(tryCatch(suppressWarnings(lavaan::lavInspect(f, "post.check")),
                      error = function(e) FALSE))
  }
  fit <- try_fit()
  if (!is_ok(fit)) {
    fit2 <- try_fit(bounds = "standard")
    if (is_ok(fit2)) fit <- fit2
  }
  bad <- data.frame(
    rho_sens = NA_real_, rho_sens_se = NA_real_,
    rho_sens_lo = NA_real_, rho_sens_hi = NA_real_,
    rho_level = NA_real_, rho_level_se = NA_real_,
    rho_level_lo = NA_real_, rho_level_hi = NA_real_,
    gap = NA_real_, gap_se = NA_real_, gap_lo = NA_real_, gap_hi = NA_real_,
    cov_rs_rl = NA_real_,
    converged = FALSE, admissible = FALSE
  )
  if (is.null(fit)) {
    return(bad)
  }
  converged <- lavaan::lavInspect(fit, "converged")
  if (!converged) {
    bad$converged <- FALSE
    return(bad)
  }
  admissible <- tryCatch(
    suppressWarnings(lavaan::lavInspect(fit, "post.check")),
    error = function(e) FALSE
  )
  pe <- lavaan::parameterEstimates(fit)
  grab <- function(lab) pe[pe$label == lab, c("est", "se", "ci.lower",
                                              "ci.upper")]
  rs <- grab("rs")
  rl <- grab("rl")
  gp <- grab("gap")
  cov_rs_rl <- tryCatch({
    v <- lavaan::lavInspect(fit, "vcov")
    if (all(c("rs", "rl") %in% rownames(v))) v["rs", "rl"] else NA_real_
  }, error = function(e) NA_real_)
  data.frame(
    rho_sens = rs$est, rho_sens_se = rs$se,
    rho_sens_lo = rs$ci.lower, rho_sens_hi = rs$ci.upper,
    rho_level = rl$est, rho_level_se = rl$se,
    rho_level_lo = rl$ci.lower, rho_level_hi = rl$ci.upper,
    gap = gp$est, gap_se = gp$se, gap_lo = gp$ci.lower, gap_hi = gp$ci.upper,
    cov_rs_rl = cov_rs_rl,
    converged = TRUE, admissible = isTRUE(admissible)
  )
}
