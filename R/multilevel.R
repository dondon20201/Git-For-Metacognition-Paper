library(lme4)

step4_control <- lme4::lmerControl(optimizer = "bobyqa",
                                   optCtrl = list(maxfun = 2e5))

step4_formula <- function(rung, drop_slope_at = NULL) {
  base <- "conf ~ 0 + domain + domain:acc_dc + "
  re <- switch(rung,
    M1 = "(1 + acc_dc | subject) + (1 + acc_dc | subject_domain)",
    M2 = "(1 + acc_dc || subject) + (1 + acc_dc || subject_domain)",
    M3 = if (identical(drop_slope_at, "subject")) {
      "(1 | subject) + (1 + acc_dc || subject_domain)"
    } else {
      "(1 + acc_dc || subject) + (1 | subject_domain)"
    }
  )
  as.formula(paste0(base, re))
}

step4_varcomps <- function(fit) {
  vc <- as.data.frame(lme4::VarCorr(fit))
  vc <- vc[is.na(vc$var2), c("grp", "var1", "vcov")]
  vc$grp <- sub("\\.[0-9]+$", "", vc$grp)
  vc$var1[is.na(vc$var1)] <- "Residual"
  stats::aggregate(vcov ~ grp + var1, data = vc, FUN = sum)
}

vc_get <- function(vc, grp, term) {
  v <- vc$vcov[vc$grp == grp & vc$var1 == term]
  if (length(v) == 1) v else NA_real_
}

step4_fit_one <- function(dat, rung, drop_slope_at = NULL) {
  f <- step4_formula(rung, drop_slope_at)
  fit <- tryCatch(
    suppressMessages(suppressWarnings(
      lme4::lmer(f, data = dat, REML = TRUE, control = step4_control)
    )),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(list(fit = NULL, converged = FALSE, singular = NA))
  }
  msgs <- fit@optinfo$conv$lme4$messages
  converged <- fit@optinfo$conv$opt == 0 &&
    (is.null(msgs) || length(msgs) == 0)
  singular <- lme4::isSingular(fit, tol = 1e-4)
  list(fit = fit, converged = converged, singular = singular)
}

step4_ladder <- function(dat, start_rung = "M1") {
  rungs <- c("M1", "M2", "M3")
  rungs <- rungs[match(start_rung, rungs):3]
  last <- NULL
  for (rung in rungs) {
    if (rung != "M3") {
      r <- step4_fit_one(dat, rung)
      if (r$converged && !isTRUE(r$singular)) {
        return(step4_result(r, rung, NULL))
      }
      last <- r
    } else {
      drop_at <- "subject_domain"
      if (!is.null(last$fit)) {
        vc <- step4_varcomps(last$fit)
        vs <- vc_get(vc, "subject", "acc_dc")
        vsd <- vc_get(vc, "subject_domain", "acc_dc")
        if (!is.na(vs) && !is.na(vsd)) {
          drop_at <- if (vs <= vsd) "subject" else "subject_domain"
        }
      }
      r <- step4_fit_one(dat, "M3", drop_at)
      return(step4_result(r, "M3", drop_at))
    }
  }
}

step4_result <- function(r, rung, drop_at) {
  ok <- !is.null(r$fit)
  vc <- if (ok) step4_varcomps(r$fit) else NULL
  v_int_s <- if (ok) vc_get(vc, "subject", "(Intercept)") else NA_real_
  v_int_sd <- if (ok) vc_get(vc, "subject_domain", "(Intercept)") else NA_real_
  v_slp_s <- if (ok) vc_get(vc, "subject", "acc_dc") else NA_real_
  v_slp_sd <- if (ok) vc_get(vc, "subject_domain", "acc_dc") else NA_real_
  v_res <- if (ok) vc_get(vc, "Residual", "Residual") else NA_real_
  icc_level <- v_int_s / (v_int_s + v_int_sd)
  icc_sens <- if (rung == "M3") NA_real_ else v_slp_s / (v_slp_s + v_slp_sd)
  list(
    fit = r$fit, rung = rung, drop_slope_at = drop_at,
    converged = isTRUE(r$converged), singular = isTRUE(r$singular),
    v_int_subject = v_int_s, v_int_subject_domain = v_int_sd,
    v_slope_subject = v_slp_s, v_slope_subject_domain = v_slp_sd,
    v_residual = v_res,
    icc_level = icc_level, icc_sens = icc_sens,
    delta_icc = icc_level - icc_sens
  )
}

add_acc_dc <- function(dat) {
  dm <- tapply(dat$accuracy, dat$domain, mean)
  dat$acc_dc <- dat$accuracy - dm[as.character(dat$domain)]
  dat
}

ipsative_exhibit <- function(dat) {
  grand <- tapply(dat$conf, dat$subject, mean)
  dommean <- stats::aggregate(conf ~ subject + domain, data = dat, FUN = mean)
  dommean$centered <- dommean$conf - grand[as.character(dommean$subject)]
  wide <- stats::reshape(dommean[, c("subject", "domain", "centered")],
                         idvar = "subject", timevar = "domain",
                         direction = "wide")
  doms <- sub("^centered\\.", "", grep("^centered\\.", names(wide),
                                       value = TRUE))
  prs <- t(combn(doms, 2))
  data.frame(
    domain_a = prs[, 1], domain_b = prs[, 2],
    r_centered = apply(prs, 1, function(p) {
      cor(wide[[paste0("centered.", p[1])]],
          wide[[paste0("centered.", p[2])]], use = "complete.obs")
    })
  )
}
