auroc2 <- function(confidence, accuracy) {
  ok <- !is.na(confidence) & !is.na(accuracy)
  conf <- confidence[ok]
  acc <- accuracy[ok]
  n1 <- as.numeric(sum(acc == 1))
  n0 <- as.numeric(sum(acc == 0))
  if (n1 == 0 || n0 == 0) {
    return(NA_real_)
  }
  r <- rank(conf)
  (sum(r[acc == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

gk_gamma <- function(confidence, accuracy) {
  ok <- !is.na(confidence) & !is.na(accuracy)
  conf <- confidence[ok]
  acc <- accuracy[ok]
  if (sum(acc == 1) == 0L || sum(acc == 0) == 0L) {
    return(NA_real_)
  }
  lv <- sort(unique(conf))
  n0 <- as.numeric(tabulate(match(conf[acc == 0], lv), nbins = length(lv)))
  n1 <- as.numeric(tabulate(match(conf[acc == 1], lv), nbins = length(lv)))
  n1_above <- rev(cumsum(rev(n1))) - n1
  n1_below <- cumsum(n1) - n1
  concordant <- sum(n0 * n1_above)
  discordant <- sum(n0 * n1_below)
  if (concordant + discordant == 0) {
    return(NA_real_)
  }
  (concordant - discordant) / (concordant + discordant)
}

boot_cor <- function(x, y, n_boot = 2000) {
  ok <- !is.na(x) & !is.na(y)
  x <- x[ok]
  y <- y[ok]
  n <- length(x)
  boots <- replicate(n_boot, {
    idx <- sample.int(n, replace = TRUE)
    cor(x[idx], y[idx])
  })
  c(
    r = cor(x, y),
    ci_lo = quantile(boots, 0.025, names = FALSE),
    ci_hi = quantile(boots, 0.975, names = FALSE)
  )
}
