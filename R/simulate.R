folded_mean <- function(mu) {
  sqrt(2 / pi) * exp(-mu^2 / 2) + mu * (1 - 2 * pnorm(-mu))
}

folded_sd <- function(mu) {
  sqrt(mu^2 + 1 - folded_mean(mu)^2)
}

conf_thresholds <- qnorm(seq_len(6) / 7)

simulate_trials <- function(n_trials, d_prime, lambda, beta) {
  mu <- d_prime / 2
  signal <- rbinom(n_trials, 1, 0.5)
  x <- rnorm(n_trials, mean = mu * (2 * signal - 1), sd = 1)
  response <- as.integer(x > 0)
  accuracy <- as.integer(response == signal)
  z <- lambda * (abs(x) - folded_mean(mu)) / folded_sd(mu) +
    sqrt(1 - lambda^2) * rnorm(n_trials) + beta
  confidence <- 1L + findInterval(z, conf_thresholds)
  data.frame(signal, response, accuracy, confidence)
}

draw_traits <- function(n_subjects, d_prime_mean = 1.2) {
  data.frame(
    d_prime = pmax(rnorm(n_subjects, d_prime_mean, 0.3), 0.2),
    lambda = runif(n_subjects, 0.15, 0.85),
    beta = rnorm(n_subjects, 0, 0.8)
  )
}

run_cohort <- function(traits, n_trials) {
  est <- lapply(seq_len(nrow(traits)), function(i) {
    tr <- simulate_trials(
      n_trials, traits$d_prime[i], traits$lambda[i], traits$beta[i]
    )
    halves <- list(all = tr, odd = tr[c(TRUE, FALSE), ], even = tr[c(FALSE, TRUE), ])
    out <- lapply(halves, function(h) {
      c(
        mean_conf = mean(h$confidence),
        accuracy = mean(h$accuracy),
        auroc2 = auroc2(h$confidence, h$accuracy),
        gamma = gk_gamma(h$confidence, h$accuracy)
      )
    })
    unlist(out)
  })
  cbind(traits, as.data.frame(do.call(rbind, est)))
}
