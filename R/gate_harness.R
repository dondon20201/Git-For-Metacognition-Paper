correlated_pair <- function(n, rho) {
  z <- matrix(rnorm(n * 2), n) %*% chol(matrix(c(1, rho, rho, 1), 2))
  list(a = z[, 1], b = z[, 2])
}

draw_correlated_traits <- function(n, rho_level, rho_sens, rho_dprime) {
  zb <- correlated_pair(n, rho_level)
  zl <- correlated_pair(n, rho_sens)
  zd <- correlated_pair(n, rho_dprime)
  mk <- function(z_beta, z_lambda, z_d) data.frame(
    d_prime = pmax(1.2 + 0.3 * z_d, 0.2),
    lambda = 0.15 + 0.7 * pnorm(z_lambda),
    beta = 0.8 * z_beta
  )
  list(a = mk(zb$a, zl$a, zd$a), b = mk(zb$b, zl$b, zd$b))
}

sim_domain <- function(traits, n_trials, k, prefix) {
  res <- lapply(seq_len(nrow(traits)), function(i) {
    tr <- simulate_trials(
      n_trials, traits$d_prime[i], traits$lambda[i], traits$beta[i]
    )
    c(parcel_measures(tr, k), dhat = d_prime_hat(tr))
  })
  m <- do.call(rbind, res)
  colnames(m) <- c(
    paste0("s", prefix, seq_len(k)), paste0("l", prefix, seq_len(k)),
    paste0("dhat", prefix)
  )
  as.data.frame(m)
}

partial_out_dprime <- function(dat, k, prefix) {
  d <- dat[[paste0("dhat", prefix)]]
  for (j in seq_len(k)) {
    col <- paste0("s", prefix, j)
    dat[[col]] <- resid(lm(dat[[col]] ~ d, na.action = na.exclude))
  }
  dat
}

k_for <- function(n_trials) if (n_trials >= 200) 4L else 2L

gate_replicate <- function(n_subjects, rhos, n_trials_a, n_trials_b, seed) {
  set.seed(seed)
  tr <- draw_correlated_traits(n_subjects, rhos[["level"]], rhos[["sens"]],
                               rhos[["dprime"]])
  k_a <- k_for(n_trials_a)
  k_b <- k_for(n_trials_b)
  dat <- cbind(sim_domain(tr$a, n_trials_a, k_a, "A"),
               sim_domain(tr$b, n_trials_b, k_b, "B"))
  dat <- partial_out_dprime(partial_out_dprime(dat, k_a, "A"), k_b, "B")
  cbind(
    data.frame(
      data_type = "SIMULATION", seed = seed,
      n_trials_a = n_trials_a, n_trials_b = n_trials_b,
      realized_level = cor(tr$a$beta, tr$b$beta),
      realized_sens = cor(tr$a$lambda, tr$b$lambda),
      realized_dprime = cor(tr$a$d_prime, tr$b$d_prime)
    ),
    fit_latent(dat, k_a, k_b)
  )
}
