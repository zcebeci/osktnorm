# Cramer-von Mises Normality Test (Monte Carlo based approach)
# Instead of using asymptotic tables or approximate p-values, 
# a Monte Carlo reference distribution is generated.
# =========================================================
cvmtest <- function(x, nsim = 10000, seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  x <- as.numeric(x)
  x <- x[is.finite(x)]
  n <- length(x)

  if (n < 5)
    stop("Sample size must be at least 5.")

  # Precompute constants
  i <- seq_len(n)
  a <- (2 * i - 1) / (2 * n)
  c0 <- 1 / (12 * n)

  # Observed statistic
  z <- sort((x - mean(x)) / sd(x))
  Fz <- pnorm(z)
  W2_obs <- c0 + sum((Fz - a)^2)

  # Monte Carlo reference
  W2_sim <- numeric(nsim)

  for (b in seq_len(nsim)) {

    y <- rnorm(n)
    z_sim <- sort((y - mean(y)) / sd(y))
    Fy <- pnorm(z_sim)

    W2_sim[b] <- c0 + sum((Fy - a)^2)
  }

  # Unbiased Monte Carlo p-value
  p_val <- (sum(W2_sim >= W2_obs) + 1) / (nsim + 1)

  structure(
    list(
      statistic = c(W2 = W2_obs),
      p.value   = p_val,
      method    = "Cramer-von Mises Test for Normality",
      data.name = deparse(substitute(x))
    ),
    class = "htest"
  )
}
