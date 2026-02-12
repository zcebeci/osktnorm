# Reweighted Maximum Likelihood Box-Cox Transformation
# Raymaekers & Rousseeuw (2024)
# ================================================================================

rewbc<- function(x, lrange = seq(-3, 3, by = 0.01), rwsteps = 2, k = 1.5) {

  x <- as.numeric(x)
  x <- x[is.finite(x)]

  if (length(x) < 3 || all(x == x[1]))
    stop("RBC not applicable: identical or insufficient values")

  # Positivity
  if (any(x <= 0))
    x <- x - min(x) + 1e-6

  n <- length(x)
  logx <- log(x)

  bc_transform <- function(lambda) {
    if (abs(lambda) < 1e-8) logx else (x^lambda - 1) / lambda
  }

  bc_loglik <- function(lambda, w = NULL) {

    xt <- bc_transform(lambda)

    if (is.null(w)) {
      mu <- mean(xt)
      sigma2 <- mean((xt - mu)^2)
    } else {
      sw <- sum(w)
      mu <- sum(w * xt) / sw
      sigma2 <- sum(w * (xt - mu)^2) / sw
    }

    sigma2 <- max(sigma2, .Machine$double.eps)

    -(n / 2) * log(sigma2) +
      (lambda - 1) * sum(logx)
  }

  # Initial MLE
  ll_vals <- vapply(lrange, bc_loglik, numeric(1))
  lambda_hat <- lrange[which.max(ll_vals)]

  # Reweighting
  w <- rep(1, n)

  for (step in seq_len(rwsteps)) {

    xt <- bc_transform(lambda_hat)

    mu <- median(xt)
    s  <- mad(xt, constant = 1)
    if (s == 0) s <- sd(xt)

    z <- (xt - mu) / s
    w <- ifelse(abs(z) <= k, 1, k / abs(z))

    ll_vals <- vapply(lrange, bc_loglik, numeric(1), w = w)
    lambda_hat <- lrange[which.max(ll_vals)]
  }

  list(
    transformed = bc_transform(lambda_hat),
    lambda = lambda_hat,
    weights = w,
    steps = rwsteps
  )
}
