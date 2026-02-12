# Reweighted Maximum Likelihood Yeo-Johnson Transformation
# Raymaekers & Rousseeuw (2024)
# ================================================================================
rewyj <- function(x, lrange = seq(-3, 3, by = 0.01), rwsteps = 2, k = 1.5) {

  x <- as.numeric(x)
  x <- x[is.finite(x)]

  if (length(x) < 3 || all(x == x[1]))
    stop("RYJ not applicable: identical or insufficient values")

  n <- length(x)

  pos <- x >= 0
  neg <- !pos

  logp <- log(x[pos] + 1)
  logn <- log(-x[neg] + 1)

  yj_transform <- function(lambda) {
    out <- numeric(n)

    if (lambda == 0) {
      out[pos] <- logp
      out[neg] <- -logn
    } else if (lambda == 2) {
      out[pos] <- ((x[pos] + 1)^lambda - 1) / lambda
      out[neg] <- -logn
    } else {
      out[pos] <- ((x[pos] + 1)^lambda - 1) / lambda
      out[neg] <- -((-x[neg] + 1)^(2 - lambda) - 1) / (2 - lambda)
    }

    out
  }

  jacobian_term <- function(lambda) {
    sum((lambda - 1) * logp) + sum((1 - lambda) * logn)
  }

  yj_loglik <- function(lambda, w = NULL) {
    xt <- yj_transform(lambda)

    if (is.null(w)) {
      mu <- mean(xt)
      sigma2 <- mean((xt - mu)^2)
    } else {
      mu <- sum(w * xt) / sum(w)
      sigma2 <- sum(w * (xt - mu)^2) / sum(w)
    }

    -(n / 2) * log(sigma2) + jacobian_term(lambda)
  }

  # Initial MLE
  ll_vals <- vapply(lrange, yj_loglik, numeric(1))
  lambda_hat <- lrange[which.max(ll_vals)]

  # Reweighting
  for (step in seq_len(rwsteps)) {

    xt <- yj_transform(lambda_hat)
    mu <- median(xt)
    s  <- mad(xt, constant = 1)
    if (s == 0) s <- sd(xt)

    z <- (xt - mu) / s
    w <- ifelse(abs(z) <= k, 1, k / abs(z))

    ll_vals <- vapply(lrange, yj_loglik, numeric(1), w = w)
    lambda_hat <- lrange[which.max(ll_vals)]
  }

  list(
    transformed = yj_transform(lambda_hat),
    lambda = lambda_hat,
    weights = w,
    steps = rwsteps
  )
}

