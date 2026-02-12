# Yeo-Johnson Tranformation
yeojohnson <- function(x, lambda = NULL, standardize = TRUE, eps = 1e-6) {

  x <- as.numeric(x)
  na_idx <- is.na(x)
  x <- x[!na_idx]

  n <- length(x)
  if (n < 5 || all(x == x[1]))
    stop("No Yeo-Johnson applied: insufficient or identical values")

  pos <- x >= 0
  neg <- !pos

  # ----------------------------
  # Yeo-Johnson transform
  # ----------------------------
  yj_trans <- function(lambda) {

    xt <- numeric(n)

    if (any(pos)) {
      xp <- x[pos] + 1
      xt[pos] <- if (abs(lambda) < eps)
        log(xp)
      else
        (xp^lambda - 1) / lambda
    }

    if (any(neg)) {
      xn <- -x[neg] + 1
      if (abs(lambda - 2) < eps) {
        xt[neg] <- -log(xn)
      } else {
        xt[neg] <- -((xn^(2 - lambda) - 1) / (2 - lambda))
      }
    }

    xt
  }

  # ----------------------------
  # Log-likelihood
  # ----------------------------
  if (is.null(lambda)) {

    c_pos <- if (any(pos)) sum(log(x[pos] + 1)) else 0
    c_neg <- if (any(neg)) sum(log(-x[neg] + 1)) else 0

    loglik <- function(lambda) {

      xt <- yj_trans(lambda)
      mu <- mean(xt)
      s2 <- mean((xt - mu)^2)

      if (!is.finite(s2) || s2 <= 0) return(-Inf)

      -n / 2 * log(s2) +
        (lambda - 1) * c_pos +
        (1 - lambda) * c_neg
    }

    lambda_grid <- seq(-3, 3, by = 0.02)
    ll <- vapply(lambda_grid, loglik, numeric(1))

    lambda <- lambda_grid[which.max(ll)]
  }

  # ----------------------------
  # Final transform
  # ----------------------------
  x_trans <- yj_trans(lambda)

  if (standardize) {
    s <- sd(x_trans)
    if (s > 0)
      x_trans <- (x_trans - mean(x_trans)) / s
  }

  out <- rep(NA_real_, length(na_idx))
  out[!na_idx] <- x_trans

  list(
    transformed = out,
    lambda = lambda,
    n = n,
    method = "Yeo-Johnson Transformation",
    standardize = standardize
  )
}


