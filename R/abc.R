# Adaptive Box-Cox (ABC) transformation - Yu et al (2022)
# ============================================================

abc <- function(x, lrange = seq(-3, 3, 0.01)) {

  x <- as.numeric(x)
  x <- x[is.finite(x)]

  n <- length(x)
  if (n < 5 || all(x == x[1]))
    stop("No ABC applied: insufficient or identical values")

  # Positive support (scale + shift)
  if (any(x <= 0)) {
    x <- x - min(x)
    x <- x / sd(x)
    x <- x + 1
  }

  # Objective: log Shapiro p-value (clipped)
  obj <- function(lambda) {

    xt <- if (abs(lambda) < 1e-6)
      log(x)
    else
      (x^lambda - 1) / lambda

    p <- tryCatch(
      shapiro.test(xt)$p.value,
      error = function(e) NA_real_
    )

    if (is.na(p) || p <= 0) return(-Inf)
    log(p)
  }

  logp <- vapply(lrange, obj, numeric(1))

  if (all(!is.finite(logp)))
    stop("No ABC applied: Shapiro failed for all lambdas")

  lambda_hat <- lrange[which.max(logp)]

  x_trans <- if (abs(lambda_hat) < 1e-6)
    log(x)
  else
    (x^lambda_hat - 1) / lambda_hat

  list(
    transformed = as.numeric(scale(x_trans)),
    lambda = lambda_hat,
    n = n,
    method = "Adaptive Box-Cox Transformation"
  )
}

