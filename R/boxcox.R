# ============================================================
# Box-Cox Transformation
# ============================================================

boxcox <- function(x,
                         lambda = NULL,
                         makepositive = FALSE,
                         eps = 1e-6) {

  x <- as.numeric(x)
  x <- x[is.finite(x)]

  n <- length(x)
  if (n < 5 || all(x == x[1]))
    stop("Box-Cox not applicable: insufficient or identical values")

  # ----------------------------
  # Positivity handling
  # ----------------------------
  shift <- 0
  if (makepositive && any(x <= 0)) {
    shift <- abs(min(x)) + 1
    warning("Data shifted to enforce positivity; lambda interpretation affected")
  }

  x0 <- x + shift
  if (any(x0 <= 0))
    stop("Box-Cox requires strictly positive data")

  logx <- log(x0)
  sumlogx <- sum(logx)

  # ----------------------------
  # Box-Cox transform
  # ----------------------------
  bc_trans <- function(lambda) {
    if (abs(lambda) < eps)
      logx
    else
      (x0^lambda - 1) / lambda
  }

  # ----------------------------
  # Log-likelihood
  # ----------------------------
  if (is.null(lambda)) {

    loglik <- function(lambda) {
      xt <- bc_trans(lambda)
      mu <- mean(xt)
      s2 <- mean((xt - mu)^2)
      if (!is.finite(s2) || s2 <= 0) return(-Inf)
      -n / 2 * log(s2) + (lambda - 1) * sumlogx
    }

    grid <- seq(-4, 4, by = 0.02)
    ll <- vapply(grid, loglik, numeric(1))
    lambda <- grid[which.max(ll)]
  }

  transformed <- bc_trans(lambda)

  list(
    transformed = transformed,
    lambda = lambda,
    shift = shift,
    n = n,
    method = "Box-Cox Transformation"
  )
}
