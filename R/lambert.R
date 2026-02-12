lambert <- function(x,
                    type = c("s", "h", "hh"),
                    standardize = TRUE,
                    tol = 1e-8,
                    maxiter = 100) {

  type <- match.arg(type)
  stopifnot(is.numeric(x))

  na_idx <- is.na(x)
  x <- x[!na_idx]
  n <- length(x)
  if (n < 5) stop("Insufficient data")

  ## Location / scale
  mu <- mean(x)
  sigma <- sd(x)
  z <- (x - mu) / sigma

  ## Lambert W (principal branch, safe)
  lambertW0 <- function(x) {
    w <- pmax(x, -1 / exp(1) + 1e-12)
    for (i in seq_len(maxiter)) {
      ew <- exp(w)
      w_new <- w - (w * ew - x) / (ew * (w + 1))
      if (max(abs(w_new - w), na.rm = TRUE) < tol) break
      w <- w_new
    }
    w
  }

  ## ---------- TYPE = "s" (skewness only) ----------
  if (type == "s") {

    delta_max <- min(abs(1 / (exp(1) * z[z != 0])))
    delta_max <- min(delta_max, 1)

    obj_s <- function(delta) {
      arg <- delta * z
      if (any(arg < -1 / exp(1))) return(1e10)
      u <- if (abs(delta) < 1e-6) z else lambertW0(arg) / delta
      u <- u - mean(u)
      abs(mean(u^3))
    }

    opt <- optimize(obj_s, c(-delta_max, delta_max))
    delta <- opt$minimum
    h <- 0

    u <- if (abs(delta) < 1e-6) z else lambertW0(delta * z) / delta
  }

  ## ---------- TYPE = "h" (kurtosis only) ----------
  if (type == "h") {

    obj_h <- function(h) {
      if (!is.finite(h) || h < 0) return(1e10)
      arg <- h * z^2
      u <- sign(z) * sqrt(lambertW0(arg) / h)
      u <- u - mean(u)
      abs(mean(u^4) - 3)
    }

    opt <- optimize(obj_h, c(1e-6, 0.5))
    h <- opt$minimum
    delta <- 0

    u <- sign(z) * sqrt(lambertW0(h * z^2) / h)
  }

  ## ---------- TYPE = "hh" (skewness + kurtosis) ----------
  if (type == "hh") {

    obj_hh <- function(par) {
      delta <- par[1]
      h <- par[2]
      if (!is.finite(delta) || !is.finite(h) || h < 0) return(1e10)

      u <- numeric(length(z))
      for (i in seq_along(z)) {
        f <- function(ui)
          ui * exp(delta * ui + 0.5 * h * ui^2) - z[i]
        rt <- try(uniroot(f, c(-8, 8))$root, silent = TRUE)
        if (inherits(rt, "try-error") || !is.finite(rt))
          return(1e10)
        u[i] <- rt
      }

      u <- u - mean(u)
      m3 <- mean(u^3)
      m4 <- mean(u^4)
      if (!is.finite(m3) || !is.finite(m4)) return(1e10)

      abs(m3) + abs(m4 - 3)
    }

    opt <- optim(
      par = c(0, 0.05),
      fn = obj_hh,
      method = "L-BFGS-B",
      lower = c(-0.3, 0),
      upper = c(0.3, 0.5)
    )

    delta <- opt$par[1]
    h <- opt$par[2]

    u <- sapply(z, function(zi) {
      f <- function(ui)
        ui * exp(delta * ui + 0.5 * h * ui^2) - zi
      uniroot(f, c(-8, 8))$root
    })
  }

  if (standardize)
    u <- as.numeric(scale(u))

  x_trans <- rep(NA_real_, length(na_idx))
  x_trans[!na_idx] <- u

  list(
    transformed = x_trans,
    delta = delta,
    h = h,
    mean = mu,
    sd = sigma,
    n = n,
    type = type,
    standardize = standardize
  )
}
