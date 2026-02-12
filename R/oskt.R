# Transformation with Skewness and Kurtosis Minimization
oskt <- function(x,
                      init_params = c(0.1, 0.1),
                      lower_bounds = c(-1, 0),
                      upper_bounds = c(1, 0.5)) {

  x <- as.numeric(x)
  if (length(x) < 5 || all(x == x[1]))
    stop("Insufficient or constant data")

  # -------------------------------------------------------
  # Pre-standardize once
  # -------------------------------------------------------
  x <- (x - mean(x)) / sd(x)
  n <- length(x)
  i <- seq_len(n)

  # Stephens correction
  stephens <- (1 + 0.75 / n + 2.25 / n^2)

  eps <- 1e-15

  # -------------------------------------------------------
  # Fast Anderson-Darling A2 wth Stephens correction
  # -------------------------------------------------------
  compA2_fast <- function(y) {

    z <- (y - mean(y)) / sd(y)
    z <- sort(z)

    Fi <- pnorm(z)
    Fi <- pmin(pmax(Fi, eps), 1 - eps)

    A2 <- -n - mean((2 * i - 1) *
              (log(Fi) + log(1 - Fi[n:1])))

    A2 * stephens
  }

  # -------------------------------------------------------
  # g-h transform (branchless style)
  # -------------------------------------------------------
  gh_transform <- function(g, h) {
    if (abs(g) > 1e-8) {
      (exp(g * x) - 1) / g * exp(0.5 * h * x^2)
    } else {
      x * exp(0.5 * h * x^2)
    }
  }

  objfun <- function(par) {
    compA2_fast(gh_transform(par[1], par[2]))
  }

  opt <- try(
    optim(
      par = init_params,
      fn = objfun,
      method = "L-BFGS-B",
      lower = lower_bounds,
      upper = upper_bounds
    ),
    silent = TRUE
  )

  if (inherits(opt, "try-error")) {
    warning("Optimization failed; returning original data")
    return(list(transformed = x, g = NA, h = NA))
  }

  g_opt <- opt$par[1]
  h_opt <- opt$par[2]

  list(
    transformed = gh_transform(g_opt, h_opt),
    g = g_opt,
    h = h_opt
  )
}
