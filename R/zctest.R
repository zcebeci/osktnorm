# =========================================================================================
# Zhang-Wu ZC Test for Normality
# Zhang (2002), Zhang & Wu (2005)
# Monte Carlo p-value with parameter re-estimation
# =========================================================================================

zctest <- function(x, nsim = 10000, eps = 1e-10, ncores = 1, seed = NULL) {

  if (!is.null(seed))
    set.seed(seed)

  # -------------------------------------------------------
  # Validate ncores (CRAN-safe)
  # -------------------------------------------------------
  if (!is.numeric(ncores) || length(ncores) != 1 || is.na(ncores))
    stop("`ncores` must be a single numeric value.", call. = FALSE)

  ncores <- as.integer(ncores)

  if (ncores < 1)
    stop("`ncores` must be >= 1.", call. = FALSE)

  max_cores <- parallel::detectCores(logical = TRUE)

  if (ncores > max_cores)
    ncores <- max_cores

  # -------------------------------------------------------
  # Data checks
  # -------------------------------------------------------
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  n <- length(x)

  if (n < 5)
    stop("Number of observations must be at least 5.")

  if (all(x == x[1]))
    stop("ZC test not applicable: identical values.")

  i <- seq_len(n)

  # -------------------------------------------------------
  # Observed ZC statistic
  # -------------------------------------------------------
  z <- sort((x - mean(x)) / sd(x))
  p <- pnorm(z)
  p <- pmin(pmax(p, eps), 1 - eps)

  T1 <- log(
    (1 / p - 1) /
    ((n - 0.5) / (i - 0.75) - 1)
  )

  ZC_obs <- sum(T1^2)

  # -------------------------------------------------------
  # Monte Carlo reference distribution
  # -------------------------------------------------------
  if (ncores > 1) {

    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl), add = TRUE)

    `%loop%` <- foreach::`%dopar%`

  } else {

    `%loop%` <- foreach::`%do%`
  }

  ZC_sim <- foreach::foreach(
    b = seq_len(nsim),
    .combine = c,
    .packages = "stats"
  ) %loop% {

    y <- rnorm(n)

    z_sim <- sort((y - mean(y)) / sd(y))
    p_sim <- pnorm(z_sim)
    p_sim <- pmin(pmax(p_sim, eps), 1 - eps)

    T1_sim <- log(
      (1 / p_sim - 1) /
      ((n - 0.5) / (i - 0.75) - 1)
    )

    sum(T1_sim^2)
  }

  # -------------------------------------------------------
  # Monte Carlo p-value (right tail, unbiased)
  # -------------------------------------------------------
  pval <- (sum(ZC_sim >= ZC_obs) + 1) / (nsim + 1)

  structure(
    list(
      statistic = c(ZC = ZC_obs),
      p.value   = pval,
      method    = "Zhang-Wu ZC Test for Normality",
      data.name = deparse(substitute(x))
    ),
    class = "htest"
  )
}
