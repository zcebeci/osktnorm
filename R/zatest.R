# Zhang-Wu ZA Test (Zhang & Wu, 2005)
# =================================================================
zatest <- function(x, nsim = 10000, eps = 1e-10, ncores = 1, seed = NULL) {

  if (!is.null(seed))
    set.seed(seed)

  # ---------------------------------------------------------------
  # Validate ncores (CRAN-safe)
  # ---------------------------------------------------------------
  if (!is.numeric(ncores) || length(ncores) != 1 || is.na(ncores))
    stop("`ncores` must be a single numeric value.", call. = FALSE)

  ncores <- as.integer(ncores)

  if (ncores < 1)
    stop("`ncores` must be >= 1.", call. = FALSE)

  max_cores <- parallel::detectCores(logical = TRUE)

  if (ncores > max_cores)
    ncores <- max_cores

  # ---------------------------------------------------------------
  # Data checks
  # ---------------------------------------------------------------
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  n <- length(x)

  if (n < 5)
    stop("Sample size must be at least 5.")

  i <- seq_len(n)

  # ---------------------------------------------------------------
  # ZA statistic
  # ---------------------------------------------------------------
  ZA_stat <- function(z) {
    z <- sort(z)
    p <- pnorm(z)
    p <- pmin(pmax(p, eps), 1 - eps)

    -sum(
      log(p) / (n - i + 0.5) +
      log(1 - p) / (i - 0.5)
    )
  }

  # ---------------------------------------------------------------
  # Observed statistic
  # ---------------------------------------------------------------
  z_obs <- (x - mean(x)) / sd(x)
  ZA_obs <- ZA_stat(z_obs)

  # ---------------------------------------------------------------
  # Monte Carlo simulation
  # ---------------------------------------------------------------
  if (ncores > 1) {

    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl), add = TRUE)

    `%loop%` <- foreach::`%dopar%`

  } else {

    `%loop%` <- foreach::`%do%`
  }

  ZA_sim <- foreach::foreach(
    b = seq_len(nsim),
    .combine = c,
    .packages = "stats"
  ) %loop% {

    y <- rnorm(n)
    z_sim <- (y - mean(y)) / sd(y)
    ZA_stat(z_sim)
  }

  # ---------------------------------------------------------------
  # Monte Carlo p-value
  # ---------------------------------------------------------------
  pval <- (sum(ZA_sim >= ZA_obs) + 1) / (nsim + 1)

  structure(
    list(
      statistic = c(ZA = ZA_obs),
      p.value   = pval,
      method    = "Zhang-Wu ZA Test for Normality",
      data.name = deparse(substitute(x))
    ),
    class = "htest"
  )
}
