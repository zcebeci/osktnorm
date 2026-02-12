# Pearson P Metric for Normality (scaled Chi-square)
# =================================================================
pearsonp <- function(x, nbins = NULL) {

  x <- na.omit(x)
  n <- length(x)

  if (n < 5)
    stop("Number of observations should be greater than 5.")

  # Number of classes (nortest default)
  if (is.null(nbins)) {
    k <- floor(n^(2/5))
  } else {
    k <- as.integer(nbins)
  }

  if (k < 3)
    stop("Number of bins must be at least 3.")

  # Degrees of freedom: k - 1 - 2 (mean and variance)
  df <- k - 3

  # Standardize
  z <- (x - mean(x)) / sd(x)

  # Equal-probability bins
  probs  <- seq(0, 1, length.out = k + 1)
  breaks <- qnorm(probs)

  obs <- as.numeric(
    table(cut(z, breaks = breaks, include.lowest = TRUE))
  )
  exp <- rep(n / k, k)

  # Pearson chi-square
  P <- sum((obs - exp)^2 / exp)

  # Scaled Pearson P metric
  Pmetric <- P / df

  structure(list(
    statistic = c(PearsonP = Pmetric),
    method = "Pearson P Metric for Normality (scaled Chi-square)",
    data.name = deparse(substitute(x)),
    df = df
  ), class = "htest")
}
