# Gel-Gastwirth Robust Jarque-Bera Test (RJB)
rjbtest <- function(x) {

  x <- as.numeric(x)
  x <- x[!is.na(x)]
  n <- length(x)

  if (n < 8) {
    stop("RJB test requires at least 8 observations")
  }

  # Quantiles
  q <- quantile(
    x,
    probs = c(0.125, 0.25, 0.375, 0.5,
              0.625, 0.75, 0.875),
    names = FALSE,
    type = 7
  )

  q125 <- q[1]; q25 <- q[2]; q375 <- q[3]
  q50  <- q[4]
  q625 <- q[5]; q75 <- q[6]; q875 <- q[7]

  # Robust skewness (Bowley)
  skew_r <- (q75 + q25 - 2 * q50) / (q75 - q25)

  # Robust kurtosis (Moors, excess)
  kurt_r <- ((q875 - q625) + (q375 - q125)) /
            (q75 - q25) - 3

  # Robust Jarque-Bera statistic
  RJB <- (n / 6) * skew_r^2 +
         (n / 24) * kurt_r^2

  pval <- 1 - pchisq(RJB, df = 2)

  structure(
    list(
      statistic = c(RJB = RJB),
      p.value = pval,
      method = "Gel-Gastwirth Robust Jarque-Bera Test",
      data.name = deparse(substitute(x))
    ),
    class = "htest"
  )
}
