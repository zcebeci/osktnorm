# Inverse Normal Transformation
int <- function(x, ties.method = "average", na.action = "keep") {
  stopifnot(is.numeric(x))

  # Compute ranks, preserving NAs
  ranks <- rank(x, ties.method = ties.method, na.last = na.action)

  # Number of non-missing observations
  n <- sum(!is.na(x))

  # Map ranks to uniform quantiles and then to normal quantiles
  p <- (ranks - 0.5) / n
  z <- qnorm(p)

  return(list(transformed = z))
}
