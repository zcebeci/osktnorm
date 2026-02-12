\name{int}
\alias{int}
\title{Inverse Normal Transformation}
\description{
Performs a rank-based inverse normal transformation (INT) that maps a numeric
vector to approximately standard normal scores.
}
\usage{
int(x, ties.method = "average", na.action = "keep")
}
\arguments{
  \item{x}{A numeric vector to be transformed. Missing values (NA) are allowed.}
  \item{ties.method}{A character string specifying how ties are handled in ranking.
    Default is \code{"average"}. See \code{\link[base]{rank}} for available options.}
  \item{na.action}{A character string indicating how missing values are treated
    during ranking. Default is \code{"keep"}.}
}
\value{
A list with one element:
\describe{
  \item{transformed}{A numeric vector of the same length as \code{x} containing
    inverse normal transformed values. Missing values are preserved.}
}
}
\details{
The inverse normal transformation (INT) is a nonparametric normalization method
based on data ranks. The procedure first maps ranks to uniform quantiles and then
applies the inverse standard normal distribution function.

For each non-missing observation \eqn{x_i}, the transformation is defined as:
\eqn{Z_i = qnorm((rank(x_i) - 0.5) / n)}
where \eqn{n} is the number of non-missing observations.
}
\examples{
x <- c(5, 2, 8, 8, 3)
res <- int(x)
res$transformed

# With missing values
x2 <- c(1.2, NA, 3.4, 2.1)
int(x2)$transformed
}
\seealso{
\code{\link[stats]{qnorm}}, \code{\link[base]{rank}}
}
\keyword{normalization}
\keyword{distribution}
\keyword{normality}
