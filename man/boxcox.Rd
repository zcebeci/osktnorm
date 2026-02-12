\name{boxcox}
\alias{boxcox}
\title{Box-Cox Transformation}
\description{
Performs a Box-Cox transformation on a numeric vector. Optionally, the data 
can be shifted to ensure all values are positive before applying the transformation. 
If lambda is not provided, it is estimated via maximum likelihood.
}
\usage{
boxcox(x, lambda = NULL, makepositive = FALSE, eps = 1e-06)
}
\arguments{
  \item{x}{A numeric vector to be transformed.}

  \item{lambda}{Optional numeric value of the Box--Cox transformation parameter.
  If \code{NULL} (default), the value of \code{lambda} is estimated by maximizing
  the profile log-likelihood.}

  \item{makepositive}{Logical. If \code{TRUE}, the data are shifted so that all
  observations are strictly positive before applying the transformation.
  This is required when \code{x} contains zero or negative values.
  Default is \code{FALSE}.}

  \item{eps}{A small positive constant used for numerical stability.
  It is added implicitly when enforcing positivity or to avoid taking the logarithm
  of zero. Default is \code{1e-06}.}
}
\details{
The Box-Cox transformation is defined as:

\deqn{
y(\lambda) = 
\begin{cases} 
\frac{x^\lambda - 1}{\lambda}, & \lambda \neq 0 \\
\log(x), & \lambda = 0
\end{cases}
}

If \code{makepositive = TRUE}, the function shifts the data by \code{abs(min(x)) + 1} 
if there are zero or negative values, to make all values positive.  

Lambda is estimated via maximum likelihood over a grid of values from -4 to 4 
(default 500 points) if not specified.
}
\value{
A list with the following components:

\item{transformed}{The transformed numeric vector.}
\item{lambda}{The lambda used in the transformation (either provided or estimated).}
\item{shift}{The amount by which the original data was shifted to make it positive (0 if no shift).}
}

\seealso{
\code{\link{abc}}, \code{\link{rewbc}}
}

\author{
 Zeynel Cebeci
}

\examples{
# Generate positively skewed example data (n > 30)
set.seed(123)
x <- rlnorm(50, meanlog = 0, sdlog = 1)

# Box-Cox with estimated lambda (MLE)
res <- boxcox(x)
head(res$transformed)
res$lambda
res$shift

# Box-Cox with specified lambda
res2 <- boxcox(x, lambda = 0.25)
head(res2$transformed)

# Box-Cox with automatic shift for nonpositive values
x2 <- x - quantile(x, 0.2)
res3 <- boxcox(x2, makepositive = TRUE)
head(res3$transformed)
res3$shift
}

\keyword{normality transformation}
\keyword{power transformation}
