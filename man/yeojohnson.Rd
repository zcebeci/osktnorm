\name{yeojohnson}
\alias{yeojohnson}
\title{Yeo-Johnson Transformation}
\description{
Performs a Yeo-Johnson transformation (Yeo & Johnson, 2000) on a numeric vector. 
The transformation estimates the optimal lambda via maximum likelihood if not provided. 
Optionally, the transformed data can be standardized.
}
\usage{
yeojohnson(x, lambda = NULL, standardize = TRUE, eps = 1e-6)
}
\arguments{
  \item{x}{A numeric vector to transform.}
  \item{lambda}{Optional numeric value of lambda for the Yeo-Johnson transformation. 
                If \code{NULL} (default), lambda is estimated via maximum likelihood.}
  \item{standardize}{Logical. If \code{TRUE} (default), the transformed values are centered 
                     and scaled to have mean 0 and standard deviation 1.}
  \item{eps}{Numeric tolerance used to handle cases where lambda is approximately 0 or 2. 
              Default is \code{1e-6}.}
}
\details{
The Yeo-Johnson transformation is a generalization of the Box-Cox transformation 
that can handle both positive and negative values:

\deqn{
y(\lambda) = 
\begin{cases} 
((x + 1)^\lambda - 1) / \lambda, & x \ge 0, \lambda \neq 0 \\
\log(x + 1), & x \ge 0, \lambda = 0 \\
-((-x + 1)^{2 - \lambda} - 1) / (2 - \lambda), & x < 0, \lambda \neq 2 \\
-\log(-x + 1), & x < 0, \lambda = 2
\end{cases}
}

If \code{lambda} is not specified, it is estimated via maximum likelihood over a 
grid of values from -3 to 3 (step 0.01). Standardization is optional and centers 
the transformed data at mean 0 with standard deviation 1.
}
\value{
A list with the following components:

\item{transformed}{The transformed numeric vector.}
\item{lambda}{The lambda value used in the transformation (either provided or estimated via MLE).}
}

\references{
Yeo, I.-K. and Johnson, R. A. (2000). A new family of power transformations to improve normality or symmetry. 
\emph{Biometrika}, 87(4), 954-959. \doi{10.1093/biomet/87.4.954}.
}

\examples{
x <- c(-2, 0.5, 2, 3, 10, 5, 3, 2)

# Yeo-Johnson with estimated lambda
res <- yeojohnson(x)
res$transformed
res$lambda

# Yeo-Johnson with specified lambda
res2 <- yeojohnson(x, lambda = 0.5)
res2$transformed

# Standardization turned off
res3 <- yeojohnson(x, standardize = FALSE)
res3$transformed
}
\seealso{
\code{\link{rewyj}}
}
