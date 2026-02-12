\name{abc}
\alias{abc}
\title{Adaptive Box-Cox (ABC) Normalization}
\description{
Adaptive Box-Cox (ABC) transformation (Yu et al, 2022) is a data transformation method to transform a non-normal 
numeric variable toward normality by tuning a power parameter based on a normality test result. 
The method selects the optimal transformation parameter by maximizing the Shapiro-Wilk normality test p-value.
}

\usage{
abc(x, lrange = seq(-3, 3, 0.01))
}

\arguments{
  \item{x}{
  A numeric vector containing the data to be transformed. If non-positive values
  are present, the data are automatically shifted to ensure positivity.
  }

  \item{lrange}{
  A numeric vector specifying the range of lambda values to be evaluated.
  The default is a sequence from -3 to 3 with step size 0.01.
  }
}

\details{
The Approximate Box-Cox (ABC) method searches over a predefined range of
transformation parameters (\eqn{\lambda}) and applies a Box-Cox-type transformation
for each candidate value.

For each \eqn{\lambda}, the transformed data are tested for normality using the
Shapiro-Wilk test. The optimal \eqn{\lambda} is selected as the value that maximizes
the logarithm of the Shapiro-Wilk p-value.

The transformation is defined as:
\deqn{
x^{(\lambda)} =
\begin{cases}
\log(x), & \lambda = 0 \\
(x^\lambda - 1)/\lambda, & \lambda \neq 0
\end{cases}
}

In this implementation, to satisfy the Box-Cox requirement of strictly positive data, the function
automatically shifts the input vector if any non-positive values are detected.
}

\value{
A list with the following components:
\item{transformed}{A numeric vector containing the transformed data.}
\item{lambda}{The selected lambda value that maximizes the Shapiro--Wilk log p-value.}
}

\references{
Yu, H., Sang, P., & Huan, T. (2022). Adaptive Box-Cox transformation: A highly flexible feature-specific 
data transformation to improve metabolomic data normality for better statistical analysis. 
\emph{Analytical Chemistry, 94}(23), 8267-8276. \doi{10.1021/acs.analchem.2c00503}.
}

\examples{
set.seed(123)
x <- rexp(100)

result <- abc(x)
result$lambda
hist(result$transformed, main = "ABC Transformed Data")
}

\author{
  Zeynel Cebeci
}

\keyword{normality transformation}

