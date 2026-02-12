\name{oskt}
\alias{oskt}
\title{Normalization via Skewness and Kurtosis Minimization of Anderson--Darling Test Statistic}
\description{
Applies the Tukey g-and-h transformation to a numeric vector and estimates optimal skewness and tail-heaviness parameters
by minimizing the Anderson-Darling normality test statistic (A\eqn{^2}).
}

\usage{
oskt(x, init_params = c(0.1, 0.1), lower_bounds = c(-1, 0), upper_bounds = c(1, 0.5))
}

\arguments{
  \item{x}{
  A numeric vector of observations to be transformed.
  The data are internally standardized to zero mean and unit variance.
  }

  \item{init_params}{
  Numeric vector of length two giving the initial values of the Tukey g-and-h parameters \code{(g, h)}.
  }

  \item{lower_bounds}{
  Numeric vector of length two specifying lower bounds for \code{g} and \code{h}.
  }

  \item{upper_bounds}{
  Numeric vector of length two specifying upper bounds for \code{g} and \code{h}.
  }
}

\details{
The Tukey g-and-h transformation is defined as:
\deqn{
Y =
\begin{cases}
\dfrac{\exp(gX) - 1}{g} \exp\!\left(\dfrac{hX^2}{2}\right), & g \neq 0, \\
X \exp\!\left(\dfrac{hX^2}{2}\right), & g = 0.
\end{cases}
}

where \eqn{g} controls skewness and \eqn{h} controls tail heaviness.

To assess normality of the transformed data, the Anderson--Darling statistic
\eqn{A^2} is computed directly from its analytical form under the standard
normal distribution. The implementation includes the Stephens correction
used when location and scale parameters are estimated from the data:
\deqn{
A^{2*} = A^2 \left(1 + \frac{0.75}{n} + \frac{2.25}{n^2}\right).
}

Optimal parameters are obtained by minimizing the corrected statistic
\eqn{A^{2*}} using constrained optimization via the \code{"L-BFGS-B"} algorithm.
}

\value{
A list with the following components:
\item{transformed}{Numeric vector containing the transformed data.}
\item{g}{Estimated skewness parameter of the Tukey g-and-h transformation.}
\item{h}{Estimated tail-heaviness parameter of the Tukey g-and-h transformation.}
}

\note{
This function does not rely on external normality testing packages.
The Anderson--Darling statistic is computed explicitly to ensure
numerical consistency and package independence.

If optimization fails, the standardized input data are returned and
\code{g} and \code{h} are set to \code{NA}.
}

\references{
Stephens, M. A. (1974). EDF statistics for goodness of fit and some comparisons.
\emph{Journal of the American Statistical Association}, 69(347), 730-737.

Tukey, J. W. (1977). \emph{Exploratory Data Analysis}. Addison-Wesley.
}

\examples{
set.seed(123)
x <- rexp(100)

restrans <- oskt(x)
restrans$g
restrans$h

hist(restrans$transformed, xlab="Transformed Values", main = "OSKT Transformed Data")

set.seed(123)
x <- rlnorm(300, meanlog = 0, sdlog = 0.75)

restrans <- oskt(x)
restrans$g
restrans$h

hist(restrans$transformed, xlab="Transformed Values", main = "OSKT Transformed Data")
}

\author{
 Zeynel Cebeci
}
