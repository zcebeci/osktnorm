\name{cvmtest}
\alias{cvmtest}
\title{Cramér-von Mises Normality Test}
\description{
Performs a Monte Carlo based Cramér-von Mises (CVM) goodness-of-fit test for normality. The p-value is approximated via simulation under the standard normal distribution.
}

\usage{
cvmtest(x, nsim = 10000, seed = NULL)
}

\arguments{
  \item{x}{
  A numeric vector of observations to be tested for normality.
  }

  \item{nsim}{
  Number of Monte Carlo simulations used to approximate the null distribution of the CVM statistic.
   Default is \code{1e4}.
  }

  \item{seed}{
  Optional integer seed for reproducibility. If \code{NULL}, the random  number generator state is not modified.
  }
}

\details{
The Cramér-von Mises statistic measures the squared distance between the empirical distribution function of the data and the theoretical cumulative distribution function under the null hypothesis.

In this implementation, the null hypothesis assumes that the data follow a standard normal distribution. The CVM statistic is computed as:
\deqn{
W^2 = \frac{1}{12n} + \sum_{i=1}^n \left( F(x_{(i)}) - \frac{2i - 1}{2n} \right)^2
}
where \eqn{F} denotes the standard normal cumulative distribution function
and \eqn{x_{(i)}} are the ordered observations.

Because the exact distribution of the statistic is not used, the p-value is estimated via Monte Carlo simulation by repeatedly generating samples from the standard normal distribution and recomputing the CVM statistic.
}

\value{
A list with the following components:
\item{statistic}{Observed Cramér-von Mises test statistic.}
\item{p.value}{Monte Carlo estimated p-value.}
}

\references{
Cramér, H. (1928). On the composition of elementary errors: First paper: Mathematical deductions. 
\emph{Scandinavian Actuarial Journal, 1928}(1), 13-74. \doi{10.1080/03461238.1928.10416862}.

von Mises, R. (1931). Wahrscheinlichkeitsrechnung und ihre Anwendung in der Statistik und 
theoretischen Physik.
}

\author{
Zeynel Cebeci
}

\examples{
set.seed(123)
x <- rnorm(50)

result <- cvmtest(x, nsim = 1000) # Increase in real apps
result$statistic
result$p.value
}

\keyword{normality test}


