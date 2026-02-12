\name{pearsonp}
\alias{pearsonp}
\title{Pearson P Statistic for Normality Check}

\description{
Computes the Pearson P metric for assessing deviation from normality.
The statistic is defined as the Pearson Chi-square goodness-of-fit statistic
divided by its degrees of freedom. This scaled form is used as a
normality metric rather than a formal hypothesis test.
}

\usage{
pearsonp(x, nbins = NULL)
}

\arguments{
  \item{x}{
    A numeric vector of observations.
    Missing values are removed prior to computation.
  }
  \item{nbins}{
    Optional integer specifying the number of equal-probability bins.
    If \code{NULL} (default), the number of bins is set to
    \eqn{\lfloor n^{2/5} \rfloor}, following the default choice in
    \code{nortest::pearson.test}.
  }
}

\details{
The data are standardized using the sample mean and standard deviation.
The standardized values are then grouped into equal-probability bins
defined by the quantiles of the standard normal distribution.

Let \eqn{P} denote the Pearson Chi-square statistic and
\eqn{df = k - 3} the degrees of freedom, accounting for estimation of
the mean and variance. The Pearson P metric is defined as
\deqn{P / df}.

Unlike \code{nortest::pearson.test}, this function does not compute
a p-value and should be interpreted as a descriptive normality metric.
Smaller values indicate closer agreement with normality.
}

\value{
An object of class \code{"htest"} containing:
  \item{statistic}{The Pearson P metric (\eqn{P / df}).}
  \item{method}{A character string describing the metric.}
  \item{data.name}{The name of the input data.}
  \item{df}{Degrees of freedom used in the scaling.}
}

\examples{
set.seed(123)
x <- rlnorm(200)

pearsonp(x)
pearsonp(x, nbins = 8)
}

\author{
Zeynel Cebeci
}