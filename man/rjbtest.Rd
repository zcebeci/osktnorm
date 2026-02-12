\name{rjbtest}
\alias{rjbtest}
\title{Gel-Gastwirth Robust Jarque-Bera Test}

\description{
Performs the Gel-Gastwirth robust version of the Jarque-Bera normality
test using quantile-based measures of skewness and kurtosis. The test is
designed to reduce sensitivity to outliers by avoiding moment-based
estimators.
}

\usage{
rjbtest(x)
}

\arguments{
  \item{x}{
  A numeric vector of observations. Missing values are removed prior to
  computation.
  }
}

\details{
The classical Jarque--Bera test relies on moment-based estimates of
skewness and kurtosis, making it highly sensitive to outliers.

The Gel--Gastwirth robust Jarque--Bera (RJB) test replaces these moments
with robust, quantile-based measures.

Robust skewness is measured using the Bowley skewness:
\deqn{
\hat{\gamma}_1^{(R)} =
\frac{Q_{0.75} + Q_{0.25} - 2 Q_{0.50}}
     {Q_{0.75} - Q_{0.25}},
}
where \eqn{Q_p} denotes the empirical \eqn{p}-quantile.

Robust kurtosis is measured using the Moors kurtosis (excess form):
\deqn{
\hat{\gamma}_2^{(R)} =
\frac{(Q_{0.875} - Q_{0.625}) + (Q_{0.375} - Q_{0.125})}
     {Q_{0.75} - Q_{0.25}} - 3.
}

The robust Jarque--Bera test statistic is defined as
\deqn{
\mathrm{RJB}
=
\frac{n}{6} \left( \hat{\gamma}_1^{(R)} \right)^2
+
\frac{n}{24} \left( \hat{\gamma}_2^{(R)} \right)^2.
}

Under the null hypothesis of normality, the statistic is asymptotically
distributed as a chi-squared distribution with 2 degrees of freedom.
}

\value{
An object of class \code{"htest"} containing the following components:
\item{statistic}{The value of the robust Jarque--Bera test statistic.}
\item{p.value}{The asymptotic p-value computed from the chi-squared
distribution with 2 degrees of freedom.}
\item{method}{A character string describing the test.}
\item{data.name}{A character string giving the name of the data.}
}

\references{
Gel, Y. R. and Gastwirth, J. L. (2008). A robust modification of the Jarque--Bera test of normality.
\emph{Economics Letters}, 99(1), 30-32.
}

\examples{
set.seed(123)
x <- rnorm(100)

rjbtest(x)
}
