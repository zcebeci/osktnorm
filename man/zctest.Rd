\name{zctest}
\alias{zctest}
\title{Zhang-Wu ZC Test for Normality}

\description{
Performs the Zhang-Wu ZC goodness-of-fit test for assessing normality.
The test is based on a likelihood-ratio type statistic proposed by
Zhang (2002)<10.1111/1467-9868.00337> and further discussed by Zhang 
and Wu (2005)<10.1016/j.csda.2004.05.034>. The p-value is obtained using 
a Monte Carlo procedure with parameter re-estimation.
}

\usage{
zctest(x, nsim = 10000, eps = 1e-10, ncores = 1, seed = NULL)
}

\arguments{
  \item{x}{
  A numeric vector of observations. Missing and non-finite values are
  removed prior to computation.
  }

  \item{nsim}{
  Number of Monte Carlo simulations used to approximate the null
  distribution of the test statistic.
  Larger values yield more accurate p-values at the expense of increased
  computation time.
  }

  \item{eps}{
  A small positive constant used to truncate normal probabilities away
  from 0 and 1 to ensure numerical stability in logarithmic computations.
  }

  \item{ncores}{
  Number of CPU cores to be used for parallel Monte Carlo simulation.
  The default value \code{1} disables parallel execution.
  Values less than 1 or greater than the number of available cores
  result in an error.
  }

  \item{seed}{
  Optional integer value used to set the random number generator seed
  for reproducibility of the Monte Carlo simulations.
  If \code{NULL}, the seed is not set.
  }
}

\details{
Let \eqn{x_1, \dots, x_n} denote the observed data.
The data are standardized using the sample mean and standard deviation,
and the ordered standardized values \eqn{z_{(i)}} are transformed to
normal scores \eqn{p_i = \Phi(z_{(i)})}, where \eqn{\Phi} denotes the
standard normal distribution function.

The ZC test statistic is defined as
\deqn{
\mathrm{ZC} = \sum_{i=1}^n
\left[
\log\left(
\frac{1/p_i - 1}
{(n - 0.5)/(i - 0.75) - 1}
\right)
\right]^2.
}

Because the finite-sample null distribution of the statistic is not
available in closed form, the p-value is computed using a Monte Carlo
procedure. In each simulation, a normal sample of size \eqn{n} is
generated, standardized using its own sample mean and standard
deviation, and the ZC statistic is recomputed.

The Monte Carlo p-value is computed using the unbiased estimator
\deqn{
p = \frac{1 + \sum_{b=1}^{B} I(T_b \ge T_{\text{obs}})}{B + 1},
}
where \eqn{T_{\text{obs}}} is the observed test statistic and
\eqn{T_b} are the simulated statistics.

To ensure numerical stability, probabilities are truncated to lie in
the interval \eqn{(\varepsilon, 1 - \varepsilon)}.
}

\value{
An object of class \code{"htest"} with the following components:
\item{statistic}{The value of the ZC test statistic.}
\item{p.value}{Monte Carlo p-value for the test.}
\item{method}{A character string describing the test.}
\item{data.name}{A character string giving the name of the data.}
}

\references{
Zhang, J. (2002). Powerful goodness-of-fit tests based on the likelihood ratio.
\emph{Journal of the Royal Statistical Society: Series B}, 64, 281-294.
\doi{10.1111/1467-9868.00337}.

Zhang, J. & Wu, Y. (2005). Likelihood-ratio tests for normality.
\emph{Computational Statistics & Data Analysis}, 49, 709-721.
\doi{10.1016/j.csda.2004.05.034}.
}

\examples{
set.seed(123)
x1 <- rnorm(50)
res1 <- zctest(x1, nsim = 100)
res1$statistic  # Test statistic
res1$p.value # P-value

set.seed(123)
x2 <- rlnorm(50, 0, 1)
res2 <- zctest(x2, nsim = 100)
res2$p.value

# Parallel execution using multiple CPU cores
# (use with care on shared or limited systems)
\dontrun{
 set.seed(123) 
 x3 <- rexp(100)
 cores <- parallel::detectCores()-1
 res3 <- zctest(x3, nsim = 10000, ncores = cores)
 res3$p.value
 }
}

\author{
Zeynel Cebeci
}

