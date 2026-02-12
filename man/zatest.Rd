\name{zatest}
\alias{zatest}
\title{Zhang-Wu ZA Test for Normality}
\description{
Performs the Zhang-Wu ZA test for assessing normality.
The test is based on a weighted empirical distribution function statistic
and uses Monte Carlo simulation to obtain p-values under the null
hypothesis of normality with unknown mean and variance.
}

\usage{
zatest(x, nsim = 10000, eps = 1e-10, ncores = 1, seed = NULL)
}

\arguments{
  \item{x}{
  A numeric vector of observations. Missing and non-finite values are
  removed prior to computation.
  }

  \item{nsim}{
  Number of Monte Carlo simulations used to approximate the null
  distribution of the test statistic.
  Larger values improve accuracy at the cost of increased computation
  time.
  }

  \item{eps}{
  A small positive constant used to truncate probability values away from
  0 and 1 to ensure numerical stability in logarithmic computations.
  }

  \item{ncores}{
  Number of CPU cores to be used for parallel Monte Carlo simulation.
  Must be a positive integer. The default \code{ncores = 1} disables
  parallel computation.
  }

  \item{seed}{
  Optional integer value used to set the random number generator seed
  for reproducibility of the Monte Carlo simulations.
  If \code{NULL}, the seed is not set.
  }
}

\details{
Let \eqn{x_{(1)} \le \cdots \le x_{(n)}} denote the ordered sample.
The data are standardized using the sample mean and standard deviation,
and the standard normal cumulative distribution function is evaluated at
the standardized observations.

The Zhang--Wu ZA test statistic is defined as
\deqn{
ZA = -\sum_{i=1}^n \left[
\frac{\log(F(x_{(i)}))}{n - i + 0.5}
+
\frac{\log(1 - F(x_{(i)}))}{i - 0.5}
\right],
}
where \eqn{F} denotes the standard normal cumulative distribution
function.

Because the null distribution of the statistic depends on estimated
parameters, asymptotic critical values are not available. Instead,
p-values are obtained via Monte Carlo simulation under the null
hypothesis by repeatedly generating samples from a normal distribution,
re-estimating the mean and variance, and recomputing the test statistic.

Parallel computation is supported via the \code{ncores} argument.
}

\value{
An object of class \code{"htest"} containing the following components:
\item{statistic}{The observed value of the ZA test statistic.}
\item{p.value}{Monte Carlo p-value for the test.}
\item{method}{A character string describing the test.}
\item{data.name}{A character string giving the name of the data.}
}

\references{
Zhang, J. and Wu, Y. (2005). Likelihood-ratio tests for normality.
\emph{Computational Statistics & Data Analysis}, 49(3), 709-721.
\doi{10.1016/j.csda.2004.05.034}.
}

\examples{
set.seed(123)

# Normal data
x <- rnorm(50)
resx <- zatest(x, nsim=100)
resx$statistic # Test statistic
resx$p.value 

# Log-normal data (non-normal)
set.seed(123)
y <- rlnorm(50, meanlog = 0, sdlog = 1)

resy <- zatest(y, nsim=100)
resy$statistic # ZA test statistic
resy$p.value 

# Exponential data (non-normal)
set.seed(123)
w <- rexp(50)
resw <- zatest(w, nsim=100, ncores=1)
resw$p.value 

# Parallel execution using multiple CPU cores
# (use with care on shared or limited systems)
\dontrun{
 cores <- parallel::detectCores()-1
 resz <- zatest(x, nsim = 10000, ncores = cores)
 resz$statistic
 resz$p.value 
 }
}

\author{
Zeynel Cebeci
}
