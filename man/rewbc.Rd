\name{rewbc}
\alias{rewbc}
\title{Box-Cox Transformation Using Reweighted Maximum Likelihood}

\description{
Performs a robust reweighted maximum likelihood estimation of the Box-Cox (RBC) 
transformation parameter for univariate data, following the methodology of
Raymaekers and Rousseeuw (2024). The procedure aims at achieving central
normality by iteratively downweighting outlying observations using
Huber-type weights.
}

\usage{
rewbc(x, lrange = seq(-3, 3, by = 0.01), rwsteps = 2, k = 1.5)
}

\arguments{
  \item{x}{
  A numeric vector of observations. Missing values are removed.
  If non-positive values are present, the data are automatically shifted
  to ensure positivity.
  }

  \item{lrange}{
  A numeric vector specifying the grid of candidate Box-Cox
  transformation parameters \eqn{\lambda} over which the (weighted)
  log-likelihood is maximized.
  }

  \item{rwsteps}{
  An integer specifying the number of iterative reweighting steps used
  in the reweighted maximum likelihood procedure.
  }

  \item{k}{
  The tuning constant for the Huber-weight function applied to the
  standardized transformed data. Smaller values lead to stronger
  downweighting of extreme observations.
  }
}

\details{
The function first computes the classical maximum likelihood estimate
(MLE) of the Box-Cox transformation parameter assuming normality of the
transformed data.

In subsequent iterations, the data are transformed using the current
estimate of \eqn{\lambda}. Robust estimates of location and scale
(median and MAD) are used to compute standardized residuals, from which
Huber-type weights are derived. These weights are then used to
re-maximize the Box-Cox log-likelihood over the specified grid of
\eqn{\lambda} values.

The weighted log-likelihood maximized at each step is given by
\deqn{
\ell(\lambda)
= -\frac{n}{2} \log(\sigma^2)
+ (\lambda - 1) \sum_{i=1}^n \log(x_i),
}
where \eqn{\sigma^2} denotes the (possibly weighted) variance of the
transformed data, and robustness enters through the estimation of
\eqn{\sigma^2} via observation weights.
}

\value{
A list with the following components:
\item{transformed}{The Box-Cox transformed data using the final estimated \eqn{\lambda}.}
\item{lambda}{The estimated Box-Cox transformation parameter.}
\item{weights}{The final Huber-type weights assigned to each observation.}
\item{steps}{The number of reweighting iterations performed.}
}

\references{
Raymaekers, J., & Rousseeuw, P. J. (2024). Transforming variables to central normality.
\emph{Machine Learning}, 113(8), 4953--4975.
}

\examples{
set.seed(123)
x <- rgamma(100, shape = 2)

res <- rewbc(x)
res$lambda

hist(res$transformed,
     main = "Reweighted Box-Cox Transformed Data",
     xlab = "Transformed values")
}
\author{
  Zeynel Cebeci
}