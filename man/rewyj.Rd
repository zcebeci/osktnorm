\name{rewyj}
\alias{rewyj}
\title{Yeo-Johnson Transformation Using Reweighted Maximum Likelihood}
\description{
Performs a robust, reweighted maximum likelihood estimation of the
Yeo-Johnson transformation parameter for univariate data.
Outliers are downweighted using Huber-type weights in an
iteratively reweighted likelihood framework.
}

\usage{
rewyj(x, lrange = seq(-3, 3, by = 0.01), rwsteps = 2, k = 1.5)
}

\arguments{
  \item{x}{
  A numeric vector of observations. Missing values are removed.
  The data may contain both positive and negative values.
  }

  \item{lrange}{
  A numeric vector specifying the grid of candidate
  \eqn{\lambda} values over which the Yeo-Johnson log-likelihood
  is maximized.
  }

  \item{rwsteps}{
  An integer specifying the number of reweighting iterations
  in the iteratively reweighted maximum likelihood procedure.
  }

  \item{k}{
  Tuning constant for the Huber weight function (Huber, 1981). Larger values
  reduce robustness, while smaller values increase downweighting
  of extreme observations.
  }
}

\details{
The function implements the reweighted maximum likelihood (RewML)
approach for the Yeo-Johnson transformation (Yeo &Johnson, 2000) 
as described by Raymaekers and Rousseeuw (2024).

In the first step, the classical maximum likelihood estimate (MLE)
of the Yeo-Johnson transformation parameter \eqn{\lambda} is obtained
under a normality assumption.

Subsequently, the algorithm iteratively:
\enumerate{
  \item Transforms the data using the current estimate of \eqn{\lambda}.
  \item Computes robust location and scale estimates using the median
        and median absolute deviation (MAD).
  \item Standardizes the transformed data and computes Huber-type weights.
  \item Re-maximizes a weighted Yeo--Johnson log-likelihood over the
        specified grid of \eqn{\lambda} values.
}

The Jacobian term of the Yeo-Johnson transformation is included
unweighted, following the formulation in Raymaekers and Rousseeuw (2024).

The weighted log-likelihood has the form
\deqn{
\ell(\lambda) =
-\frac{n}{2} \log(\sigma^2)
+ \sum_{i=1}^n g_\lambda(x_i),
}
where \eqn{\sigma^2} is the weighted variance of the transformed data
and \eqn{g_\lambda(x)} denotes the Jacobian contribution of the
Yeo-Johnson transformation.
}

\value{
A list with the following components:
\item{transformed}{The Yeo-Johnson transformed data using the estimated \eqn{\lambda}.}
\item{lambda}{The estimated Yeo--Johnson transformation parameter.}
\item{weights}{Final robust weights assigned to each observation.}
\item{steps}{Number of reweighting iterations performed.}
}

\references{
Raymaekers, J., and Rousseeuw, P. J. (2024). Transforming variables to central normality.
\emph{Machine Learning}, 113(8), 4953-4975.

Yeo, I.-K. and Johnson, R. A. (2000). A new family of power transformations to improve normality or symmetry. 
\emph{Biometrika}, 87(4), 954-959. \doi{10.1093/biomet/87.4.954}.

Huber, P. J. (1981). \emph{Robust Statistics}. Wiley.
}

\examples{
set.seed(123)
x <- c(rnorm(90), rnorm(10, mean = 5))

res <- rewyj(x)
res$lambda
hist(res$transformed, main = "Reweighted Yeo-Johnson Transformed Data")
}

\author{
Zeynel Cebeci
}
