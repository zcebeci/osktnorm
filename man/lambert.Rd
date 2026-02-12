\name{lambert}
\alias{lambert}
\title{Lambert W x F Normalization (s, h, hh types)}

\description{
Applies Lambert W x F transformations to a numeric vector in order to reduce
skewness and/or heavy tails and to approximate a Gaussian distribution.
This function implements the skewness-only (\code{"s"}), heavy-tail-only
(\code{"h"}), and combined skewness + heavy-tail (\code{"hh"}) models described
by Goerg (2011, 2014).
}

\usage{
lambert(
  x,
  type = c("s", "h", "hh"),
  standardize = TRUE,
  tol = 1e-8,
  maxiter = 100
)
}

\arguments{
\item{x}{
A numeric vector to be Gaussianized. Missing values are allowed and are
preserved in the output.
}

\item{type}{
Character string specifying the transformation type:
\describe{
  \item{\code{"s"}}{Skewness-only Lambert W transformation.}
  \item{\code{"h"}}{Heavy-tail-only (Tukey's h) Lambert W transformation.}
  \item{\code{"hh"}}{Combined skewness and heavy-tail Lambert W transformation.}
}
}

\item{standardize}{
Logical. If \code{TRUE}, the transformed values are centered and scaled to
have mean zero and unit variance.
}

\item{tol}{
Numerical tolerance used in the Newton iterations for evaluating the Lambert W
function.
}

\item{maxiter}{
Maximum number of iterations for the Lambert W Newton solver.
}
}

\details{
Let \eqn{X} denote the observed variable and define the standardized variable
\eqn{Z = (X - \mu) / \sigma}, where \eqn{\mu} and \eqn{\sigma} are the sample mean
and standard deviation.

The Lambert W x F framework assumes that \eqn{Z} is generated from a latent
Gaussian variable \eqn{U} via one of the following models.

\strong{Type = "s" (skewness only):}
\deqn{
Z = U \exp(\delta U)
}
which implies
\deqn{
U = \frac{1}{\delta} W(\delta Z),
}
where \eqn{W(\cdot)} is the principal branch of the Lambert W function.
The parameter \eqn{\delta} is estimated by minimizing the absolute skewness of
\eqn{U}.

\strong{Type = "h" (heavy tails only):}
\deqn{
Z = U \exp\left(\frac{h}{2} U^2\right)
}
which implies
\deqn{
U = \mathrm{sign}(Z)\sqrt{\frac{1}{h} W(h Z^2)}.
}
The parameter \eqn{h} controls tail heaviness and is estimated by minimizing
excess kurtosis.

\strong{Type = "hh" (skewness and heavy tails):}
\deqn{
Z = U \exp(\delta U)\exp\left(\frac{h}{2} U^2\right).
}
In this case no closed-form inverse exists and \eqn{U} is recovered by numerical
root finding. Both \eqn{\delta} and \eqn{h} are estimated jointly by minimizing
skewness and excess kurtosis.

This implementation follows the moment-based estimation strategy described in
Goerg (2011, 2014). It is algorithmically consistent with
\code{LambertW::Gaussianize(type = "s"/"h"/"hh")}, but does not rely on that
package.
}

\value{
A list with the following components:
\describe{
  \item{transformed}{The Gaussianized version of \code{x}.}
  \item{delta}{Estimated skewness parameter.}
  \item{h}{Estimated heavy-tail parameter.}
  \item{mean}{Sample mean of \code{x} prior to transformation.}
  \item{sd}{Sample standard deviation of \code{x} prior to transformation.}
  \item{n}{Number of non-missing observations.}
  \item{type}{Transformation type used.}
  \item{standardize}{Logical indicating whether standardization was applied.}
}
}

\note{
The \code{"s"} transformation is one-to-one and numerically stable for most
applications. The \code{"h"} and \code{"hh"} transformations are particularly
useful for extremely heavy-tailed data, but may fail to extrapolate beyond the
range of the observed sample.
}

\references{
Goerg, G. M. (2011).
Lambert W random variables: A new family of generalized skewed distributions
with applications to risk estimation.
\emph{Annals of Applied Statistics}, 5(3), 2197--2230.

Goerg, G. M. (2014).
The Lambert Way to Gaussianize heavy-tailed data with the inverse of Tukey's h
transformation as a special case.
\emph{The Scientific World Journal}.

Goerg, G. M. (2016).
LambertW: An R package for Lambert W x F random variables.
R package version 0.6.4.
}

\examples{
set.seed(1)

x <- rgamma(200, shape = 1)

fit_s  <- lambert(x, type = "s")
fit_h  <- lambert(x, type = "h")
fit_hh <- lambert(x, type = "hh")

hist(fit_s$transformed, main = "Type s")
hist(fit_h$transformed, main = "Type h")
hist(fit_hh$transformed, main = "Type hh")
}
