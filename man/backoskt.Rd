\name{backoskt}
\alias{backoskt}
\title{Inverse OSKT Normality Transformation}
\description{
Performs numerical back-transformation (inverse transformation) for the
Optimized Skewness and Kurtosis Transformation (OSKT) based on the Tukey
\eqn{g-h} family (Tukey, 1977). The function recovers the original-scale data from the
OSKT-transformed values using either a Newton--Raphson algorithm or a
bracketing root-finding method.
}

\usage{
backoskt(Z, X_mean, X_sd, g, h, tol = 1e-10, maxiter = 1e6, method = c("ur", "nr"))
}

\arguments{
  \item{Z}{
A numeric vector of OSKT-transformed and standardized observations.
}

  \item{X_mean}{
Sample mean of the original data used during standardization.
}

  \item{X_sd}{
Sample standard deviation of the original data used during standardization.
}

  \item{g}{
Skewness parameter of the Tukey \eqn{g-h} transformation.
}

  \item{h}{
Kurtosis (tail heaviness) parameter of the Tukey \eqn{g-h} transformation.
}

  \item{tol}{
Numerical tolerance for convergence of the root-finding algorithms.
Default is \code{1e-10}.
}

  \item{maxiter}{
Maximum number of iterations allowed for the \code{uniroot} method.
Default is \code{1e6}.
}

  \item{method}{
Character string specifying the numerical inversion method.
Either \code{"nr"} for Newton--Raphson (fast but potentially unstable)
or \code{"ur"} for a bracketing root-finding method (robust but slower).
}
}

\details{
The OSKT transformation is based on the Tukey \eqn{g-h} transformation applied
to standardized data. Since the inverse transformation has no closed-form
solution, numerical methods are required.

Two inversion strategies are provided:
\itemize{
  \item{\code{"nr"}: Newton--Raphson iteration initialized at \eqn{x = z},
        offering fast convergence when the derivative is well-behaved.}
  \item{\code{"ur"}: A safe bracketing method using \code{uniroot},
        ensuring convergence at the expense of computational speed.}
}

After inversion, the results are de-standardized using the supplied mean
and standard deviation to recover the original data scale.
}

\value{
A list with the following components:
\itemize{
  \item{\code{X_orig}: Back-transformed observations on the original scale.}
  \item{\code{X_s}: Back-transformed standardized values.}
  \item{\code{time_seconds}: Total computation time in seconds.}
  \item{\code{method}: The numerical inversion method used.}
}
}

\references{
Tukey, J. W. (1977). \emph{Exploratory Data Analysis}. Addison-Wesley.

Headrick, T. C., & Pant, M. D. (2012). Characterizing Tukey \eqn{g} and \eqn{h} distributions through their moments. 
\emph{Journal of Statistical Distributions and Applications}, 1(1), 1--20.
}

\seealso{
\code{\link{oskt}}, \code{\link{osktfast}}, \code{\link{backosktfast}}
}

\examples{
set.seed(123)
x <- rt(200, df = 3)

# Example parameters (typically estimated via oskt)
g <- 0.2
h <- 0.15

# Standardize and apply forward g-h transformation
x_s <- scale(x)
z <- ((exp(g * x_s) - 1) / g) * exp(0.5 * h * x_s^2)

# Back-transformation
res <- backoskt(
  Z = z,
  X_mean = mean(x),
  X_sd = sd(x),
  g = g,
  h = h,
  method = "nr"
)

head(res$X_orig)
}

\author{
  Zeynel Cebeci
}

\keyword{back transformation}
