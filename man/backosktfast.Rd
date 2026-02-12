\name{backosktfast}
\alias{backosktfast}
\title{Fast Reverse OSKT Transformation}

\description{
Computes the inverse of the Optimized Skewness and Kurtosis Transformation (OSKT)
using high-performance numerical root-finding algorithms implemented in C++.
The function efficiently recovers original-scale observations from OSKT-transformed
values by solving a nonlinear equation for each observation.
}

\usage{
backosktfast(
  Z, X_mean, X_sd,
  g, h,
  method = "auto",
  tol = 1e-10,
  maxiter_nr = 1000,
  maxiter_brent = 2000
)
}

\arguments{
\item{Z}{
Numeric vector of OSKT-transformed values to be inverted.
Missing values (\code{NA}) are allowed and are propagated to the output.
}

\item{X_mean}{
Numeric scalar. Mean of the original data before standardization.
}

\item{X_sd}{
Numeric scalar. Standard deviation of the original data before standardization.
}

\item{g}{
Numeric scalar. Optimized skewness parameter returned by the OSKT transformation function.
Values close to zero are handled using a numerically stable limiting form.
}

\item{h}{
Numeric scalar. Optimized kurtosis parameter returned by the OSKT transformation function.
Must be non-negative (\eqn{h \ge 0}).
}

\item{method}{
Character string specifying the numerical root-finding strategy:
\describe{
  \item{\code{"auto"}}{(Default) Attempts Newton–Raphson first and falls back to
  Brent–Dekker if convergence fails. Recommended for most use cases.}
  \item{\code{"nr"}}{Pure Newton–Raphson method only. Fastest but no fallback.}
  \item{\code{"brent"}}{Pure Brent–Dekker method only. Most robust but slower.}
}
}

\item{tol}{
Positive numeric scalar specifying the convergence tolerance for the root-finding
algorithms.
}

\item{maxiter_nr}{
Positive integer. Maximum number of iterations allowed for the Newton–Raphson phase.
}

\item{maxiter_brent}{
Positive integer. Maximum number of iterations allowed for the Brent–Dekker phase.
}
}

\details{
The Optimized Skewness and Kurtosis Transformation (OSKT) is defined as
\deqn{
T_{g,h}(x_s) =
\frac{e^{g x_s} - 1}{g} \, e^{\frac{1}{2} h x_s^2},
}
where \eqn{x_s = (X - \mu)/\sigma} is the standardized variable.

When \eqn{g = 0}, the transformation is defined by the continuous limit
\deqn{
T_{0,h}(x_s) = x_s \, e^{\frac{1}{2} h x_s^2}.
}

This function numerically solves the nonlinear equation
\deqn{
T_{g,h}(x_s) = Z
}
for \eqn{x_s}, and then applies the inverse standardization
\deqn{
X = x_s \sigma + \mu.
}

\strong{Numerical Methods:}
\itemize{
  \item \strong{Newton–Raphson} uses analytic derivatives and exhibits quadratic
  convergence near the solution but may fail for extreme values or poor initial guesses.
  \item \strong{Brent–Dekker} is a robust bracketing algorithm combining bisection,
  secant, and inverse quadratic interpolation (Brent, 1973). Convergence is guaranteed 
  if a root is bracketed.
  \item \strong{Auto mode} combines both approaches, achieving high performance
  while retaining robustness.
}

All heavy numerical computations are implemented in C++ via \pkg{Rcpp}.
}

\value{
A list with the following components:
\describe{
  \item{X_orig}{
Numeric vector of inverse-transformed values on the original scale.
Entries are \code{NA} where inversion failed or input values were \code{NA}.
  }
  \item{method_used}{
Character vector of the same length as \code{Z}, indicating which method
succeeded for each observation:
\describe{
  \item{\code{"failed"}}{Root-finding failed to converge.}
  \item{\code{"nr"}}{Newton–Raphson succeeded (when \code{method = "nr"}).}
  \item{\code{"brent"}}{Brent–Dekker succeeded (when \code{method = "brent"}).}
  \item{\code{"auto-nr"}}{Newton–Raphson succeeded in auto mode.}
  \item{\code{"auto-brent"}}{Brent–Dekker fallback succeeded in auto mode.}
}
}
}
}

\author{
Zeynel Cebeci
}

\references{
Brent, R. P. (1973). \emph{Algorithms for Minimization without Derivatives}.
Prentice-Hall, Englewood Cliffs, NJ.
}

\seealso{
\code{\link{oskt}} for the forward OSKT transformation. \cr
\code{\link{osktfast}} for the fast forward OSKT transformation. \cr
\code{\link{backoskt}} for the forward OSKT transformation using uniroot in R. \cr
\code{\link[stats]{uniroot}} for R's base root-finding routine.
}

\examples{
# Auto mode
set.seed(123)
Z <- c(-50, -10, 0, 10, 50)
res <- backosktfast(Z, X_mean = 0, X_sd = 1, g = 0.5, h = 0.1)
res$X_orig
table(res$method_used)

# Pure Newton–Raphson
res_nr <- backosktfast(Z, 0, 1, g = 0.5, h = 0.1, method = "nr")
table(res_nr$method_used)

# Pure Brent–Dekker
res_br <- backosktfast(Z, 0, 1, g = 0.5, h = 0.1, method = "brent")
table(res_br$method_used)

# Handling missing values
Z_na <- c(-10, 0, 10, NA)
backosktfast(Z_na, 0, 1, g = 0.3, h = 0.05)$X_orig
}

\keyword{normality transformation}
\keyword{back transformation}
\keyword{skewness}
\keyword{kurtosis}
\keyword{numerical root finding}
\keyword{Brent-Dekker algorithm}

