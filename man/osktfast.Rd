\name{osktfast}
\alias{osktfast}
\title{Optimized Skewness-Kurtosis Transformation (OSKT)}
\description{
Performs a g-and-h transformation to reduce skewness and kurtosis, 
minimizing the Anderson-Darling A2 statistic. The transformation aims to 
normalize non-normal data by targeting skewness and kurtosis simultaneously.
}
\usage{
osktfast(x, 
         init_params = c(0.1, 0.1), 
         lower_bounds = c(-1, 0), 
         upper_bounds = c(1, 0.5), 
         maxiter = 200)
}
\arguments{
  \item{x}{Numeric vector of input data.}
  \item{init_params}{Numeric vector of length 2: initial values for g and h. Default is c(0.1, 0.1).}
  \item{lower_bounds}{Numeric vector of length 2: lower bounds for g and h. Default is c(-1, 0).}
  \item{upper_bounds}{Numeric vector of length 2: upper bounds for g and h. Default is c(1, 0.5).}
  \item{maxiter}{Integer: maximum number of iterations for the optimizer. Default is 200.}
}
\value{
A list containing:
\item{transformed}{The transformed numeric vector.}
\item{g}{Estimated g parameter of the g-and-h transformation.}
\item{h}{Estimated h parameter of the g-and-h transformation.}
\item{value}{Anderson-Darling A2 statistic at the optimum.}
}
\details{
This function uses a pure C++ implementation of the L-BFGS-B algorithm to optimize the
g-and-h transformation parameters for normality. The objective function is the 
modified Anderson-Darling A2 statistic with Stephen's correction for small samples. 
It is suitable for moment-targeting normalization is desired.
}
\examples{
set.seed(123)
x <- rnorm(100, mean=5, sd=2) # Simulate non-normal data
res <- osktfast(x)

# Check transformed data
head(res$transformed)
res$g
res$h
res$value
}
\author{
Zeynel Cebeci
}
\keyword{transform}
\keyword{kurtosis}
\keyword{skewness}
\keyword{optimization}
\keyword{normality}
