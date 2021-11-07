\name{lqm.fit.gs}
\alias{lqm.fit.gs}

\title{
Quantile Regression Fitting by Gradient Search
}
\description{
This function controls the arguments to be passed to routines written in C for LQM estimation. The optimization algorithm is based on the gradient of the Laplace log--likelihood (Bottai, Orsini and Geraci, 2013).
}
\usage{
lqm.fit.gs(theta, x, y, weights, tau, control)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theta}{
starting values for the regression coefficients.
}
  \item{x}{
the model matrix.
}
  \item{y}{
the model response.
}
  \item{weights}{
the weights used in the fitting process.
}
  \item{tau}{
the quantile to be estimated.
}
  \item{control}{
list of control parameters used for optimization (see \code{\link{lqmControl}}).
}
}

\details{
See argument \code{fit} in \code{\link{lqm}} for generating a list of arguments to be called by this function.
}


\value{
An object of class \code{list} containing the following components:

\item{theta}{a vector of coefficients.}
\item{scale}{the scale parameter.}
\item{gradient}{the gradient.}
\item{logLik}{the log--likelihood.}
\item{opt}{number of iterations when the estimation algorithm stopped.}.
}

\references{
Bottai M, Orsini N, Geraci M (2014). A Gradient Search Maximization Algorithm for the Asymmetric Laplace Likelihood, Journal of Statistical Computation and Simulation, 85, 1919-1925.
}

\author{
Marco Geraci
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{lqm}}
}
\examples{

set.seed(123)
n <- 500
test <- data.frame(x = runif(n,0,1))
test$y <- 30 + test$x + rnorm(n)
lqm.ls <- lqm(y ~ x, data = test, fit = FALSE)

do.call("lqm.fit.gs", lqm.ls)
}

\keyword{fitting}
