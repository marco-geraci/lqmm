\name{lqm.counts}
\alias{lqm.counts}
\title{
Quantile Regression for Counts}
\description{
This function is used to fit a quantile regression model when the response is a count variable.
}
\usage{
lqm.counts(formula, data, weights = NULL, offset = NULL, contrasts = NULL,
	tau = 0.5, M = 50, zeta = 1e-05, B = 0.999, cn = NULL, alpha = 0.05,
	control = list())
}

\arguments{
  \item{formula}{
an object of class \code{\link{formula}}: a symbolic description of the model to be fitted.
}
  \item{data}{
an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment(formula), typically the environment from which lqm is called.
}
  \item{weights}{
an optional vector of weights to be used in the fitting process.
}

  \item{offset}{
an optional offset to be included in the model frame.  
}

  \item{contrasts}{
an optional list. See the \code{contrasts.arg} of \code{\link{model.matrix.default}}.
}
  \item{tau}{
quantile to be estimated.
}
  \item{M}{
number of dithered samples.
}
  \item{zeta}{
small constant (see References).
}

  \item{B}{
right boundary for uniform random noise U[0,B] to be added to the response variable (see References).
}

 \item{cn}{
small constant to be passed to \code{\link{F.lqm}} (see References).
}

 \item{alpha}{
significance level.
}
  \item{control}{
list of control parameters of the fitting process. See \code{\link{lqmControl}}.
}

}

\details{A linear quantile regression model if fitted to the log--transformed response. Additional tranformation functions will be implemented. The notation used here follows closely that of Machado and Santos Silva (2005).}

\value{
an object of class "lqm.counts" containing the following components

\item{tau}{the estimated quantile.}
\item{theta}{regression quantile (on the log--scale).}
\item{fitted}{predicted quantile (on the response scale).}
\item{tTable}{coefficients, standard errors, etc.}
\item{x}{the model matrix.}
\item{y}{the model response.}
\item{offset}{offset.}
\item{nobs}{the number of observations.}
\item{M}{specified number of dithered samples for standard error estimation.}
\item{Mn}{actual number of dithered samples used for standard error estimation that gave an invertible D matrix (Machado and Santos Silva, 2005).}
\item{term.labels}{names for theta.}
\item{terms}{the terms object used.}
\item{rdf}{the number of residual degrees of freedom.}
\item{InitialPar}{starting values for theta.}
\item{control}{list of control parameters used for optimization (see \code{\link{lqmControl}}).}
}

\references{
Machado JAF and Santos Silva JMC (2005). Quantiles for counts. Journal of the American Statistical Association, 100(472), 1226--1237. 
}
\author{
Marco Geraci
}



\examples{

n <- 100
x <- runif(n)
test <- data.frame(x = x, y = rpois(n, 2*x))
lqm.counts(y ~ x, data = test, M = 50)


}
\keyword{quantiles for counts}
