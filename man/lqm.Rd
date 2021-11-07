\name{lqm}
\alias{lqm}

\title{Fitting Linear Quantile Models}

\description{\code{lqm} is used to fit linear quantile models based on the asymmetric Laplace distribution.}

\usage{
lqm(formula, data, subset, na.action, weights = NULL, tau = 0.5,
	contrasts = NULL, control = list(), fit = TRUE)}


\arguments{
  \item{formula}{
an object of class \code{\link{formula}} for fixed effects: a symbolic description of the model to be fitted.}

  \item{data}{
an optional data frame, list or environment (or object coercible by \code{\link{as.data.frame}} to a data frame) containing the variables in the model. If not found in data, the variables are taken from \code{environment(formula)}, typically the environment from which \code{lqm} is called.}

  \item{subset}{
an optional vector specifying a subset of observations to be used in the fitting process.
}
  \item{na.action}{
a function which indicates what should happen when the data contain NAs. The default is set by the \code{na.action} setting of \code{\link{options}}.
}
  \item{weights}{An optional vector of weights to be used in the fitting process.}
  
  \item{tau}{
the quantile(s) to be estimated. This must be a number between 0 and 1, otherwise the execution is stopped. If more than one quantile is specified, rounding off to the 4th decimal must give non--duplicated values of \code{tau}, otherwise the execution is stopped.
}
  \item{contrasts}{
an optional list. See the contrasts.arg of \code{\link{model.matrix.default}}.}

  \item{control}{
list of control parameters of the fitting process. See \code{\link{lqmControl}}.
}
  \item{fit}{
logical flag. If \code{FALSE} the function returns a list of arguments to be passed to \code{lqm.fit.gs}.}
}

\details{
The function computes an estimate on the tau-th quantile function of the response, conditional on the covariates, as specified by the formula argument. The quantile predictor is assumed to be linear. The function maximizes the (log)likelihood of a Laplace regression which is equivalent to the minimization of the weighted sum of absolute residuals (Koenker and Bassett, 1978). The optimization algorithm is based on the gradient of the Laplace log--likelihood (Bottai, Orsini and Geraci, 2013).
}

\value{
\code{lqm} returns an object of \code{\link{class}} \code{lqm}.

The function \code{summary} is used to obtain and print a summary of the results.

An object of class \code{lqm} is a list containing the following components:

\item{theta}{a vector of coefficients. \code{theta} is a named matrix of coefficients when \code{tau} is a vector of values.}
\item{scale}{the scale parameter.}
\item{gradient}{the gradient.}
\item{logLik}{the log--likelihood.}
\item{opt}{details on optimization (see \code{\link{lqm.fit.gs}}).}
\item{call}{the matched call.}
\item{term.labels}{names for theta.}
\item{terms}{the terms object used.}
\item{nobs}{the number of observations.}
\item{edf,dim_theta}{the length of theta.}
\item{rdf}{the number of residual degrees of freedom.}
\item{tau}{the estimated quantile(s).}
\item{x}{the model matrix.}
\item{y}{the model response.}
\item{weights}{the weights used in the fitting process (a vector of 1's if \code{weights} = NULL).}
\item{InitialPar}{starting values for theta.}
\item{control}{list of control parameters used for optimization (see \code{\link{lqmControl}}).}
}

\references{
Bottai M, Orsini N, Geraci M (2015). A Gradient Search Maximization Algorithm for the Asymmetric Laplace Likelihood, Journal of Statistical Computation and Simulation, 85(10), 1919-1925.

Chen C (2007). A finite smoothing algorithm for quantile regression. Journal of Computational and Graphical Statistics, 16(1), 136-164. 

Koenker R and Bassett G (1978). Regression Quantiles. Econometrica 46(1), 33--50.
}

\author{
Marco Geraci}

\note{
Updates/FAQ/news are published here \url{https://marcogeraci.wordpress.com/}. New versions are usually published here \url{https://github.com/marco-geraci/lqmm/} before going on CRAN.
}

\seealso{
\code{\link{summary.lqm}, \link{coef.lqm}, \link{predict.lqm}, \link{residuals.lqm}}
}
\examples{

set.seed(123)
n <- 500
p <- 1:3/4
test <- data.frame(x = runif(n,0,1))
test$y <- 30 + test$x + rnorm(n)
fit.lqm <- lqm(y ~ x, data = test, tau = p,
	control = list(verbose = FALSE, loop_tol_ll = 1e-9), fit = TRUE)
fit.lqm
}

\keyword{quantile regression}
