\name{summary.lqm}
\alias{summary.lqm}

\title{Summary for an \code{lqm} Object}

\description{
Summary method for class \code{lqm}.}

\usage{
\method{summary}{lqm}(object, method = "boot", alpha = 0.05, covariance = FALSE, ...)
}

\arguments{
  \item{object}{
an object of \code{\link{class}} \code{lqm}}

  \item{method}{
specifies the method used to compute standard errors: "boot" for bootstrap (default), "nid" for large sample approximations under \emph{nid} assumptions.}
  \item{alpha}{
significance level.
}
  \item{covariance}{
logical flag. If \code{TRUE} the covariance matrix of the bootstrap estimates is provided.
}
  \item{\dots}{
see \code{\link{boot.lqm}} for additional arguments.
}
}

\details{
\code{print.summary.lqm} formats the coefficients, standard errors, etc. and additionally gives `significance stars'.}

\value{
an object of class \code{summary.lqm}. The function \code{summary.lqm} computes and returns a list of summary statistics of the fitted linear quantile mixed model given in \code{object}, using the components (list elements) from its argument, plus

\item{Cov}{the covariance matrix obtained from the bootstrapped estimates (if \code{covariance = TRUE}).}
\item{tTable}{a matrix with estimates, standard errors, etc.}
}

\source{
The code for the "nid" method has been adapted from the function \code{summary.rq} in package \code{quantreg}. It depends on the function \code{bandwidth.rq}. 

Roger Koenker (2016). quantreg: Quantile Regression. R package version 5.29. \url{https://CRAN.R-project.org/package=quantreg}

}

\author{
Marco Geraci}

\seealso{
\code{\link{print.summary.lqm}}
\code{\link{lqm}}
}

\examples{

set.seed(12356)
n <- 200
p <- 1:3/4
test <- data.frame(x = runif(n,0,1))
test$y <- 30 + test$x + rnorm(n)
fit.lqm <- lqm(y ~ x, data = test, tau = p)
summary(fit.lqm, R = 50)

}

\keyword{bootstrap}
\keyword{standard errors}
