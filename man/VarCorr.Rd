\name{VarCorr.lqmm}
\alias{VarCorr}
\alias{VarCorr.lqmm}

\title{
Extract Variance-Covariance Matrix
}
\description{
This function extracts the variance-covariance matrix of the random effects from a fitted \code{lqmm} object.
}
\usage{
\method{VarCorr}{lqmm}(x, sigma = NULL, ...)
}

\arguments{
  \item{x}{
an object of \code{\link{class}} "lqmm".
}

 \item{sigma}{
not used.
}

 \item{...}{
not used.
}

}

\details{
This function returns the variance or the variance-covariance matrix of the random effects. It calls \code{\link{covHandling}} to manage the output of \code{\link{lqmm.fit.gs}} or \code{\link{lqmm.fit.df}}. A post-fitting approximation to the nearest positive (semi)definite matrix (Higham, 2002) is applied if necessary. The generic function \code{VarCorr} is imported from the \code{nlme} package (Pinheiro et al, 2014).
}

\references{
Higham N (2002). Computing the Nearest Correlation Matrix - A Problem from Finance. IMA Journal of Numerical Analysis, 22, 329-343.

Pinheiro J, Bates D, DebRoy S, Sarkar D and R Core Team (2014). nlme: Linear and Nonlinear Mixed Effects Models. R package version 3.1-117, \url{https://CRAN.R-project.org/package=nlme}.
}

\author{
Marco Geraci
}

\seealso{
\code{\link{lqmm}}
\code{\link{coef.lqmm}}
}

\keyword{covariance}
\keyword{coefficients}
