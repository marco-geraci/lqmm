\name{residuals.lqm}
\alias{residuals.lqm}
\alias{residuals.lqm.counts}

\title{
Residuals from an LQM Objects
}

\description{
This function computes the residuals from a fitted linear quantile model.
}

\usage{
\method{residuals}{lqm}(object, ...)
}

\arguments{
  \item{object}{
an \code{lqm} or \code{lqm.counts} object.
}

  \item{\dots}{
not used.
}
}


\value{
a vector or matrix of residuals.
}


\author{
Marco Geraci
}


\seealso{
\code{\link{lqm}}, \code{\link{lqm.counts}}, \code{\link{predict.lqm}}, \code{\link{coef.lqm}}
}

\keyword{residuals}
