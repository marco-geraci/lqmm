\name{coef.lqm}
\alias{coef.lqm}
\alias{coef.lqm.counts}

\title{
Extract LQM Coefficients
}
\description{
\code{coef} extracts model coefficients from \code{lqm}, \code{lqm.counts} objects.
}
\usage{
\method{coef}{lqm}(object, ...)
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
a vector for single quantiles or a matrix for multiple quantiles. 
}

\author{
Marco Geraci
}

\seealso{
\code{\link{lqm}}
\code{\link{summary.lqm}}
\code{\link{lqm.counts}}
}

\keyword{coefficients}

