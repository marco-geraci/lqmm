\name{coef.lqmm}
\alias{coef.lqmm}

\title{
Extract LQMM Coefficients
}
\description{
\code{coef} extracts model coefficients from \code{lqmm} objects.
}
\usage{
\method{coef}{lqmm}(object, ...)
}

\arguments{
  \item{object}{
a fitted object of \code{\link{class}} "lqmm".
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
\code{\link{lqmm}}
\code{\link{summary.lqmm}}
}

\keyword{coefficients}
