\name{covHandling}
\alias{covHandling}

\title{
Variance-Covariance Matrix 
}

\description{
This is an auxiliary function. 
}

\usage{
covHandling(theta, n, cov_name, quad_type)
}

\arguments{
  \item{theta}{
unique parameters of the variance-covariance matrix of the random effects as returned by \code{\link{lqmm}} in \code{theta_z}.
}
  \item{n}{
dimension of the vector of random effects.
}
  \item{cov_name}{
see argument \code{covariance} in \code{\link{lqmm}}.
}
  \item{quad_type}{
type of quadrature "c("normal","robust")".
}
}


\author{
Marco Geraci
}


\seealso{
\code{\link{VarCorr.lqmm}}
}
\keyword{covariance}

