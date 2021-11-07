\name{ranef.lqmm}
\alias{ranef}
\alias{ranef.lqmm}

\title{
Extract Random Effects
}
\description{
This function computes random effects for a linear quantile mixed model.
}
\usage{
\method{ranef}{lqmm}(object, ...)
}

\arguments{
  \item{object}{
an object of \code{\link{class}} \code{lqmm}.
}
  \item{\dots}{
not used.
}

}

\details{
The prediction of the random effects is done via estimated best linear prediction (Geraci and Bottai, 2014). The generic function \code{ranef} is imported from the \code{nlme} package (Pinheiro et al, 2014).
}

\value{
a data frame or a list of data frames of predicted random effects.
}

\references{
Geraci M and Bottai M (2014). Linear quantile mixed models. Statistics and Computing, 24(3), 461--479. doi: 10.1007/s11222-013-9381-9.

Pinheiro J, Bates D, DebRoy S, Sarkar D and R Core Team (2014). nlme: Linear and Nonlinear Mixed Effects Models. R package version 3.1-117, \url{https://CRAN.R-project.org/package=nlme}.
}

\author{
Marco Geraci
}

\seealso{
\code{\link{lqmm}}, \code{\link{coef.lqmm}}
}

\keyword{random effects}
\keyword{coefficients}
