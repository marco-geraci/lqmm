\name{residuals.lqmm}
\alias{residuals.lqmm}

\title{
Residuals from an \code{lqmm} Object
}

\description{
The residuals at level 0 correspond to population residuals (based only on the fixed effects estimates). The residuals at level 1 are obtained by adding the best linear predictions of the random effects to the predictions at level 0 and the subtracting these from the model response.
}

\usage{
\method{residuals}{lqmm}(object, level = 0, ...)
}


\arguments{
  \item{object}{
an \code{lqmm} object.
}
  \item{level}{
an optional integer vector giving the level of grouping to be used in obtaining the predictions. Level zero corresponds to the population residuals.
}
  \item{\dots}{
not used.
}
}

\value{
a matrix of residuals.
}

\references{
Geraci M and Bottai M (2014). Linear quantile mixed models. Statistics and Computing, 24(3), 461--479. doi: 10.1007/s11222-013-9381-9.
}

\author{
Marco Geraci}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{lqmm}}, \code{\link{predict.lqmm}}, \code{\link{coef.lqmm}}, \code{\link{ranef.lqmm}},
}

\keyword{residuals}
