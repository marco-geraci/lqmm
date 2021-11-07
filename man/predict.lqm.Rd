\name{predict.lqm}
\alias{predict.lqm}
\alias{predict.lqm.counts}

\title{
Predictions from LQM Objects
}

\description{
This function computes predictions based on fitted linear quantile model.
}

\usage{
\method{predict}{lqm}(object, newdata, interval = FALSE,
	level = 0.95, na.action = na.pass, ...)
\method{predict}{lqm.counts}(object, newdata,
	na.action = na.pass, ...)
}

\arguments{
  \item{object}{
an \code{lqm} or \code{lqm.counts} object.
}
  \item{newdata}{
an optional data frame in which to look for variables with which to predict. If omitted, the fitted values are used.
}
  \item{interval}{
logical flag. If \code{TRUE}, bootstrap percentile intervals for predictions are provided. This argument is for \code{lqm} objects only.
}
  \item{level}{
confidence level. This argument is for \code{lqm} objects only.
}
  \item{na.action}{
function determining what should be done with missing values in \code{newdata}. The default is to predict \code{NA}.}

  \item{\dots}{
further arguments passed to \code{boot.lqm}.
}

}

\value{
a vector or a matrix or an array of predictions.
}


\author{
Marco Geraci
}


\seealso{
\code{\link{residuals.lqm}}, \code{\link{residuals.lqm.counts}}, \code{\link{lqm}}, \code{\link{lqm.counts}}, \code{\link{coef.lqm}}, \code{\link{boot.lqm}}
}

\keyword{prediction}
