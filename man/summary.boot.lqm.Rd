\name{summary.boot.lqm}
\alias{summary.boot.lqm}

\title{
Summary for a \code{boot.lqm} Object
}
\description{
Summary method for class \code{boot.lqm}.
}
\usage{
\method{summary}{boot.lqm}(object, alpha = 0.05, digits = max(3, getOption("digits") - 3), ...)
}

\arguments{
  \item{object}{
an object of \code{\link{class}} \code{lqm}.
}

  \item{alpha}{
numeric value for the interval confidence level (\code{1-alpha}).
}

  \item{digits}{
a non-null value for digits specifies the minimum number of significant digits to be printed in values.}

  \item{\dots}{
not used.
}
}


\author{
Marco Geraci
}


\seealso{
\code{\link{boot.lqm}}, \code{\link{lqm}},
}

\keyword{summary}
\keyword{bootstrap}
