\name{mleAL}
\alias{mleAL}

\title{
Maximum Likelihood Estimation of Asymmetric Laplace Distribution}

\description{
This function estimates the parameters of an asymmetric Laplace distribution for a sample.
}

\usage{
mleAL(x)
}

\arguments{
\item{x}{
a numeric vector.
} 
}


\value{

an object of class \code{list} containing the following components:

\item{m}{location parameter}
\item{sigma}{scale parameter}
\item{tau}{skewness parameter}
\item{r}{number of iterations}
}

\references{
Yu K and Zhang J (2005). A three-parameter asymmetric Laplace distribution and its extension. Communications in Statistics-Theory and Methods 34, 1867--1879.
}

\author{
Marco Geraci
}

\seealso{
\code{\link{dal}}, \code{\link{meanAL}}
}

\keyword{asymmetric Laplace distribution}
\keyword{maximum likelihood estimation}

