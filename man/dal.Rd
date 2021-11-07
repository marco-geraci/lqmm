\name{dal}
\alias{pal}
\alias{qal}
\alias{ral}
\alias{dal}

\title{
The Asymmetric Laplace Distribution
}
\description{
Density, distribution function, quantile function and random generation for the asymmetric Laplace distribution.
}
\usage{
dal(x, mu = 0, sigma = 1, tau = 0.5, log = FALSE)
pal(x, mu = 0, sigma = 1, tau = 0.5)
qal(x, mu = 0, sigma = 1, tau = 0.5)
ral(n, mu = 0, sigma = 1, tau = 0.5)
}

\arguments{
\item{x}{
vector of quantiles (\code{dal}, \code{pal}) or probabilities (\code{qal}).
}
  \item{n}{
number of observations.
  }
  \item{mu}{
location parameter.
  }
  \item{sigma}{
positive scale parameter.
  }
  \item{tau}{
skewness parameter (0,1).
 }
 \item{log}{
logical; if \code{TRUE}, probabilities are log--transformed.
 }
 
}

\details{

The asymmetric Laplace distribution with parameters (mu, sigma, tau) has density  
\deqn{f(x) = \tau(1-\tau)/\sigma e^{-1/(2\sigma) (\theta max(x,0) + (1 - \theta) max(-x,0))}}

}

\author{
Marco Geraci
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{lqmm}}, \code{\link{lqm}}
}
\keyword{asymmetric Laplace distribution}
