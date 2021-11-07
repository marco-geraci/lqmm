\name{lqmmControl}
\alias{lqmmControl}

\title{
Control parameters for lqmm estimation
}
\description{
A list of parameters for controlling the fitting process.
}

\usage{
lqmmControl(method = "gs", LP_tol_ll = 1e-5, LP_tol_theta = 1e-5,
	check_theta = FALSE, LP_step = NULL, beta = 0.5, gamma = 1,
	reset_step = FALSE, LP_max_iter = 500, UP_tol = 1e-4,
	UP_max_iter = 20, startQR = FALSE, verbose = FALSE)
}

\arguments{
  \item{method}{
character vector that specifies the estimation method: "gs" for gradient search (default) and "df" for Nelder-Mead.
}
 
 \item{LP_tol_ll}{
tolerance expressed as absolute change of the log-likelihood.
}

  \item{LP_tol_theta}{
tolerance expressed as absolute change of \code{theta}
}

  \item{check_theta}{
logical flag. If TRUE the algorithm performs an additional check on the change in the estimates.
}

  \item{LP_step}{
step size (default standard deviation of response).
 }

 \item{beta}{
decreasing step factor for line search (0,1).
}

  \item{gamma}{
nondecreasing step factor for line search (>= 1).
}

  \item{reset_step}{
logical flag. If \code{TRUE} the step size is reset to the initial value at each iteration.
}

 \item{LP_max_iter}{
maximum number of iterations
 }
  \item{UP_tol}{
tolerance expressed as absolute change of the \code{scale} parameter.
}
  \item{UP_max_iter}{
maximum number of iterations.
}
  \item{startQR}{
logical flag. If \code{FALSE} (default) the least squares estimate of the fixed effects is used as starting value of \code{theta_x} and \code{scale}. If \code{TRUE} the \code{\link{lqm}} estimate is used.
 }
  \item{verbose}{
logical flag.
}

}

\details{
\code{LP} (lower loop) refers to the estimation of regression coefficients and variance-covariance parameters. \code{UP} (upper loop) refers to the estimation of the scale parameter.
}

\value{
a list of control parameters.
}

\author{
Marco Geraci}

\seealso{
\code{\link{lqmm}}
}
\keyword{fitting}
\keyword{control}
