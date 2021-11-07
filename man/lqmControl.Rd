\name{lqmControl}
\alias{lqmControl}

\title{
Control parameters for lqm estimation
}
\description{
A list of parameters for controlling the fitting process.
}
\usage{
lqmControl(method = "gs1", loop_tol_ll = 1e-5, loop_tol_theta = 1e-3,
	check_theta = FALSE, loop_step = NULL, beta = 0.5, gamma = 1.25,
	reset_step = FALSE, loop_max_iter = 1000, smooth = FALSE,
	omicron = 0.001, verbose = FALSE)
}

\arguments{
  \item{method}{
character vector that specifies which code to use for carrying out the gradient search algorithm: "gs1" (default) based on C code and "gs2" based on R code. Method "gs3" uses a smoothed loss function. See details.
}
  \item{loop_tol_ll}{
tolerance expressed as relative change of the log-likelihood.
}
  \item{loop_tol_theta}{
tolerance expressed as relative change of the estimates.
}
  \item{check_theta}{
logical flag. If \code{TRUE} the algorithm performs a check on the change in the estimates in addition to the likelihood.
}

  \item{loop_step}{
step size (default standard deviation of response).
}

  \item{beta}{
decreasing step factor for line search (0,1).
}

  \item{gamma}{
nondecreasing step factor for line search (>= 1).
}

  \item{reset_step}{
logical flag. If \code{TRUE} the step size is re-setted to the initial value at each iteration.
}

  \item{loop_max_iter}{
maximum number of iterations.
}

  \item{smooth}{
logical flag. If \code{TRUE} the standard loss function is replaced with a smooth approximation.
}

  \item{omicron}{
small constant for smoothing the loss function when using \code{smooth = TRUE}. See details.
}

  \item{verbose}{
logical flag.
}
}

\value{
a list of control parameters.

}

\details{
The methods "gs1" and "gs2" implement the same algorithm (Bottai et al, 2015). The former is based on C code, the latter on R code. While the C code is faster, the R code seems to be more efficient in handling large datasets. For method "gs2", it is possible to replace the classical non-differentiable loss function with a smooth version (Chen, 2007).
}

\references{
Bottai M, Orsini N, Geraci M (2015). A Gradient Search Maximization Algorithm for the Asymmetric Laplace Likelihood, Journal of Statistical Computation and Simulation, 85(10), 1919-1925.

Chen C (2007). A finite smoothing algorithm for quantile regression. Journal of Computational and Graphical Statistics, 16(1), 136-164. 
}

\author{
Marco Geraci}

\seealso{
\code{\link{lqm}}
}

\keyword{fitting}
\keyword{control}
