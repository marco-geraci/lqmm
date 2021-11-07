\name{lqmm}
\alias{lqmm}

\title{Fitting Linear Quantile Mixed Models}

\description{
\code{lqmm} is used to fit linear quantile mixed models based on the asymmetric Laplace distribution.}

\usage{
lqmm(fixed, random, group, covariance = "pdDiag", tau = 0.5,
	nK = 7, type = "normal", rule = 1, data = sys.frame(sys.parent()),
	subset, weights, na.action = na.fail, control = list(),
	contrasts = NULL, fit = TRUE) 

}
\arguments{
  \item{fixed}{
an object of class \code{\link{formula}} for fixed effects: a symbolic description of the model to be fitted.
}
  \item{random}{
a one-sided formula of the form \code{~x1 + x2 + ... + xn} for random effects: a symbolic description of the model to be fitted.
}
  \item{group}{
grouping factor.
}
  \item{covariance}{
variance--covariance matrix of the random effects. Default is \code{pdDiag} (see details).
}
  \item{tau}{
the quantile(s) to be estimated.
}
  \item{nK}{
number of quadrature knots.
}
  \item{type}{
type of quadrature "c("normal","robust")" (see details).
}
  \item{rule}{
quadrature rule (see details).
}

  \item{data}{
an optional data frame containing the variables named in
   \code{fixed}, \code{random} and \code{group}. By default the variables are taken from the environment from which \code{lqmm} is called.
}
  \item{subset}{
an optional vector specifying a subset of observations to be used in the fitting process.}

  \item{weights}{
an optional vector of weights to be used in the fitting process of the same length as the number of rows of \code{data}. Weights are given to clusters, therefore units within the same cluster receive the same weight (see details).
}
 \item{na.action}{a function that indicates what should happen when the
   data contain \code{NA}s.  The default action (\code{na.fail}) causes
   \code{lqmm} to print an error message and terminate if there are any
   incomplete observations.}

  \item{control}{
  list of control parameters of the fitting process. See \code{\link{lqmmControl}}.
}

  \item{contrasts}{
  not yet implemented.
}

 \item{fit}{
logical flag. If FALSE the function returns a list of arguments to be passed to \code{lqmm.fit}.
}
}


\details{
The function computes an estimate on the tau-th quantile function of the response, conditional on the covariates, as specified by the \code{formula} argument, and on random effects, as specified by the \code{random} argument. The quantile predictor is assumed to be linear. The function maximizes the (log)likelihood of the Laplace regression proposed by Geraci and Bottai (2014). The likelihood is numerically integrated via Gaussian quadrature techniques. The optimization algorithm is based on the gradient of the Laplace log--likelihood (\code{control = list(method = "gs")}). An alternative optimization algorithm is based on a Nelder-Mead algorithm (\code{control = list(method = "df")}) via \code{\link{optim}}. The scale parameter is optimized in a refinement step via \code{\link{optimize}}.

Quadrature approaches include Gauss-Hermite (\code{type = "normal"}) and Gauss-Laguerre (\code{type = "robust"}) quadrature. The argument \code{rule} takes one of the following: 1 (product rule quadrature), 2 (sparse grid quadrature), 3 (nested quadrature rule - only for \code{type = "normal"}), 4 (quadrature rule with the smallest number of nodes between rules 1 or 2). Rules 2 and 3 have not yet been tested extensively.

Different standard types of positive--definite matrices for the random effects can be specified: \code{pdIdent} multiple of an identity; \code{pdCompSymm} compound symmetry structure (constant diagonal and constant off--diagonal elements); \code{pdDiag} diagonal; \code{pdSymm} general positive--definite matrix, with no additional structure.

Weights are given to clusters, therefore it is expected that these are constant within cluster. When the weights are specified in the main call, then the first value by \code{group} in the vector \code{weights} will be replicated for the same length of each group. Alternatively, different weights within the same cluster can be introduced with a direct call to \code{\link{lqmm.fit.gs} or \link{lqmm.fit.df}}.

The \code{lqmm} vignette can be accessed by typing \code{help(package = "lqmm")} and then following the link 'User guides, package vignettes and other documentation'.
}

\value{
\code{lqmm} returns an object of \code{\link{class}} \code{lqmm}.

The function \code{summary} is used to obtain and print a summary of the results.

An object of class \code{lqmm} is a list containing the following components:

\item{theta}{a vector containing fixed regression coefficients and parameters of the variance--covariance matrix of the random effects. See \code{\link{VarCorr.lqmm}} to extract the variance--covariance of the random effects from an "lqmm" object.}
\item{theta_x,theta_z}{partition of \code{theta}: fixed regression coefficients (\code{theta_x}) and unique variance--covariance parameters (\code{theta_z}).}
\item{scale}{the scale parameter.}
\item{gradient}{the gradient (\code{control = list(method = "gs")}).}
\item{logLik}{the log--likelihood.}
\item{opt}{details on optimization (see \code{\link{lqmm.fit.gs}} and \code{\link{lqmm.fit.df}}).}
\item{call}{the matched call.}
\item{nn}{column names of \code{mmf}.}
\item{mm}{column names of \code{mmr}.}
\item{nobs}{the number of observations.}
\item{dim_theta}{the number of columns in \code{mmf} and \code{mmr}.}
\item{dim_theta_z}{the length of \code{theta_z}.}
\item{edf}{length of \code{theta}.}
\item{rdf}{the number of residual degrees of freedom.}
\item{df}{edf + 1 (scale parameter).}
\item{tau}{the estimated quantile(s).}
\item{mmf}{the model matrix -- fixed effects.}
\item{mmr}{the model matrix -- random effects.}
\item{y}{the model response.}
\item{revOrder}{original order of observations (now ordered according to \code{group}).}
\item{weights}{the likelihood weights used in the fitting process (a vector of 1's if \code{weights} is missing or \code{NULL}).}
\item{group}{the grouping factor.}
\item{ngroups}{the number of groups.}
\item{QUAD}{quadrature nodes and weights.}
\item{type}{the type of quadrature.}
\item{rule}{quadrature rule.}
\item{InitialPar}{starting values for theta.}
\item{control}{list of control parameters used for optimization (see \code{\link{lqmmControl}}).}
\item{cov_name}{class of variance-covariance matrix for the random effects.}
\item{mfArgs}{arguments for \code{\link{model.frame}} to return the full data frame.}
}

\references{
Genz A, and Keister BD (1996). Fully symmetric interpolatory rules for multiple integrals over infinite regions with Gaussian weight. Journal of Computational and Applied Mathematics, 71(2), 299--309. <doi:10.1016/0377-0427(95)00232-4>

Geraci M (2014). Linear quantile mixed models: The lqmm package for Laplace quantile regression. Journal of Statistical Software, 57(13), 1--29. <doi:10.18637/jss.v057.i13>

Geraci M and Bottai M (2007). Quantile regression for longitudinal data using the asymmetric Laplace distribution. Biostatistics 8(1), 140--154. <doi:10.1093/biostatistics/kxj039>

Geraci M and Bottai M (2014). Linear quantile mixed models. Statistics and Computing, 24(3), 461--479. <doi:10.1007/s11222-013-9381-9>.

Heiss F, and Winschel V (2008). Likelihood approximation by numerical integration on sparse grids. Journal of Econometrics, 144(1), 62--80. <doi:10.1016/j.jeconom.2007.12.004>
}

\author{
Marco Geraci
}

\note{
Updates/FAQ/news are published here \url{http://marcogeraci.wordpress.com/}. New versions are usually published here \url{https://r-forge.r-project.org/R/?group_id=1396} before going on CRAN.
}

\seealso{
\code{\link{lqm}, \link{summary.lqmm}, \link{coef.lqmm}, \link{VarCorr.lqmm}, \link{predict.lqmm}, \link{residuals.lqmm}}
}

\examples{

# Test example
set.seed(123)

M <- 50
n <- 10
test <- data.frame(x = runif(n*M,0,1), group = rep(1:M,each=n))
test$y <- 10*test$x + rep(rnorm(M, 0, 2), each = n) + rchisq(n*M, 3)
fit.lqmm <- lqmm(fixed = y ~ x, random = ~ 1, group = group,	
	data = test, tau = 0.5, nK = 11, type = "normal")
fit.lqmm

#Call: lqmm(fixed = y ~ x, random = ~1, group = group, tau = 0.5, nK = 11, 
#    type = "normal", data = test)
#Quantile 0.5 

#Fixed effects:
#(Intercept)            x  
#      3.443        9.258  

#Covariance matrix of the random effects:
#(Intercept) 
#      3.426 

#Residual scale parameter: 0.8697 (standard deviation 2.46)
#Log-likelihood: -1178 

#Number of observations: 500 
#Number of groups: 50 


## Orthodont data
data(Orthodont)

# Random intercept model
fitOi.lqmm <- lqmm(distance ~ age, random = ~ 1, group = Subject,
	tau = c(0.1,0.5,0.9), data = Orthodont)
coef(fitOi.lqmm)

# Random slope model
fitOs.lqmm <- lqmm(distance ~ age, random = ~ age, group = Subject,
	tau = c(0.1,0.5,0.9), cov = "pdDiag", data = Orthodont)

# Extract estimates
VarCorr(fitOs.lqmm)
coef(fitOs.lqmm)
ranef(fitOs.lqmm)

# AIC
AIC(fitOi.lqmm)
AIC(fitOs.lqmm)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{quantile regression}
