\name{predict.lqmm}
\alias{predint}
\alias{predint.lqmm}
\alias{predict.lqmm}

\title{
Predictions from an \code{lqmm} Object
}

\description{
The predictions at level 0 correspond to predictions based only on the fixed effects estimates. The predictions at level 1 are obtained by adding the best linear predictions of the random effects to the predictions at level 0. See details for interpretation. The function \code{predint} will produce 1-alpha confidence intervals based on bootstrap centiles.
}

\usage{
\method{predict}{lqmm}(object, newdata, level = 0,
	na.action = na.pass, ...)
\method{predint}{lqmm}(object, level = 0, alpha = 0.05,
	R = 50, seed = round(runif(1, 1, 10000)))
}

\arguments{
  \item{object}{
an \code{lqmm} object.
}
  \item{newdata}{
an optional data frame in which to look for variables with which to predict. If omitted, the fitted values are produced.  
}
  \item{level}{
an optional integer vector giving the level of grouping to be used in obtaining the predictions.
}
  \item{na.action}{
function determining what should be done with missing values in \code{newdata}. The default is to predict \code{NA}.
}

  \item{alpha}{
1-\code{alpha} is the confidence level.  
}
  \item{R}{
number of bootstrap replications.
  }
  \item{seed}{
optional random number generator seed.
  }
  \item{\dots}{
not used.
}
}

\details{

As discussed by Geraci and Bottai (2014), integrating over the random effects will give "weighted averages" of the cluster-specific quantile effects. These may be interpreted strictly as population regression quantiles for the median (\code{tau=0.5}) only. Therefore, predictions at the population level (\code{code=0}) should be interpreted analogously.

}

\value{
a vector or a matrix of predictions for \code{predict.lqmm}. A data frame or a list of data frames for \code{predint.lqmm} containing predictions, lower and upper bounds of prediction intervals, and standard errors.
}

\references{
Geraci M and Bottai M (2014). Linear quantile mixed models. Statistics and Computing, 24(3), 461--479.
}

\author{
Marco Geraci
}


\seealso{
\code{\link{lqmm}}, \code{\link{ranef.lqmm}}, \code{\link{coef.lqmm}}
}

\examples{
## Orthodont data
data(Orthodont)

# Random intercept model
fitOi.lqmm <- lqmm(distance ~ age, random = ~ 1, group = Subject,
	tau = c(0.1,0.5,0.9), data = Orthodont)

# Predict (y - Xb)	
predict(fitOi.lqmm, level = 0)

# Predict (y - Xb - Zu)
predict(fitOi.lqmm, level = 1)

# 95\% confidence intervals
predint(fitOi.lqmm, level = 0, alpha = 0.05)

}

\keyword{prediction}
