\name{labor}
\alias{labor}
\title{Labor Pain Data}
\description{
  The \code{labor} data frame has 358 rows and 4 columns of the
  change in pain over time for several 83 women in labor.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{subject}{
      an ordered factor indicating the subject on which the
      measurement was made.  The levels are labelled \code{1}
      to \code{83}.
    }
    \item{pain}{
      a numeric vector of self--reported pain scores on a 100mm line.
    }
    \item{treatment}{
      a dummy variable with values
      \code{1} for subjects who received a pain medication and
      \code{0} for subjects who received a placebo. 
    }
    \item{time}{
      a numeric vector of times (minutes since randomization) at which \code{pain} was measured.
    }

  }
}

\details{
The labor pain data were reported by Davis (1991) and successively analyzed by Jung (1996) and Geraci and Bottai (2007). The data
set consists of repeated measurements of self--reported amount of pain on N = 83 women in labor, of which 43 were randomly assigned to a pain medication group and 40 to a placebo group. The response was measured every 30 min on a 100--mm line, where 0 means no pain and 100 means extreme pain. A nearly monotone pattern of missing data was found for the response variable and the maximum number of measurements for each woman was six.
}
\references{
Geraci M and Bottai M (2007). Quantile regression for longitudinal data using the asymmetric Laplace distribution. Biostatistics 8(1), 140--154.

Jung S (1996). Quasi--likelihood for median regression models. Journal of the American Statistical Association 91, 251--7.
}

\source{
Davis CS (1991). Semi--parametric and non--parametric methods for the analysis of repeated measurements with applications to clinical trials. Statistics in Medicine 10, 1959--80.
}


\keyword{datasets}
