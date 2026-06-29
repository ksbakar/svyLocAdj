\name{summary.svyLocAdj}
\alias{summary}
\title{Summary for svyLocAdj Objects}
\description{
Produces a summary of a fitted \code{svyLocAdj} model, including
fixed-effect parameter estimates, variability parameters, and model
fit statistics based on the Watanabe--Akaike Information Criterion (WAIC).
}

\usage{
\method{summary}{svyLocAdj}(object, ...)
}

\arguments{
  \item{object}{
    An object of class \code{"svyLocAdj"} returned by
    \code{\link{svyLocAdj}}.
  }

  \item{\dots}{
    Additional arguments passed to or from other methods.
    Currently ignored.
  }
}

\details{
The summary method displays:

\itemize{
  \item The model family used for fitting the DHS spatial location displacement model.
  \item Estimated fixed parameters.
  \item Estimated variability parameters.
  \item The log pointwise predictive density (lppd).
  \item An approximation to the Watanabe--Akaike Information Criterion (WAIC).
}

The method is intended for a concise inspection of model estimates and
overall model fit.
}

\value{
The function is called for its side effect of printing a summary to the
console. It returns \code{NULL} invisibly.
}

\author{
Your Name
}

\seealso{
\code{\link{svymodel}},
\code{\link[svyLocAdj]{plot}}
}

\examples{
## Not run:
fit <- svymodel(
  formula = y ~ x1 + x2,
  data = dat
)

summary(fit)
## End(Not run)
}

\keyword{models}
\keyword{methods}
