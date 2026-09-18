\name{summary}
\alias{summary.svyLocAdj}
\title{Summary for svyLocAdj Objects}
\description{
Produces a summary of a fitted \code{svyLocAdj} model, including
fixed-effect parameter estimates, variability parameters, and model
fit statistics based on the Watanabe-Akaike Information Criterion (WAIC).
}
\usage{
\method{summary}{svyLocAdj}(object, ...)
}
\arguments{
\item{object}{An object of class \code{"svyLocAdj"}.}
\item{\dots}{Additional arguments passed to or from other methods. Currently ignored.}
}
\details{
The summary method displays:

\itemize{
  \item The model family used for fitting the DHS spatial location displacement model.
  \item Estimated fixed parameters.
  \item Estimated variability parameters.
  \item The log pointwise predictive density (lppd).
  \item An approximation to the Watanabe-Akaike Information Criterion (WAIC).
}

The method is intended for a concise inspection of model estimates and
overall model fit.
}

\value{
  Summary output.
}

\seealso{
\code{\link{svymodel}},
\code{\link{plot.svyLocAdj}}
}

\examples{
\donttest{

## Not run:

head(sample_data)

fit <- svymodel(
  formula = ch_stunt ~ age_of_child,
  data = sample_data
)

summary(fit)

## End(Not run)
}
}

\keyword{models}
\keyword{methods}
