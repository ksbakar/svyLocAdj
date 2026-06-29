\name{plot.svyLocAdj}
\alias{plot}

\title{Plot Method for svyLocAdj Objects}

\description{
Produces graphical summaries of a fitted \code{svyLocAdj} model.
Depending on the selected plot type, the method displays either a
spatial map of adjusted survey locations or a forest plot of fixed-effect
parameter estimates and their uncertainty intervals.
}

\usage{
\method{plot}{svyLocAdj}(x, type = c("map", "para"), title = "Forest Plot",
                         drop_intercept = FALSE, ...)
}

\arguments{
\item{x}{
An object of class \code{"svyLocAdj"} returned by
\code{\link{svyLocAdj}}.
}

\item{type}{
Character string specifying the type of plot to produce.
Supported values are:
\describe{
\item{\code{"map"}}{Displays the spatial map stored in the fitted model object.}
\item{\code{"para"}}{Displays a forest plot of fixed-effect parameter estimates and their uncertainty intervals.}
}
}

\item{title}{
Title used for the forest plot when \code{type = "para"}.
}

\item{drop_intercept}{
TRUE or FALSE argument, if TRUE then plot intercept.
}.


\item{\dots}{
Additional arguments passed to or from other methods.
Currently ignored.
}
}

\details{
The plotting method provides two complementary visualizations of a fitted
\code{svyLocAdj} model:

\itemize{
\item \strong{Map plot} (\code{type = "map"}): displays the spatial
output associated with the fitted model, typically showing adjusted
survey cluster locations and related geographic information.

\item \strong{Parameter plot} (\code{type = "para"}): displays a forest
plot of the fixed-effect parameter estimates. Point estimates are shown
together with their corresponding uncertainty intervals on the odds-ratio
scale.
}

The parameter plot is useful for assessing the magnitude and direction
of covariate effects, while the map plot provides a visual summary of
the spatial component of the model.
}

\value{
Called for its side effect of producing a graphical display.

The function returns the generated plot object invisibly.
}

\author{
Your Name
}

\seealso{
\code{\link{svymodel}},
\code{\link[svyLocAdj]{summary}}
}

\examples{

## Not run:

fit <- svymodel(
formula = y ~ age_of_child + precipitation,
data = dat
)

# Spatial map

plot(fit)

# Forest plot of fixed effects

plot(fit, type = "para")

# Forest plot with custom title

plot(fit, type = "para", title = "Adjusted Odds Ratios")

## End(Not run)

}

\keyword{models}
\keyword{methods}
\keyword{hplot}
