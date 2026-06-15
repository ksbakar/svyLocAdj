\name{runDHS}
\alias{runDHS}
\title{Fit a Bayesian DHS Location Displacement Adjustment Model}

\description{
Fits a Bayesian model that accounts for the spatial displacement of
Demographic and Health Survey (DHS) cluster locations and obtains
Maximum A Posteriori (MAP) estimates of model parameters.
The function supports Gaussian, binomial, and Poisson response models.
}

\usage{
runDHS(formula, data, family = "binomial", id_var, coord_var,
       cluster_var, urban_rural_var, grid_x, grid_y, sigma_u_prior,
       sigma_r_prior, sigma_nu_prior, beta_prior, zeta_prior,
       digits=4, draws = 10)
}

\arguments{

\item{formula}{
A model formula specifying the response variable and covariates.
}

\item{data}{
A data frame containing the response variable, covariates, cluster
identifiers, coordinates, and other variables required by the model.
}

\item{family}{
The response distribution. Supported values are
\code{"gaussian"}, \code{"binomial"}, and \code{"poisson"}.
}

\item{id_var}{
Character string specifying the cluster identifier variable in
\code{data}.
}

\item{coord_var}{
Character vector of length two specifying the longitude and latitude
variables in \code{data}, e.g., \code{c("lon", "lat")}.
}

\item{cluster_var}{
Character vector specifying cluster-level predictor variables included
in the displacement model.
}

\item{urban_rural_var}{
Character string specifying the variable that identifies urban and rural
DHS clusters.
}

\item{grid_x}{
Number of grid points used along the x-axis when constructing the
spatial grid.
}

\item{grid_y}{
Number of grid points used along the y-axis when constructing the
spatial grid.
}

\item{sigma_u_prior}{
Numeric vector of length two specifying the shape and rate parameters of
the Gamma prior distribution for the urban-area variability parameter,
e.g., \code{c(20, 20)}.
}

\item{sigma_r_prior}{
Numeric vector of length two specifying the shape and rate parameters of
the Gamma prior distribution for the rural-area variability parameter,
e.g., \code{c(10, 25)}.
}

\item{sigma_nu_prior}{
Numeric vector of length two specifying the shape and rate parameters of
the Gamma prior distribution for the \eqn{\nu} variability parameter,
e.g., \code{c(2, 1)}.
}

\item{beta_prior}{
Numeric vector of length two specifying the mean and standard deviation
of the Normal prior distribution for regression coefficients,
e.g., \code{c(0, 2)}.
}

\item{zeta_prior}{
Numeric vector of length two specifying the mean and standard deviation
of the Normal prior distribution for displacement-model coefficients,
e.g., \code{c(0, 2)}.
}

\item{digits}{
Integer specifying the number of decimal places displayed in printed
output.
}

\item{draws}{
Number of posterior draws used when obtaining MAP estimates and related
summaries.
}
}

\details{
The model explicitly accounts for the spatial perturbation applied to
DHS cluster coordinates for confidentiality protection. Spatial grids
are constructed over the study region, and uncertainty arising from
location displacement is incorporated into the estimation procedure.

The returned object is of class \code{"DHSadjust"} and includes methods
for summarization and visualization.
}

\value{
An object of class \code{"DHSadjust"} containing:

\item{results}{
Estimated model parameters, uncertainty measures, and model fit
statistics.
}

\item{data}{
The data used for model fitting.
}

\item{sp_data}{
Spatial information used during estimation, including the generated
spatial grid.
}
}

\references{

Bakar, K. S., et al. (2026).
\emph{A Bayesian model adjusted for location perturbations in sample surveys:
with application to the analysis of demographic and health survey data}.
Journal of the Royal Statistical Society Series A.
doi:10.1093/jrsssa/qnag068

}

\seealso{
\code{\link{summary.DHSadjust}},
\code{\link{plot.DHSadjust}}
}

\examples{
\dontrun{


}
}

\keyword{spatial}
\keyword{bayes}
\keyword{models}
