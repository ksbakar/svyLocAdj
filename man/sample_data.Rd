\name{sample_data}
\docType{data}
\alias{sample_data}
\title{Sample survey data with displaced locations.}
\description{
       This data set contains values of daily
	 8-hour maximum average ozone concentrations (parts per billion (ppb)),
	 maximum  temperature (in degree Celsius),
	 wind speed (knots), and
	 relative humidity,
	 obtained from 28 monitoring sites of New York, USA.

NYgrid: This dataset contains total 6200 rows for 62 days of observations for 10x10 = 100 grid points.

}
\usage{sample_data}
\format{Columns: each contains 435 observations.
	\itemize{
	\item 1st col = DHSCLUST,
	\item 2nd col = maternal_age,
	\item 3rd col = maternal_education,
	\item	4th col = wealth_index,
	\item	5th col = breastfeeding,
	\item	6th col = sex_of_child,
	\item	7th col = age_of_child,
	\item	8th col = hw70,
	\item	9th col = URBAN_RURA.
	\item	10th col = LONGNUM.
	\item	11th col = LATNUM.
	\item	12th col = Precipitation_2020.
	\item	13th col = Mean_Temperature_2020.
	\item	14th col = ch_stunt.
  }
}
\seealso{
\code{\link{svymodel}},
\code{\link{summary.svyLocAdj}},
\code{\link{plot.svyLocAdj}}
}
\examples{
\dontrun{
##

head(sample_data)

##
}
}
\keyword{datasets}
