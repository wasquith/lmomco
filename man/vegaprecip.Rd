\name{vegaprecip}
\docType{data}
\alias{vegaprecip}
\title{Annual Maximum Precipitation Data for Vega, Texas}
\description{
Annual maximum precipitation data for Vega, Texas
}
\usage{
data(vegaprecip)
}
\format{
An \R \code{data.frame} with
\describe{
  \item{YEAR}{The calendar year of the annual maxima.}
  \item{DEPTH}{The depth of 7-day annual maxima rainfall in inches.}
}
}
\references{
Asquith, W.H., 1998, Depth-duration frequency of precipitation for
Texas: U.S. Geological Survey Water-Resources Investigations Report
98--4044, 107 p.
}
\examples{
data(vegaprecip)
summary(vegaprecip)
}
\keyword{datasets}
