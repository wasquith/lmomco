\name{USGSsta06766000dvs}
\docType{data}
\alias{USGSsta06766000dvs}
\title{Daily Mean Streamflow Data for U.S. Geological Survey Streamflow-Gaging Station 06766000}
\description{
Daily mean streamflow data for U.S. Geological Survey streamflow-gaging station 06766000 PLATTE RIVER AT BRADY, NE. The qualification code \code{X01_00060_00003_cd} values are:
\describe{
\item{A}{Approved for publication --- Processing and review completed.}
\item{1}{Daily value is write protected without any remark code to be printed.}
}
}
\usage{
data(USGSsta06766000dvs)
}
\format{
An \R \code{data.frame} with
\describe{
  \item{agency_cd}{The agency code \code{USGS}.}
  \item{site_no}{The station identification number.}
  \item{datetime}{The date and time of the data.}
  \item{X01_00060_00003}{The daily mean streamflow data in cubic feet per second.}
  \item{X01_00060_00003_cd}{A code on the data value.}
}
}
\examples{
data(USGSsta06766000dvs)
\dontrun{plot(USGSsta06766000dvs)}
}
\keyword{datasets}
