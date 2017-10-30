\name{USGSsta14321000peaks}
\docType{data}
\alias{USGSsta14321000peaks}
\title{Annual Peak Streamflow Data for U.S. Geological Survey Streamflow-Gaging Station 14321000}
\description{
Annual peak streamflow data for U.S. Geological Survey streamflow-gaging station 14321000. The peak streamflow-qualification codes \code{Flag} are:
\describe{
\item{1}{Discharge is a Maximum Daily Average}
  \item{2}{Discharge is an Estimate}
  \item{3}{Discharge affected by Dam Failure}
  \item{4}{Discharge less than indicated value, which is Minimum Recordable Discharge at this site}
  \item{5}{Discharge affected to unknown degree by Regulation or Diversion}
  \item{6}{Discharge affected by Regulation or Diversion}
  \item{7}{Discharge is an Historic Peak}
  \item{8}{Discharge actually greater than indicated value}
  \item{9}{Discharge due to Snowmelt, Hurricane, Ice-Jam or Debris Dam breakup}
  \item{A}{Year of occurrence is unknown or not exact}
  \item{B}{Month or Day of occurrence is unknown or not exact}
  \item{C}{All or part of the record affected by Urbanization, Mining, Agricultural changes, Channelization, or other}
  \item{D}{Base Discharge changed during this year}
  \item{E}{Only Annual Maximum Peak available for this year}
}

The gage height qualification codes \code{Flag.1} are:
\describe{
  \item{1}{Gage height affected by backwater}
  \item{2}{Gage height not the maximum for the year}
  \item{3}{Gage height at different site and(or) datum}
  \item{4}{Gage height below minimum recordable elevation}
  \item{5}{Gage height is an estimate}
  \item{6}{Gage datum changed during this year}
}
}
\usage{
data(USGSsta14321000peaks)
}
\format{
An \R \code{data.frame} with
\describe{
  \item{Date}{The date of the annual peak streamflow.}
  \item{Streamflow}{Annual peak streamflow data in cubic feet per second.}
  \item{Flags}{Qualification flags on the streamflow data.}
  \item{Stage}{Annual peak stage (gage height, river height) in feet.}
  \item{Flags.1}{Qualification flags on the gage height data.}
}
}
\examples{
data(USGSsta14321000peaks)
\dontrun{plot(USGSsta14321000peaks)}
}
\keyword{datasets}

