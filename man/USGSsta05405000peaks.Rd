\name{USGSsta05405000peaks}
\docType{data}
\alias{USGSsta05405000peaks}
\title{Annual Peak Streamflow Data for U.S. Geological Survey Streamflow-Gaging Station 05405000}
\description{
Annual peak streamflow data for U.S. Geological Survey streamflow-gaging station 05405000. The peak streamflow-qualification codes \code{Flag} are:
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
data(USGSsta05405000peaks)
}
\format{
An \R \code{data.frame} with
\describe{
  \item{agency_cd}{Agency code.}
  \item{site_no}{Agency station number.}
  \item{peak_dt}{The date of the annual peak streamflow.}
  \item{peak_tm}{Time of the peak streamflow.}
  \item{peak_va}{Annual peak streamflow data in cubic feet per second.}
  \item{peak_cd}{Qualification flags on the streamflow data.}
  \item{gage_ht}{Annual peak stage (gage height, river height) in feet.}
  \item{gage_ht_cd}{Qualification flags on the gage height data.}
  \item{year_last_pk}{Peak streamflow reported is the highest since this year.}
  \item{ag_dt}{Date of maximum gage-height for water year (if not concurrent with peak).}
  \item{ag_tm}{Time of maximum gage-height for water year (if not concurrent with peak).}
  \item{ag_gage_ht}{Maximum gage height for water year in feet (if not concurrent with peak).}
  \item{ag_gage_ht_cd}{Maximum gage height code.}
}
}
\examples{
data(USGSsta05405000peaks)
\dontrun{plot(USGSsta05405000peaks)}
}
\keyword{datasets}
