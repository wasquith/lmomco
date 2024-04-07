\encoding{utf8}
\name{.lmomcohash}
\docType{data}
\alias{.lmomcohash}
\title{Storage of Lookup Tables for the lmomco Package}
\description{
This is a hidden data object contained in the \code{R/sysdata.rda} file of the \pkg{lmomco} package. The system files \code{inst/doc/}\code{SysDataBuilder01.R} and \code{SysDataBuilder02.R} of the package are responsible for the construction of these data with the exception of the Eta-Mu and Kappa-Mu distribution content.
}
\format{
An \R \code{environment} with entries:
\describe{
  \item{AEPkh2lmrTable}{A \code{data.frame} of asymmetric exponential power (4-parameter) relations between its two shape parameters, numerical, and theoretical L-skew and L-kurtosis. The table stems from \code{inst/doc/SysDataBuilder01.R}. (See also \code{\link{paraep4}})}
  \item{EMU_lmompara_byeta}{A \code{data.frame} of pre-computed table of relations between the parameters and L-moments of the Eta-Mu distribution. (See also \code{\link{lmomemu}}, \code{\link{paremu}})}
  \item{KMU_lmompara_bykappa}{A \code{data.frame} of pre-computed table of relations between the parameters and L-moments of the Kappa-Mu distribution. (See also \code{\link{lmomkmu}}, \code{\link{parkmu}})}
  \item{RiceTable}{A \code{data.frame} with coefficient of L-variation, signal to noise ratio, a parameter G, and L-skew and L-kurtosis of the Rice distribution. This is useful for quick parameter estimation. The table stems from \code{inst/doc/SysDataBuilder01.R}. (See also \code{\link{lmomrice}}, \code{\link{parrice}})}
  \item{RiceTable.maxLCV}{Maximum coefficient of L-variation representable (or apparently so) within \R. The value stems from \code{inst/doc/SysDataBuilder01.R}.}
  \item{RiceTable.minLCV}{Minimum coefficient of L-variation representable (or apparently so) within \R. The value stems from \code{inst/doc/SysDataBuilder01.R}.}
  \item{tau46list}{Various relations of Tau4-Tau6 for symmetrical distributions and used to support the access layer provided by \code{\link{lmrdia46}} for Tau4-Tau6 L-moment ratio diagrams. The tables in the list stem from \code{inst/doc/SysDataBuilder02.R}, which is designed to be run after the \code{SysDataBuilder01.R}.}
}
}
\keyword{datasets}
