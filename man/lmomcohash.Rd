\encoding{utf8}
\name{.lmomcohash}
\docType{data}
\alias{.lmomcohash}
\title{Storage of Lookup Tables for the lmomco Package}
\description{
This is a hidden data object contained in the \code{R/sysdata.rda} file of the \pkg{lmomco} package. The system file \code{'inst/doc/SysDataBuilder.R'} of the package is responsible for the construction of these data.
}
\format{
An \R \code{environment} with entries:
\describe{
  \item{AEPkh2lmrTable}{A \code{data.frame} of asymmetric exponential power (4-parameter) relations between its two shape parameters, numerical, and theoretical L-skew and L-kurtosis. (See also \code{\link{paraep4}})}
  \item{EMU_lmompara_byeta}{A \code{data.frame} of pre-computed table of relations between the parameters and L-moments of the Eta-Mu distribution. (See also \code{\link{lmomemu}}, \code{\link{paremu}})}
  \item{KMU_lmompara_bykappa}{A \code{data.frame} of pre-computed table of relations between the parameters and L-moments of the Kappa-Mu distribution. (See also \code{\link{lmomkmu}}, \code{\link{parkmu}})}
  \item{RiceTable}{A \code{data.frame} with coefficient of L-variation, signal to noise ratio, a parameter G, and L-skew and L-kurtosis of the Rice distribution. This is useful for quick parameter estimation. (See also \code{\link{lmomrice}}, \code{\link{parrice}})}
  \item{RiceTable.maxLCV}{Maximum coefficient of L-variation representable (or apparently so) within \R.}
  \item{RiceTable.minLCV}{Minimum coefficient of L-variation representable (or apparently so) within \R.}
}
}
\keyword{datasets}
