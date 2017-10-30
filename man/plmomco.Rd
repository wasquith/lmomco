\name{plmomco}
\alias{plmomco}
\title{Cumulative Distribution Function of the Distributions}
\description{
This function acts as an alternative front end to \code{\link{par2cdf}}. The nomenclature of the \code{\link{plmomco}} function is to mimic that of built-in \R functions that interface with distributions.
}
\usage{
plmomco(x, para)
}
\arguments{
  \item{x}{A real value.}
  \item{para}{The parameters from \code{\link{lmom2par}} or similar.}
}
\value{
  Nonexceedance probability (\eqn{0 \le F \le 1}) for \code{x}.
}
\author{ W.H. Asquith}
\seealso{\code{\link{dlmomco}}, \code{\link{qlmomco}}, \code{\link{rlmomco}},  \code{\link{slmomco}}, \code{\link{add.lmomco.axis}}
}
\examples{
para <- vec2par(c(0,1),type='nor') # Standard Normal parameters
nonexceed <- plmomco(1,para) # percentile of one standard deviation
}
\keyword{cumulative distribution function}
\concept{cumulative distribution}
\keyword{The lmomco functions}
\keyword{The lmomco function mimics of R nomenclature}
