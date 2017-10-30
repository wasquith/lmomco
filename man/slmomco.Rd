\name{slmomco}
\alias{slmomco}
\title{Reversed Cumulative Distribution Function (Survival Function) of the Distributions}
\description{
This function acts as an alternative front end to \code{\link{par2cdf}} but reverses the probability to form the survival function. Conceptually, \eqn{S(F) = 1 - F(x)} where \eqn{F(x)} is \code{\link{plmomco}} (implemented by \code{\link{par2cdf}}). The nomenclature of the \code{\link{slmomco}} function is to mimic that of built-in \R functions that interface with distributions.
}
\usage{
slmomco(x, para)
}
\arguments{
  \item{x}{A real value.}
  \item{para}{The parameters from \code{\link{lmom2par}} or similar.}
}
\value{
  Exceedance probability (\eqn{0 \le S \le 1}) for \code{x}.
}
\author{ W.H. Asquith}
\seealso{\code{\link{dlmomco}}, \code{\link{plmomco}}, \code{\link{qlmomco}}, \code{\link{rlmomco}}, \code{\link{add.lmomco.axis}} }
\examples{
para <- vec2par(c(0,1),type='nor') # Standard Normal parameters
exceed <- slmomco(1, para) # percentile of one standard deviation
}
\keyword{reversed cumulative distribution function}
\keyword{survival function}
\concept{survivial distribution}
\concept{survival function}
\keyword{The lmomco functions}
\keyword{The lmomco function mimics of R nomenclature}
