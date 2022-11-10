\name{pdfkap}
\alias{pdfkap}
\title{Probability Density Function of the Kappa Distribution}
\description{
This function computes the probability density of the Kappa distribution given parameters (\eqn{\xi}, \eqn{\alpha}, \eqn{\kappa}, and \eqn{h})  computed by \code{\link{parkap}}. The probability density function is
\deqn{f(x) = \alpha^{-1} [1-\kappa(x - \xi)/\alpha]^{1/k-1} \times [F(x)]^{1-h}}
where \eqn{f(x)} is the probability density for quantile \eqn{x}, \eqn{F(x)} is the cumulative distribution function  or nonexceedance probability at \eqn{x}, \eqn{\xi} is a location parameter, \eqn{\alpha} is a scale parameter, and \eqn{\kappa} is a shape parameter.
}
\usage{
pdfkap(x, para)
}
\arguments{
  \item{x}{A real value vector.}
  \item{para}{The parameters from \code{\link{parkap}} or \code{\link{vec2par}}.}
}

\value{
  Probability density (\eqn{f}) for \eqn{x}.
}
\references{
Hosking, J.R.M. and Wallis, J.R., 1997, Regional frequency analysis---An
approach based on L-moments: Cambridge University Press.

Sourced from written communication with Dr. Hosking in October 2007.
}
\author{W.H. Asquith}
\seealso{\code{\link{cdfkap}}, \code{\link{quakap}}, \code{\link{lmomkap}}, \code{\link{parkap}}}
\examples{
kap <- vec2par(c(1000,15000,0.5,-0.4),type='kap')
F <- nonexceeds()
x <- quakap(F,kap)
check.pdf(pdfkap,kap,plot=TRUE)
}
\keyword{distribution}
\keyword{probability density function}
\keyword{Distribution: Kappa}
