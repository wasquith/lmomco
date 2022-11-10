\name{pdfwak}
\alias{pdfwak}
\title{Probability Density Function of the Wakeby Distribution}
\description{
This function computes the probability density
of the Wakeby distribution given parameters (\eqn{\xi}, \eqn{\alpha}, \eqn{\beta}, \eqn{\gamma}, and \eqn{\delta})  computed by \code{\link{parwak}}. The probability density function is
\deqn{f(x) = (\alpha[1-F(x)]^{\beta - 1} + \gamma[1-F(x)]^{-\delta - 1})^{-1}\mbox{,}}
where \eqn{f(x)} is the probability density for quantile \eqn{x}, \eqn{F(x)} is the cumulative distribution function or nonexceedance probability at \eqn{x}, \eqn{\xi} is a location parameter, \eqn{\alpha} and \eqn{\beta} are scale parameters, and \eqn{\gamma}, and \eqn{\delta} are shape parameters. The five returned parameters from \code{\link{parwak}} in order are \eqn{\xi}, \eqn{\alpha}, \eqn{\beta}, \eqn{\gamma}, and \eqn{\delta}.
}
\usage{
pdfwak(x, para)
}
\arguments{
  \item{x}{A real value vector.}
  \item{para}{The parameters from \code{\link{parwak}} or \code{\link{vec2par}}.}
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
\seealso{\code{\link{cdfwak}}, \code{\link{quawak}}, \code{\link{lmomwak}}, \code{\link{parwak}}}
\examples{
\dontrun{
lmr <- vec2lmom(c(1,0.5,.4,.3,.15))
wak <- parwak(lmr)
F <- nonexceeds()
x <- quawak(F,wak)
check.pdf(pdfwak,wak,plot=TRUE)
}
}
\keyword{distribution}
\keyword{probability density function}
\keyword{Distribution: Wakeby}
