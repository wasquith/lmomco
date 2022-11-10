\name{pdfwei}
\alias{pdfwei}
\title{Probability Density Function of the Weibull Distribution}
\description{
This function computes the probability density of the Weibull distribution given parameters (\eqn{\zeta}, \eqn{\beta}, and \eqn{\delta})  computed by \code{\link{parwei}}. The probability density function is
\deqn{f(x) = \delta Y^{\delta-1} \exp(-Y^\delta)/\beta }
where \eqn{f(x)} is the probability density, \eqn{Y = (x-\zeta)/\beta}, quantile \eqn{x},
\eqn{\zeta} is a location parameter, \eqn{\beta} is a scale parameter, and
\eqn{\delta} is a shape parameter.

The Weibull distribution is a reverse Generalized Extreme Value distribution.  As result, the Generalized Extreme Value algorithms are used for implementation of the Weibull in \pkg{lmomco}. The relations between the Generalized Extreme Value parameters (\eqn{\xi}, \eqn{\alpha}, and \eqn{\kappa}) are \eqn{\kappa = 1/\delta}, \eqn{\alpha = \beta/\delta}, and \eqn{\xi = \zeta - \beta}. These relations are available in Hosking and Wallis (1997).

In \R, the probability distribution function of the Weibull distribution is \code{pweibull}. Given a Weibull parameter object \code{para}, the \R syntax is \code{pweibull(x+para$para[1],} \code{para$para[3],}
\cr \code{scale=para$para[2])}. For the \pkg{lmomco} implmentation, the reversed Generalized Extreme Value distribution \code{\link{pdfgev}} is used and again in \R syntax is \code{pdfgev(-x,para)}.
}
\usage{
pdfwei(x, para)
}
\arguments{
  \item{x}{A real value vector.}
  \item{para}{The parameters from \code{\link{parwei}} or \code{\link{vec2par}}.}
}
\value{
  Probability density (\eqn{f}) for \eqn{x}.
}
\references{
Hosking, J.R.M. and Wallis, J.R., 1997, Regional frequency analysis---An
approach based on L-moments: Cambridge University Press.
}
\author{W.H. Asquith}
\seealso{\code{\link{cdfwei}}, \code{\link{quawei}}, \code{\link{lmomwei}}, \code{\link{parwei}}}
\examples{
  # Evaluate Weibull deployed here and built-in function (pweibull)
  lmr <- lmoms(c(123,34,4,654,37,78))
  WEI <- parwei(lmr)
  F1  <- cdfwei(50,WEI)
  F2  <- pweibull(50+WEI$para[1],shape=WEI$para[3],scale=WEI$para[2])
  if(F1 == F2) EQUAL <- TRUE
\dontrun{
  # The Weibull is a reversed generalized extreme value
  Q <- sort(rlmomco(34,WEI)) # generate Weibull sample
  lm1 <- lmoms( Q)   # regular L-moments
  lm2 <- lmoms(-Q)   # L-moment of negated (reversed) data
  WEI <- parwei(lm1) # parameters of Weibull
  GEV <- pargev(lm2) # parameters of GEV
  F <- nonexceeds()  # Get a vector of nonexceedance probabilities
  plot(pp(Q),Q)
  lines(cdfwei(Q,WEI),Q,lwd=5,col=8)
  lines(1-cdfgev(-Q,GEV),Q,col=2) # line overlaps previous distribution
}
}
\keyword{distribution}
\keyword{probability density function}
\keyword{Distribution: Weibull}
