\name{cdfwei}
\alias{cdfwei}
\title{Cumulative Distribution Function of the Weibull Distribution}
\description{
This function computes the cumulative probability or nonexceedance probability of the Weibull distribution given parameters (\eqn{\zeta}, \eqn{\beta}, and \eqn{\delta}) of the distribution computed by \code{\link{parwei}}. The cumulative distribution function is
\deqn{F(x) = 1 - \exp(Y^\delta) \mbox{,}}
where \eqn{Y} is
\deqn{Y = -\frac{x+\zeta}{\beta}\mbox{,}}
where \eqn{F(x)} is the nonexceedance probability for quantile \eqn{x},
\eqn{\zeta} is a location parameter, \eqn{\beta} is a scale parameter, and
\eqn{\delta} is a shape parameter.

The Weibull distribution is a reverse Generalized Extreme Value distribution.  As result, the Generalized Extreme Value algorithms are used for implementation of the Weibull in this package. The relations between the Generalized Extreme Value parameters (\eqn{\xi}, \eqn{\alpha}, and \eqn{\kappa}) are
\deqn{\kappa = 1/\delta \mbox{,}}
\deqn{\alpha = \beta/\delta \mbox{, and}}
\deqn{\xi = \zeta - \beta \mbox{,}}
which are taken from Hosking and Wallis (1997).

In \R, the cumulative distribution function of the Weibull distribution is \code{pweibull}. Given a Weibull parameter object \code{para}, the \R syntax is \code{pweibull(x+para$para[1], para$para[3],} \cr \code{scale=para$para[2])}. For the current implementation for this package, the reversed Generalized Extreme Value distribution is used \code{1-cdfgev(-x,para)}.
}
\usage{
cdfwei(x, para)
}
\arguments{
  \item{x}{A real value vector.}
  \item{para}{The parameters from \code{\link{parwei}} or  \code{\link{vec2par}}.}
}

\value{
  Nonexceedance probability (\eqn{F}) for \eqn{x}.
}
\references{
Hosking, J.R.M., and Wallis, J.R., 1997, Regional frequency analysis---An
approach based on L-moments: Cambridge University Press.
}
\author{W.H. Asquith}
\seealso{\code{\link{pdfwei}}, \code{\link{quawei}}, \code{\link{lmomwei}}, \code{\link{parwei}}}
\examples{
  # Evaluate Weibull deployed here and within R (pweibull)
  lmr <- lmoms(c(123,34,4,654,37,78))
  WEI <- parwei(lmr)
  F1  <- cdfwei(50,WEI)
  F2  <- pweibull(50+WEI$para[1],shape=WEI$para[3],scale=WEI$para[2])
  if(F1 == F2) EQUAL <- TRUE

  # The Weibull is a reversed generalized extreme value
  Q <- sort(rlmomco(34,WEI)) # generate Weibull sample
  lm1 <- lmoms(Q)    # regular L-moments
  lm2 <- lmoms(-Q)   # L-moment of negated (reversed) data
  WEI <- parwei(lm1) # parameters of Weibull
  GEV <- pargev(lm2) # parameters of GEV
  F <- nonexceeds()  # Get a vector of nonexceedance probs
  plot(pp(Q),Q)
  lines(cdfwei(Q,WEI),Q,lwd=5,col=8)
  lines(1-cdfgev(-Q,GEV),Q,col=2) # line overlaps previous
}
\keyword{distribution}
\keyword{cumulative distribution function}
\keyword{Distribution: Weibull}
