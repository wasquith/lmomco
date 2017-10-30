\name{quawei}
\alias{quawei}
\title{Quantile Function of the Weibull Distribution}
\description{
This function computes the quantiles of the Weibull distribution given parameters (\eqn{\zeta}, \eqn{\beta}, and \eqn{\delta}) computed by \code{\link{parwei}}. The quantile function  is
\deqn{x(F) = \beta[- \log(1-F)]^{1/\delta} - \zeta \mbox{,}}
where \eqn{x(F)} is the quantile for nonexceedance probability \eqn{F},
\eqn{\zeta} is a location parameter, \eqn{\beta} is a scale parameter, and \eqn{\delta} is a shape parameter.

The Weibull distribution is a reverse Generalized Extreme Value distribution.  As result, the Generalized Extreme Value algorithms are used for implementation of the Weibull in \pkg{lmomco}. The relationS between the Generalized Extreme Value distribution parameters (\eqn{\xi}, \eqn{\alpha}, \eqn{\kappa}) are
\eqn{\kappa}) is \eqn{\kappa = 1/\delta}, 
\eqn{\alpha = \beta/\delta}, and
\eqn{\xi = \zeta - \beta}.
These relations are taken from Hosking and Wallis (1997).

In \R, the quantile function of the Weibull distribution is \code{qweibull}. Given a Weibull parameter object \code{p}, the \R syntax is \code{qweibull(f, p$para[3], scale=p$para[2]) - p$para[1]}. For the current implementation for this package, the reversed Generalized Extreme Value distribution \code{\link{quagev}} is used and the implementation is \code{-quagev((1-f),para)}.
}
\usage{
quawei(f, para, paracheck=TRUE)
}
\arguments{
  \item{f}{Nonexceedance probability (\eqn{0 \le F \le 1}).}
  \item{para}{The parameters from \code{\link{parwei}} or \code{\link{vec2par}}.}
  \item{paracheck}{A logical controlling whether the parameters are checked for validity. Overriding of this check might be extremely important and needed for use of the quantile function in the context of TL-moments with nonzero trimming.}
}
\value{
  Quantile value for nonexceedance probability \eqn{F}.
}
\references{
Hosking, J.R.M., and Wallis, J.R., 1997, Regional frequency analysis---An
approach based on L-moments: Cambridge University Press.
}
\author{ W.H. Asquith}
\seealso{\code{\link{cdfwei}}, \code{\link{pdfwei}}, \code{\link{lmomwei}}, \code{\link{parwei}}}
\examples{
  # Evaluate Weibull deployed here and within R (qweibull)
  lmr <- lmoms(c(123,34,4,654,37,78))
  WEI <- parwei(lmr)
  Q1  <- quawei(0.5,WEI)
  Q2  <- qweibull(0.5,shape=WEI$para[3],scale=WEI$para[2])-WEI$para[1]
  if(Q1 == Q2) EQUAL <- TRUE

  # The Weibull is a reversed generalized extreme value
  Q <- sort(rlmomco(34,WEI)) # generate Weibull sample
  lm1 <- lmoms(Q)    # regular L-moments
  lm2 <- lmoms(-Q)   # L-moment of negated (reversed) data
  WEI <- parwei(lm1) # parameters of Weibull
  GEV <- pargev(lm2) # parameters of GEV
  F <- nonexceeds()  # Get a vector of nonexceedance probs
  plot(pp(Q),Q)
  lines(F,quawei(F,WEI))
  lines(F,-quagev(1-F,GEV),col=2) # line over laps previous
}
\keyword{distribution}
\keyword{quantile function}
\keyword{Distribution: Weibull}
