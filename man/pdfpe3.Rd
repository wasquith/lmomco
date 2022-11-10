\name{pdfpe3}
\alias{pdfpe3}
\title{Probability Density Function of the Pearson Type III Distribution}
\description{
This function computes the probability density of the Pearson Type III distribution given parameters (\eqn{\mu}, \eqn{\sigma}, and \eqn{\gamma})  computed by \code{\link{parpe3}}. These parameters are equal to the product moments (\code{\link{pmoms}}): mean, standard deviation, and skew. The probability density function for \eqn{\gamma \ne 0} is
\deqn{f(x) = \frac{Y^{\alpha -1} \exp({-Y/\beta})}
                  {\beta^\alpha\, \Gamma(\alpha)} \mbox{,}}
where \eqn{f(x)} is the probability density for quantile \eqn{x}, \eqn{\Gamma} is the complete gamma function in \R as \code{gamma},
\eqn{\xi} is a location parameter, \eqn{\beta} is a scale parameter,
\eqn{\alpha} is a shape parameter, and \eqn{Y = x - \xi} for \eqn{\gamma > 0} and \eqn{Y = \xi - x} for \eqn{\gamma < 0}. These three \dQuote{new} parameters are related to the product moments (\eqn{\mu}, mean; \eqn{\sigma}, standard deviation; \eqn{\gamma}, skew) by
\deqn{\alpha = 4/\gamma^2 \mbox{,}}
\deqn{\beta  = \frac{1}{2}\sigma |\gamma| \mbox{,\ and}}
\deqn{\xi = \mu - 2\sigma/\gamma \mbox{.}}
If \eqn{\gamma = 0}, the distribution is symmetrical and simply is the probability density Normal distribution with mean and standard deviation of \eqn{\mu} and \eqn{\sigma}, respectively. Internally, the \eqn{\gamma = 0} condition is implemented by \R function \code{dnorm}. The \pkg{PearsonDS} package supports the Pearson distribution system including the Type III (see \bold{Examples}).
}
\usage{
pdfpe3(x, para)
}
\arguments{
  \item{x}{A real value vector.}
  \item{para}{The parameters from \code{\link{parpe3}} or \code{\link{vec2par}}.}
}
\value{
  Probability density (\eqn{f}) for \eqn{x}.
}
\references{
Hosking, J.R.M., 1990, L-moments---Analysis and estimation of
distributions using linear combinations of order statistics: Journal
of the Royal Statistical Society, Series B, v. 52, pp. 105--124.

Hosking, J.R.M., and Wallis, J.R., 1997, Regional frequency analysis---An
approach based on L-moments: Cambridge University Press.
}
\author{W.H. Asquith}
\seealso{\code{\link{cdfpe3}}, \code{\link{quape3}}, \code{\link{lmompe3}}, \code{\link{parpe3}}}
\examples{
  lmr <- lmoms(c(123,34,4,654,37,78))
  pe3 <- parpe3(lmr)
  x <- quape3(0.5,pe3)
  pdfpe3(x,pe3)
\dontrun{
# Demonstrate Pearson Type III between lmomco and PearsonDS
qlmomco.pearsonIII <- function(f, para) {
   MU    <- para$para[1] # product moment mean
   SIGMA <- para$para[2] # product moment standard deviation
   GAMMA <- para$para[3] # product moment skew
   L <- para$para[1] - 2*SIGMA/GAMMA # location
   S <- (1/2)*SIGMA*abs(GAMMA)       # scale
   A <- 4/GAMMA^2                    # shape
   return(PearsonDS::qpearsonIII(f, A, L, S)) # shape comes first!
}
FF <- nonexceeds(); para <- vec2par(c(6,.4,.7), type="pe3")
plot( FF, qlmomco(FF, para), xlab="Probability", ylab="Quantile", cex=3)
lines(FF, qlmomco.pearsonIII(FF, para), col=2, lwd=3) # }

\dontrun{
# Demonstrate forced Pearson Type III parameter estimation via PearsonDS package
para <- vec2par(c(3, 0.4, 0.6), type="pe3"); X <- rlmomco(105, para)
lmrpar <- lmom2par(lmoms(X), type="pe3")
mpspar <- mps2par(X, type="pe3"); mlepar <- mle2par(X, type="pe3")
PDS <- PearsonDS:::pearsonIIIfitML(X) # force function exporting
if(PDS$convergence != 0) {
  warning("convergence failed"); PDS <- NULL # if null, rerun simulation [new data]
} else {
  # This is a list() mimic of PearsonDS::pearsonFitML()
  PDS   <- list(type=3, shape=PDS$par[1], location=PDS$par[2], scale=PDS$par[3])
  skew  <- sign(PDS$shape) * sqrt(4/PDS$shape)
  stdev <-    2*PDS$scale  / abs(skew); mu <- PDS$location + 2*stdev/skew
  PDS <- vec2par(c(mu,stdev,skew), type="pe3") # lmomco form of parameters
}
print(lmrpar$para); print(mpspar$para); print(mlepar$para); print(PDS$para)
#        mu     sigma     gamma 
# 2.9653380 0.3667651 0.5178592 # L-moments (by lmomco, of course)
# 2.9678021 0.3858198 0.4238529 # MPS by lmomco
# 2.9653357 0.3698575 0.4403525 # MLE by lmomco
# 2.9653379 0.3698609 0.4405195 # MLE by PearsonDS
# So we can see for this simulation that the two MLE approaches are similar.}
}
\keyword{distribution}
\keyword{probability density function}
\keyword{Distribution: Pearson Type III}
