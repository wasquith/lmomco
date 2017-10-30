\name{cdfaep4}
\alias{cdfaep4}
\title{Cumulative Distribution Function of the 4-Parameter Asymmetric Exponential Power Distribution}
\description{
This function computes the cumulative probability or nonexceedance probability of the 4-parameter Asymmetric Exponential Power distribution given parameters (\eqn{\xi}, \eqn{\alpha}, \eqn{\kappa}, and \eqn{h}) computed by \code{\link{paraep4}}. The cumulative distribution function is
\deqn{F(x) = \frac{\kappa^2}{(1+\kappa^2)} \; \gamma([(\xi - x)/(\alpha\kappa)]^h,\; 1/h)\mbox{,}}
for \eqn{x < \xi} and
\deqn{F(x) = 1 - \frac{1}{(1+\kappa^2)} \; \gamma([\kappa(x - \xi)/\alpha]^h,\; 1/h)\mbox{,} }
for \eqn{x \ge \xi}, where \eqn{F(x)} is the nonexceedance probability for quantile \eqn{x}, \eqn{\xi} is a location parameter, \eqn{\alpha} is a scale parameter, \eqn{\kappa} is a shape parameter, \eqn{h} is another shape parameter, and \eqn{\gamma(Z, s)} is the upper tail of the incomplete gamma function for the two arguments. The upper tail of the incomplete gamma function is \code{pgamma(Z, shape, lower.tail=FALSE)} in \R and mathematically is
\deqn{\gamma(Z, a) = \int_Z^\infty y^{a-1} \exp(-y)\, \mathrm{d}y \, /\, \Gamma(a)\mbox{.}}
}
\usage{
cdfaep4(x, para, paracheck=TRUE)
}
\arguments{
  \item{x}{A real value vector.}
  \item{para}{The parameters from \code{\link{paraep4}} or \code{\link{vec2par}}.}
  \item{paracheck}{A logical controlling whether the parameters and checked for validity.}
}

\value{
  Nonexceedance probability (\eqn{F}) for \eqn{x}.
}
\author{ W.H. Asquith}
\references{
Asquith, W.H., 2014, Parameter estimation for the 4-parameter asymmetric exponential power distribution by the method of L-moments using R: Computational Statistics and Data Analysis, v. 71, pp. 955--970.

Delicado, P., and Goria, M.N., 2008, A small sample comparison of maximum likelihood,
moments and L-moments methods for the asymmetric exponential power distribution:
Computational Statistics and Data Analysis, v. 52, no. 3, pp. 1661--1673.
}
\seealso{\code{\link{pdfaep4}}, \code{\link{quaaep4}}, \code{\link{lmomaep4}}, \code{\link{paraep4}}}
\examples{
x <- -0.1
para <- vec2par(c(0, 100, 0.5, 4), type="aep4")
FF <- cdfaep4(-.1,para)
cat(c("F=",FF,"  and estx=",quaaep4(FF, para),"\n"))
\dontrun{
delx <- .1
x <- seq(-20,20, by=delx);
K <- 1;
PAR <- list(para=c(0,1, K, 0.5), type="aep4");
plot(x,cdfaep4(x, PAR), type="n",ylim=c(0,1), xlim=range(x),
     ylab="NONEXCEEDANCE PROBABILITY");
lines(x,cdfaep4(x,PAR), lwd=4);
lines(quaaep4(cdfaep4(x,PAR),PAR), cdfaep4(x,PAR), col=2)
PAR <- list(para=c(0,1, K, 1), type="aep4");
lines(x,cdfaep4(x, PAR), lty=2, lwd=4);
lines(quaaep4(cdfaep4(x,PAR),PAR), cdfaep4(x,PAR), col=2)
PAR <- list(para=c(0,1, K, 2), type="aep4");
lines(x,cdfaep4(x, PAR), lty=3, lwd=4);
lines(quaaep4(cdfaep4(x,PAR),PAR), cdfaep4(x,PAR), col=2)
PAR <- list(para=c(0,1, K, 4), type="aep4");
lines(x,cdfaep4(x, PAR), lty=4, lwd=4);
lines(quaaep4(cdfaep4(x,PAR),PAR), cdfaep4(x,PAR), col=2)}
}
\keyword{distribution}
\keyword{cumulative distribution function}
\keyword{Distribution: Asymmetric Exponential Power}
