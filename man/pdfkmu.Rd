\name{pdfkmu}
\alias{pdfkmu}
\title{Probability Density Function of the Kappa-Mu Distribution}
\description{
This function computes the probability density
of the Kappa-Mu (\eqn{\kappa:\mu}) distribution given parameters (\eqn{\kappa} and \eqn{\mu})  computed by \code{\link{parkmu}}. The probability density function  is
\deqn{
f(x) = \frac{2\mu(1+\kappa)^{(\mu + 1)/2}}{\kappa^{(\mu-1)/2}\mathrm{exp}(\mu\kappa)}\,x^\mu\,\exp(-\mu(1+\kappa)x^2)\,I_{\mu - 1}(2\mu\sqrt{\kappa(1+\kappa)}x)\mbox{,}
}
where \eqn{f(x)} is the nonexceedance probability for quantile \eqn{x}, and the modified Bessel function of the first kind is \eqn{I_k(x)}, and define \eqn{m} as
\deqn{m = \frac{\mu(1+\kappa)^2}{1+2\kappa}\mbox{.}}
and for a given \eqn{m}, the new parameter \eqn{\mu} must lie in the range
\deqn{0 \le \mu \le m\mbox{.}}
The definition of \eqn{I_k(x)} is seen under \code{\link{pdfemu}}. Lastly, if \eqn{\kappa = \infty}, then there is a Dirca Delta function of probability at \eqn{x=0}.
}
\usage{
pdfkmu(x, para, paracheck=TRUE)
}
\arguments{
  \item{x}{A real value vector.}
  \item{para}{The parameters from \code{\link{parkmu}} or \code{\link{vec2par}}.}
  \item{paracheck}{A logical controlling whether the parameters and checked for validity.}
}
\value{
  Probability density (\eqn{f}) for \eqn{x}.
}
\references{
Yacoub, M.D., 2007, The kappa-mu distribution and the eta-mu distribution: IEEE Antennas and Propagation Magazine, v. 49, no. 1, pp. 68--81
}
\author{ W.H. Asquith}
\seealso{\code{\link{cdfkmu}}, \code{\link{quakmu}}, \code{\link{lmomkmu}}, \code{\link{parkmu}}}
\examples{
\dontrun{
x <- seq(0,4, by=.1)
para <- vec2par(c(.5, 1.4), type="kmu")
F <- cdfkmu(x, para)
X <- quakmu(F, para, quahi=pi)
plot(F, X, type="l", lwd=8)
lines(F, x, col=2)
}
\dontrun{
# Note that in this example very delicate steps are taken to show
# how one interacts with the Dirac Delta function (x=0) when the m
# is known but mu == 0. For x=0, the fraction of total probability
# is recorded, but when one is doing numerical summation to evaluate
# whether the total probability under the PDF is unity some algebraic
# manipulations are needed as shown for the conditional when kappa
# is infinity.

delx <- 0.001
x <- seq(0,3, by=delx)

plot(c(0,3), c(0,1), xlab="RHO", ylab="pdfkmu(RHO)", type="n")
m <- 1.25
mus <- c(0.25, 0.50, 0.75, 1, 1.25, 0)
for(mu in mus) {
   kappa <- m/mu - 1 + sqrt((m/mu)*((m/mu)-1))
   para <- vec2par(c(kappa, mu), type="kmu")
   if(! is.finite(kappa)) {
      para <- vec2par(c(Inf,m), type="kmu")
      density <- pdfkmu(x, para)
      lines(x, density, col=2, lwd=3)
      dirac <- 1/delx - sum(density[x != 0])
      cumulant <- (sum(density) + density[1]*(1/delx - 1))*delx
      density[x == 0] <- rep(dirac, length(density[x == 0]))
      message("Total integrated probability is ", cumulant, "\n")
   }
   lines(x, pdfkmu(x, para))
}
mtext("Yacoub (2007, figure 3)")
}
}
\keyword{distribution}
\keyword{probability density function}
\keyword{Distribution: Kappa-Mu}
