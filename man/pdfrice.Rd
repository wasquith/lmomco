\name{pdfrice}
\alias{pdfrice}
\title{Probability Density Function of the Rice Distribution}
\description{
This function computes the probability density
of the Rice distribution given parameters  (\eqn{\nu} and \eqn{\mathrm{SNR}})  computed by \code{\link{parrice}}. The probability density function is
\deqn{
f(x) = \frac{x}{\alpha^2}\,\exp\!\left(\frac{-(x^2+\nu^2)}{2\alpha^2}\right)\,I_0(x\nu/\alpha^2)\mbox{,}
}
where \eqn{f(x)} is the nonexceedance probability for quantile \eqn{x},
\eqn{\nu} is a parameter, and \eqn{\nu/\alpha} is a form of signal-to-noise ratio \eqn{\mathrm{SNR}}, and \eqn{I_k(x)} is the modified Bessel function of the first kind, which for integer \eqn{k=0} is seen under \code{\link{LaguerreHalf}}. If \eqn{\nu=0}, then the Rayleigh distribution results and \code{\link{pdfray}} is used. If \eqn{24 < \mathrm{SNR} < 52} is used, then the Normal distribution functions are used with appropriate parameter estimation for \eqn{\mu} and \eqn{\sigma} that include the Laguerre polynomial \code{\link{LaguerreHalf}}. If \eqn{\mathrm{SNR} > 52}, then the Normal distribution functions continue to be used with \eqn{\mu=\alpha\times\mathrm{SNR}} and \eqn{\sigma = A}.
}
\usage{
pdfrice(x, para)
}
\arguments{
  \item{x}{A real value vector.}
  \item{para}{The parameters from \code{\link{parrice}} or \code{\link{vec2par}}.}
}
\value{
  Probability density (\eqn{f}) for \eqn{x}.
}
\note{
The \pkg{VGAM} package provides a pdf of the Rice for reference:
\preformatted{
"drice" <- function(x, vee, sigma, log = FALSE) { # From the VGAM package
    if(!is.logical(log.arg <- log)) stop("bad input for argument 'log'")
    rm(log)
    N = max(length(x), length(vee), length(sigma))
    x = rep(x, len=N); vee = rep(vee, len=N); sigma = rep(sigma, len=N)
    logdensity = rep(log(0), len=N)
    xok = (x > 0)
    x.abs = abs(x[xok]*vee[xok]/sigma[xok]^2)
    logdensity[xok] = log(x[xok]) - 2 * log(sigma[xok]) +
                      (-(x[xok]^2+vee[xok]^2)/(2*sigma[xok]^2)) +
                      log(besselI(x.abs, nu=0, expon.scaled = TRUE)) + x.abs
    logdensity[sigma <= 0] <- NaN; logdensity[vee < 0] <- NaN
    if(log.arg) logdensity else exp(logdensity)
}
}
}
\author{W.H. Asquith}
\references{
Asquith, W.H., 2011, Distributional analysis with L-moment statistics using the R environment for statistical computing: Createspace Independent Publishing Platform, ISBN 978--146350841--8.
}
\seealso{\code{\link{cdfrice}}, \code{\link{quarice}}, \code{\link{lmomrice}}, \code{\link{parrice}}}
\examples{
lmr <- lmoms(c(10, 43, 27, 26, 49, 26, 62, 39, 51, 14))
rice <- parrice(lmr)
x <- quarice(nonexceeds(),rice)
plot(x,pdfrice(x,rice), type="b")


# For SNR=v/a > 24 or 240.001/10 > 24, the Normal distribution is
# used by the Rice as implemented here.
rice1 <- vec2par(c(239.9999,10), type="rice")
rice2 <- vec2par(c(240.0001,10), type="rice")
x <- 200:280
plot( x, pdfrice(x, rice1), type="l", lwd=5, lty=3) # still RICIAN code
lines(x, dnorm(  x, mean=240.0001, sd=10), lwd=3, col=2) # NORMAL obviously
lines(x, pdfrice(x, rice2), lwd=1, col=3) # NORMAL distribution code is triggered

# For SNR=v/a > 52 or 521/10 > 52, the Normal distribution
# used by the Rice as implemented here with simple parameter estimation
# because this high of SNR is beyond limits of Bessel function in Laguerre
# polynomial
rice1 <- vec2par(c(520,10), type="rice")
rice2 <- vec2par(c(521,10), type="rice")
x <- 10^(log10(520) - 0.05):10^(log10(520) + 0.05)
plot( x, pdfrice(x, rice1), type="l", lwd=5, lty=3)
lines(x, pdfrice(x, rice2), lwd=1, col=3) # NORMAL code triggered
}
\keyword{distribution}
\keyword{probability density function}
\keyword{Distribution: Rice}
\keyword{Package: VGAM}
