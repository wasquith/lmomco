\name{cdfemu}
\alias{cdfemu}
\title{Cumulative Distribution Function of the Eta-Mu Distribution}
\description{
This function computes the cumulative probability or nonexceedance probability of the Eta-Mu (\eqn{\eta:\mu}) distribution given parameters (\eqn{\eta} and \eqn{\mu}) computed by \code{\link{parkmu}}. The cumulative distribution function is complex and numerical integration of the probability density function \code{\link{pdfemu}} is used or the Yacoub (2007) \eqn{Y_\nu(a,b)} integral. The cumulative distribution function in terms of this integral is
\deqn{
F(x) = 1- Y_\nu\biggl( \frac{H}{h},\, x\sqrt{2h\mu} \biggr)\mbox{,}
}
where
\deqn{
Y_\nu(a,b) = \frac{2^{3/2 - \nu}\sqrt{\pi}(1-a^2)^\nu}{a^{\nu - 1/2} \Gamma(\nu)} \int_b^\infty x^{2\nu}\,\mathrm{exp}(-x^2)\,I_{\nu-1/2}(ax^2) \; \mathrm{d}x\mbox{,}
}
where \eqn{I_{\nu}(a)}  is the ``\eqn{\nu}th-order modified Bessel function of the first kind.''

}
\usage{
cdfemu(x, para, paracheck=TRUE, yacoubsintegral=TRUE)
}
\arguments{
  \item{x}{A real value vector.}
  \item{para}{The parameters from \code{\link{paremu}} or  \code{\link{vec2par}}.}
  \item{paracheck}{A logical controlling whether the parameters and checked for validity.}
  \item{yacoubsintegral}{A logical controlling whether the integral by Yacoub (2007) is used instead of numerical integration of \code{\link{pdfemu}}.}
}
\value{
  Nonexceedance probability (\eqn{F}) for \eqn{x}.
}
\references{
Yacoub, M.D., 2007, The kappa-mu distribution and the eta-mu distribution: IEEE Antennas and Propagation Magazine, v. 49, no. 1, pp. 68--81
}
\author{ W.H. Asquith}
\seealso{\code{\link{pdfemu}}, \code{\link{quaemu}}, \code{\link{lmomemu}}, \code{\link{paremu}}}
\examples{
para <- vec2par(c(0.5, 1.4), type="emu")
cdfemu(1.2, para, yacoubsintegral=TRUE)
cdfemu(1.2, para, yacoubsintegral=FALSE)
\dontrun{
delx <- 0.01; x <- seq(0,3, by=delx)
nx <- 20*log10(x)
plot(c(-30,10), 10^c(-3,0), log="y", xaxs="i", yaxs="i",
     xlab="RHO", ylab="cdfemu(RHO)", type="n")
m <- 0.75
mus <- c(0.7425, 0.7125, 0.675, 0.6, 0.5, 0.45)
for(mu in mus) {
   eta <- sqrt((m / (2*mu))^-1 - 1)
   lines(nx, cdfemu(x, vec2par(c(eta, mu), type="emu")))
}
mtext("Yacoub (2007, figure 8)")

# Now add some last boundary lines
mu <- m; eta <- sqrt((m / (2*mu))^-1 - 1)
lines(nx, cdfemu(x, vec2par(c(eta, mu), type="emu")),  col=8, lwd=4)
mu <- m/2; eta <- sqrt((m / (2*mu))^-1 - 1)
lines(nx, cdfemu(x, vec2par(c(eta, mu), type="emu")), col=4, lwd=2, lty=2)


delx <- 0.01; x <- seq(0,3, by=delx)
nx <- 20*log10(x)
m <- 0.75; col <- 4; lty <- 2
plot(c(-30,10), 10^c(-3,0), log="y", xaxs="i", yaxs="i",
     xlab="RHO", ylab="cdfemu(RHO)", type="n")
for(mu in c(m/2,seq(m/2+0.01,m,by=0.01), m-0.001, m)) {
   if(mu > 0.67) { col <- 2; lty <- 1 }
   eta <- sqrt((m / (2*mu))^-1 - 1)
   lines(nx, cdfemu(x, vec2par(c(eta, mu), type="emu")),
         col=col, lwd=.75, lty=lty)
}}
}
\keyword{distribution}
\keyword{cumulative distribution function}
\keyword{Distribution: Eta-Mu}
