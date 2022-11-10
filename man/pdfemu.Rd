\name{pdfemu}
\alias{pdfemu}
\title{Probability Density Function of the Eta-Mu Distribution}
\description{
This function computes the probability density of the Eta-Mu (\eqn{\eta:\mu}) distribution given parameters (\eqn{\eta} and \eqn{\mu})  computed by \code{\link{paremu}}. The probability density function is
\deqn{
f(x) = \frac{4\sqrt{\pi}\mu^{\mu - 1/2}h^\mu}{\gamma(\mu)H^{\mu - 1/2}}\,x^{2\mu}\,\exp(-2\mu h x^2)\,I_{\mu-1/2}(2\mu H x^2)\mbox{,}
}
where \eqn{f(x)} is the nonexceedance probability for quantile \eqn{x}, and the modified Bessel function of the first kind is \eqn{I_k(x)}, and the \eqn{h} and \eqn{H} are
\deqn{
h = \frac{1}{1-\eta^2}\mbox{,}
}
and
\deqn{
H = \frac{\eta}{1-\eta^2}\mbox{,}
}
for \dQuote{Format 2}  as described by Yacoub (2007). This format is exclusively used in the algorithms of the \pkg{lmomco} package.

If \eqn{\mu=1}, then the Rice distribution results, although \code{\link{pdfrice}} is not used. If \eqn{\kappa \rightarrow 0}, then the exact Nakagami-m density function results with a close relation to the Rayleigh distribution.

Define \eqn{m} as
\deqn{m = 2\mu\biggl[1 + {\biggr(\frac{H}{h}\biggl)}^2 \biggr]\mbox{,}}
where for a given \eqn{m}, the parameter \eqn{\mu} must lie in the range
\deqn{m/2 \le \mu \le m\mbox{.}}

The \eqn{I_k(x)} for real \eqn{x > 0} and noninteger \eqn{k} is
\deqn{I_k(x) = \frac{1}{\pi} \int_0^\pi \exp(x\cos(\theta)) \cos(k \theta)\; \mathrm{d}\theta - \frac{\sin(k\pi)}{\pi}\int_0^\infty \exp(-x \mathrm{cosh}(t) - kt)\; \mathrm{d}t\mbox{.}}
}
\usage{
pdfemu(x, para, paracheck=TRUE)
}
\arguments{
  \item{x}{A real value vector.}
  \item{para}{The parameters from \code{\link{paremu}} or \code{\link{vec2par}}.}
  \item{paracheck}{A logical controlling whether the parameters and checked for validity.}
}
\value{
  Probability density (\eqn{f}) for \eqn{x}.
}
\references{
Yacoub, M.D., 2007, The kappa-mu distribution and the eta-mu distribution: IEEE Antennas and Propagation Magazine, v. 49, no. 1, pp. 68--81
}
\author{W.H. Asquith}
\seealso{\code{\link{cdfemu}}, \code{\link{quaemu}}, \code{\link{lmomemu}}, \code{\link{paremu}}}
\examples{
\dontrun{
x <- seq(0,4, by=.1)
para <- vec2par(c(.5, 1.4), type="emu")
F <- cdfemu(x, para);         X <- quaemu(F, para)
plot(F, X, type="l", lwd=8);  lines(F, x, col=2)

delx <- 0.005
x <- seq(0,3, by=delx)
plot(c(0,3), c(0,1), xaxs="i", yaxs="i",
     xlab="RHO", ylab="pdfemu(RHO)", type="n")
mu <- 0.6
# Note that in order to produce the figure correctly using the etas
# shown in the figure that it must be recognized that these are the etas
# for format1, but all of the algorithms in lmomco are built around
# format2
etas.format1 <- c(0, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 1)
etas.format2 <- (1 - etas.format1)/(1+etas.format1)
H <- etas.format2 / (1 - etas.format2^2)
h <-            1 / (1 - etas.format2^2)
for(eta in etas.format2) {
   lines(x, pdfemu(x, vec2par(c(eta, mu), type="emu")),
         col=rgb(eta^2,0,0))
}
mtext("Yacoub (2007, figure 5)")

plot(c(0,3), c(0,2), xaxs="i", yaxs="i",
     xlab="RHO", ylab="pdfemu(RHO)", type="n")
eta.format1 <- 0.5
eta.format2 <- (1 - eta.format1)/(1 + eta.format1)
mus <- c(0.25, 0.3, 0.5, 0.75, 1, 1.5, 2, 3)
for(mu in mus) {
   lines(x, pdfemu(x, vec2par(c(eta, mu), type="emu")))
}
mtext("Yacoub (2007, figure 6)")

plot(c(0,3), c(0,1), xaxs="i", yaxs="i",
     xlab="RHO", ylab="pdfemu(RHO)", type="n")
m <- 0.75
mus <- c(0.7425, 0.75, 0.7125, 0.675, 0.45, 0.5, 0.6)
for(mu in mus) {
   eta <- sqrt((m / (2*mu))^-1 - 1)
   print(eta)
   lines(x, pdfemu(x, vec2par(c(eta, mu), type="emu")))
}
mtext("Yacoub (2007, figure 7)") #}
}
\keyword{distribution}
\keyword{probability density function}
\keyword{Distribution: Eta-Mu}

