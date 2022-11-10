\name{cdfkmu}
\alias{cdfkmu}
\title{Cumulative Distribution Function of the Kappa-Mu Distribution}
\description{
This function computes the cumulative probability or nonexceedance probability of the Kappa-Mu (\eqn{\kappa:\mu})  distribution given parameters (\eqn{\kappa} and \eqn{\mu}) computed by \code{\link{parkmu}}. The cumulative distribution function is complex and numerical integration of the probability density function \code{\link{pdfkmu}} is used. Alternatively, the cumulative distribution function may be defined in terms of the Marcum Q function
\deqn{
F(x) = 1 - Q_\nu\biggl(\sqrt{2\kappa\mu},\, x\sqrt{2(1+\kappa)\mu}\biggr)\mbox{,}
}
where \eqn{F(x)} is the nonexceedance probability for quantile \eqn{x} and \eqn{Q_v(a,b)} is the Marcum Q function defined by
\deqn{
Q_\nu(a,b) = \frac{1}{\alpha^{\nu-1}}\int_b^\infty t^\nu \, \exp(-(t^2 + a^2)/2) \, I_{\nu-1}(at)\; \mathrm{d}t\mbox{,}
}
which can be numerically difficult to work with and particularly so with real number values for \eqn{\nu}. \eqn{I_\nu(a)}  is the ``\eqn{\nu}th-order modified Bessel function of the first kind.''

Following an apparent breakthrough(?) by Shi (2012), \eqn{\nu} can be written as \eqn{\nu = n + \Delta} where \eqn{n} is an integer and \eqn{0 < \Delta \le 1}. The author of \pkg{lmomco} refers to this alternative formulation as the \dQuote{delta nu method}. The Marcum Q function for \eqn{\nu > 0} (\eqn{n = 1,2,3, \cdots)} is
\deqn{
Q_\nu(a,b) = Q_\Delta(a,b) + \exp(-(a^2 + b^2)/2) \, \sum_{i=0}^{n-1}\biggl(\frac{b}{a}\biggr)^{i+\Delta} \, I_{i+\Delta}(ab)\mbox{,}
}
and the function for \eqn{\nu \le 0} (\eqn{n=-1,-2,-3,\cdots}) is
\deqn{
Q_\nu(a,b) = Q_\Delta(a,b) - \mathrm{exp}(-(a^2 + b^2)/2) \times \sum_{i=n}^{-1}\biggl(\frac{b}{a}\biggr)^{i+\Delta} \mathrm{I}_{i+\Delta}(ab)\mbox{,}
}
and the function for \eqn{\nu = 0} is
\deqn{
Q_\nu(a,b) = Q_\Delta(a,b) + \mathrm{exp}(-(a^2 + b^2)/2)\mbox{.}
}

Shi (2012) concludes that the \dQuote{merit} of these two expressions is that the evaulation of the Marcum Q function is reduced to the numerical evaluation of \eqn{Q_\Delta(a,b)}. This difference can result in measurably faster computation times (confirmed by limited analysis by the author of \pkg{lmomco}) and possibly better numerical performance.

Shi (2012) uses notation and text that implies evaluation of the far-right additive term (the summation) for \eqn{n=0} as part of the condition \eqn{\nu > 0}. To clarify, Shi (2012) implies for \eqn{\nu > 0; n = 0} (but \eqn{n=0} occurs also for \eqn{-1 < \nu <= 0}) the following computation
\deqn{
Q_\nu(a,b) = Q_\Delta(a,b) + \mathrm{exp}(-(a^2 + b^2)/2) \times \biggl[\biggl(\frac{b}{a}\biggr)^{\Delta} \mathrm{I}_{\Delta}(ab) + \biggl(\frac{b}{a}\biggr)^{\Delta-1} \mathrm{I}_{\Delta-1}(ab)\biggr]
}
This result produces incompatible cumulative distribution functions of the distribution using \eqn{Q_\nu(a,b)} for \eqn{-1 < \nu < 1}. Therefore, the author of \pkg{lmomco} concludes that Shi (2012) is in error (or your author misinterprets the summation notation) and that the specific condition for \eqn{\nu = 0} shown above and lacking \eqn{\sum} is correct; there are three individual and separate conditions to support the Marcum Q function using the \dQuote{delta nu method}: \eqn{\nu \le -1}, \eqn{-1 < \nu < 1}, and \eqn{\nu \ge -1}.
}
\usage{
cdfkmu(x, para, paracheck=TRUE, getmed=TRUE, qualo=NA, quahi=NA,
                marcumQ=TRUE, marcumQmethod=c("chisq", "delta", "integral"))
}
\arguments{
  \item{x}{A real value vector.}
  \item{para}{The parameters from \code{\link{parkmu}} or  \code{\link{vec2par}}.}
  \item{paracheck}{A logical controlling whether the parameters and checked for validity.}
  \item{getmed}{Numerical problems rolling onto the distribution from the right can result in erroneous \eqn{F} being integrated of \code{\link{pdfkmu}}. This option is used to interrupt recurrsion, but if \code{TRUE}, then the median will be computed and for those \eqn{x} values less than the median and \eqn{F} initially computing as greater than 50 percent, are reset to 0.  Users are unlikely to need this option changed. But the hack can be turned off by setting \code{getmed=FALSE} as the user level.}
  \item{qualo}{A lower limit of the range of \eqn{x} to look for a \code{uniroot} of \eqn{F(x) = 0.5} to estimate the median quantile that is used to mitigate for erroneous numerical results. This argument is passed along to \code{\link{quakmu}} but also used as a truncation point for which \eqn{F=1} is returned if \eqn{x <} \code{qualo}. Lastly, see the last example below.}
    \item{quahi}{An upper limit of the range of \eqn{x} to look for a \code{uniroot} of \eqn{F(x) = 0.5} to estimate the median quantile that is used to mitigate for erroneous numerical results. This argument is passed along to \code{\link{quakmu}} but also used as a truncation point for which \eqn{F=1} is returned if \eqn{x >} \code{quahi}. Lastly, see the last example below.}
  \item{marcumQ}{A logical controlling whether the Marcum Q function is used instead of numerical integration of \code{\link{pdfkmu}}.}
  \item{marcumQmethod}{Which method for Marcum Q computation is to be used (see source code).}
}
\value{
  Nonexceedance probability (\eqn{F}) for \eqn{x}.
}
\note{
Code developed from Weinberg (2006). The biascor feature is of my own devise and this Poisson method does not seem to accommodate nu < 1 although Chornoboy claims valid for non-negative integer. The example implementation here will continue to use real values of nu.

\preformatted{
See NEWS file and entries for version 2.0.1 for this "R Marcum"
"marcumq" <- function(a, b, nu=1) {
	      pchisq(b^2, df=2*nu, ncp=a^2, lower.tail=FALSE) }

"marcumq.poissons" <-
   function(a,b, nu=NULL, nsim=10000, biascor=0.5) {
   asint <- as.logical(nu %% as.integer(nu))
   biascor <- ifelse(! asint, 0, biascor)
   marcumQint <- marcumq(a, b, nu=nu)
   B <- rpois(nsim, b^2/2)
   A <- nu - 1 + biascor + rpois(nsim, a^2/2)
   L <- B <= A
   marcumQppois <- length(L[L == TRUE])/nsim
   z <- list(MarcumQ.by.usingR = marcumQint,
             MarcumQ.by.poisson = marcumQppois)
   return(z)
}
x <- y <- vector()
for(i in 1:10000) {
   nu <- i/100
   z <- marcumq.poissons(12.4, 12.5, nu=nu)
   x[i] <- z$MarcumQ.by.usingR
   y[i] <- z$MarcumQ.by.poisson
}
plot(x,y, pch=16, col=rgb(x,0,0,.2),
     xlab="Marcum Q-function using R (ChiSq distribution)",
     ylab="Marcum Q-function by two Poisson random variables")
abline(0,1, lty=2)
}
}
\references{
Shi, Q., 2012, Semi-infinite Gauss-Hermite quadrature based approximations to the generalized Marcum and Nuttall Q-functions and further applications: First IEEE International Conference on Communications in China---Communications Theory and Security (CTS), pp. 268--273, ISBN 978--1--4673--2815--9,12.

Weinberg, G.V., 2006, Poisson representation and Monte Carlo estimation of generalized Marcum Q-function: IEEE Transactions on Aerospace and Electronic Systems, v. 42, no. 4, pp. 1520--1531.

Yacoub, M.D., 2007, The kappa-mu distribution and the eta-mu distribution: IEEE Antennas and Propagation Magazine, v. 49, no. 1, pp. 68--81.
}
\author{W.H. Asquith}
\seealso{\code{\link{pdfkmu}}, \code{\link{quakmu}}, \code{\link{lmomkmu}}, \code{\link{parkmu}}}
\examples{
\dontrun{
x <- seq(0,3, by=0.5)
para <- vec2par(c(0.69, 0.625), type="kmu")
cdfkmu(x, para, marcumQ=TRUE, marcumQmethod="chisq")
cdfkmu(x, para, marcumQ=TRUE, marcumQmethod="delta")
cdfkmu(x, para, marcumQ=FALSE) # about 3 times slower}
\dontrun{
para <- vec2par(c(0.69, 0.625), type="kmu")
quahi <- supdist(para, delexp=.1)$support[2]
cdfkmu(quahi, para, quahi=quahi)
}
\dontrun{
delx <- 0.01
x <- seq(0,3, by=delx)

plot(c(0,3), c(0,1), xlab="RHO", ylab="cdfkmu(RHO)", type="n")
para <- list(para=c(0, 0.75), type="kmu")
cdf <- cdfkmu(x, para)
lines(x, cdf, col=2, lwd=4)
para <- list(para=c(1, 0.5625), type="kmu")
cdf <- cdfkmu(x, para)
lines(x, cdf, col=3, lwd=4)

kappas <- c(0.00000001, 0.69, 1.37,  2.41, 4.45, 10.48, 28.49)
mus    <- c(0.75, 0.625,  0.5,  0.375, 0.25,  0.125, 0.05)
for(i in 1:length(kappas)) {
   kappa <- kappas[i]
   mu    <- mus[i]
   para <- list(para=c(kappa, mu), type="kmu")
   cdf <- cdfkmu(x, para)
   lines(x, cdf, col=i)
}
}
\dontrun{
delx <- 0.005
x <- seq(0,3, by=delx)
nx <- 20*log10(x)
plot(c(-30,10), 10^c(-4,0), log="y", xaxs="i", yaxs="i",
     xlab="RHO", ylab="cdfkmu(RHO)", type="n")
m <- 1.25
mus <- c(0.25, 0.50, 0.75, 1, 1.25, 0)
for(mu in mus) {
   col <- 1
   kappa <- m/mu - 1 + sqrt((m/mu)*((m/mu)-1))
   para <- vec2par(c(kappa, mu), type="kmu")
   if(! is.finite(kappa)) {
      para <- vec2par(c(Inf,m), type="kmu")
      col <- 2
   }
   lines(nx, cdfkmu(x, para), col=col)
}
mtext("Yacoub (2007, figure 4)")
}
\dontrun{
# The Marcum Q use for the CDF avoid numerical integration of pdfkmu(), but
# below is an example for which there is some failure that remains to be found.
para <- vec2par(c(10, 23), type="kmu")
# The following are reliable but slower as they avoid the Marcum Q function
# and use traditional numerical integration of the PDF function.
A <- cdfkmu(c(0.10, 0.35, 0.9, 1, 1.16), para, marcumQ=FALSE)
# Continuing, the first value in c() has an erroneous value for the next call.
B <- cdfkmu(c(0.10, 0.35, 0.9, 1, 1.16), para, marcumQ=TRUE)
# But this distribution is tightly peaks and well away from the origin, so in
# order to snap the erroneous value to zero, we need a successful median
# computation.  We can try again using the qualo argument to pass through to
# quakmu() like the following:
C <- cdfkmu(c(0.10, 0.35, 0.9, 1, 1.16), para, marcumQ=TRUE, qualo=0.4)
# The existance of the median for the last one also triggers a truncation of
# the CDF to 0 when negative solution results for the 0.35, although the
# negative is about -1E-14.
}
\dontrun{
# Does the discipline of the signal litature just "know" about the apparent
# upper support of the Kappa-Mu being quite near or even at pi?
"simKMU" <- function() {
   km <- 10^runif(2, min=-3, max=3)
   f <- cdfkmu(pi, vec2par(km, type="kmu"))
   return(c(km, f))
}
EndStudy <- sapply(1:1000, function(i) { simKMU() } )
boxplot(EndStudy[3,])
}
}
\keyword{distribution}
\keyword{cumulative distribution function}
\keyword{Distribution: Kappa-Mu}
