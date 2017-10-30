\name{fpds2f}
\alias{fpds2f}
\title{Conversion of Partial-Duration Nonexceedance Probability to Annual Nonexceedance Probability}
\description{
This function takes partial duration series nonexceedance probability and converts it to a  an annual exceedance probability through a simple assumption that the Poisson distribution is appropriate. The relation between the cumulative distribution function \eqn{F(x)} for the annual series is related to the cumulative distribution function \eqn{G(x)} of the partial-duration series by
\deqn{F(x) = \mathrm{exp}(-\eta[1 - G(x)])\mathrm{.}}
The core assumption is that successive events in the partial-duration series can be considered as \emph{independent}. The \eqn{\eta} term is the arrival rate of the events. For example, suppose that 21 events have occurred in 15 years, then \eqn{\eta = 21/15 = 1.4} events per year.

The example documented here provides a comprehensive demonstration of the function along with a partnering function \code{\link{f2fpds}}. That function performs the opposite conversion. Lastly, the cross reference to \code{\link{x2xlo}} is made because the example contained therein provides another demonstration of partial-duration and annual series frequency analysis.
}
\usage{
fpds2f(fpds, rate=NA)
}
\arguments{
  \item{fpds}{A vector of partial-duration nonexceedance probabilities.}
  \item{rate}{The number of events per year.}
}
\value{
  A vector of converted nonexceedance probabilities.
}
\references{
Stedinger, J.R., Vogel, R.M., Foufoula-Georgiou, E., 1993, Frequency analysis of extreme events: \emph{in} Handbook of Hydrology, ed. Maidment, D.R., McGraw-Hill, Section 18.6 Partial duration series, mixtures, and censored data, pp. 18.37--18.39.
}
\author{ W.H. Asquith}
\seealso{ \code{\link{f2fpds}}, \code{\link{x2xlo}}, \code{\link{f2flo}}, \code{\link{flo2f}} }
\examples{
\dontrun{
stream <- "A Stream in West Texas"
Qpds    <- c(61.8, 122, 47.3, 71.1, 211, 139, 244, 111, 233, 102)
Qann <- c(61.8, 122, 71.1, 211, 244, 0, 233)
years  <- length(Qann)  # gage has operated for about 7 years
visits <- 27  # number of visits or "events"
rate   <- visits/years
Z <- rep(0, visits-length(Qpds))
Qpds <- c(Qpds,Z) # The creation of a partial duration series
# that will contain numerous zero values.

Fs <- seq(0.001,.999, by=.005) # used to generate curves

type <- "pe3" # The Pearson type III distribution
PPpds <- pp(Qpds); Qpds <- sort(Qpds) # plotting positions (partials)
PPann <- pp(Qann); Qann <- sort(Qann) # plotting positions (annuals)
parann <- lmom2par(lmoms(Qann), type=type) # parameter estimation (annuals)
parpsd <- lmom2par(lmoms(Qpds), type=type) # parameter estimation (partials)

Fsplot    <- qnorm(Fs) # in order to produce normal probability paper
PPpdsplot <- qnorm(fpds2f(PPpds, rate=rate)) # ditto
PPannplot <- qnorm(PPann) # ditto

# There are many zero values in this particular data set that require leaving
# them out in order to achieve appropriate curvature of the Pearson type III
# distribution. Conditional probability adjustments will be used.
Qlo <- x2xlo(Qpds) # Create a left out object with an implied threshold of zero
parlo <- lmom2par(lmoms(Qlo$xin), type=type) # parameter estimation for the
# partial duration series values that are greater than the threshold, which
# defaults to zero.

plot(PPpdsplot, Qpds, type="n", ylim=c(0,400), xlim=qnorm(c(.01,.99)),
     xlab="STANDARD NORMAL VARIATE", ylab="DISCHARGE, IN CUBIC FEET PER SECOND")
mtext(stream)
points(PPannplot, Qann, col=3, cex=2, lwd=2, pch=0)
points(qnorm(fpds2f(PPpds, rate=rate)), Qpds, pch=16, cex=0.5 )
points(qnorm(fpds2f(flo2f(pp(Qlo$xin), pp=Qlo$pp), rate=rate)),
       sort(Qlo$xin), col=2, lwd=2, cex=1.5, pch=1)
points(qnorm(fpds2f(Qlo$ppout, rate=rate)),
       Qlo$xout, pch=4, col=4)

lines(qnorm(fpds2f(Fs, rate=rate)),
      qlmomco(Fs, parpsd), lwd=1, lty=2)
lines(Fsplot, qlmomco(Fs, parann), col=3, lwd=2)
lines(qnorm(fpds2f(flo2f(Fs, pp=Qlo$pp), rate=rate)),
      qlmomco(Fs, parlo), col=2, lwd=3)

# The following represents a subtle application of the probability transform
# functions. The show how one starts with annual recurrence intervals
# converts into conventional annual nonexceedance probabilities as well as
# converting these to the partial duration nonexceedance probabilities.
Tann <- c(2, 5, 10, 25, 50, 100)
Fann <- T2prob(Tann); Gpds <- f2fpds(Fann, rate=rate)
FFpds <- qlmomco(f2flo(Gpds, pp=Qlo$pp), parlo)
FFann <- qlmomco(Fann, parann)
points(qnorm(Fann), FFpds, col=2, pch=16)
points(qnorm(Fann), FFann, col=3, pch=16)

legend(-2.4,400, c("True annual series (with one zero year)",
                "Partial duration series (including 'visits' as 'events')",
                "Partial duration series (after conditional adjustment)",
                "Left-out values (<= zero) (trigger of conditional probability)",
                "PE3 partial-duration frequency curve (PE3-PDS)",
                "PE3 annual-series frequency curve (PE3-ANN)",
                "PE3 partial-duration frequency curve (zeros removed) (PE3-PDSz)",
                "PE3-ANN  T-year event: 2, 5, 10, 25, 50, 100 years",
                "PE3-PDSz T-year event: 2, 5, 10, 25, 50, 100 years"),
       bty="n", cex=.75,
       pch=c(0,  16, 1, 4, NA, NA, NA, 16, 16),
       col=c(3,  1, 2,  4,  1,  3,  2,  3, 2),
       pt.lwd=c(2,1,2,1), pt.cex=c(2, 0.5, 1.5, 1, NA, NA, NA, 1, 1),
       lwd=c(0,0,0,0,1,2,3), lty=c(0,0,0,0,2,1,1))
}
}
\keyword{utility (probability)}
