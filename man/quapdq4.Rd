\encoding{utf8}
\name{quapdq4}
\alias{quapdq4}
\title{Quantile Function of the Polynomial Density-Quantile4 Distribution}
\description{
This function computes the quantiles of the Polynomial Density-Quantile4
distribution (PDQ4) given parameters (\eqn{\xi}, \eqn{\alpha}, and
\eqn{\kappa}) computed by \code{\link{parpdq4}}. The quantile function
for \eqn{0 < \kappa < 1} is
\deqn{x(F) = \xi + \alpha \biggl[\log\biggl(\frac{F}{1-F}\biggr) -
  2\kappa\;\mathrm{atanh}(\kappa[2F-1])\biggr] \mbox{\ and}}
for \eqn{-\infty < \kappa < 0} is
\deqn{x(F) = \xi + \alpha \biggl[\log\biggl(\frac{F}{1-F}\biggr) + 2\kappa\;\mathrm{atan}(\kappa[2F-1])\biggr] \mbox{,}}
where \eqn{x(F)} is the quantile for nonexceedance probability \eqn{F},
\eqn{\xi} is a location parameter, \eqn{\alpha} is a scale parameter,
and \eqn{\kappa} is a shape parameter. The range of the distribution is \eqn{-\infty < x < \infty}.
}
\usage{
quapdq4(f, para, paracheck=TRUE)
}
\arguments{
  \item{f}{Nonexceedance probability (\eqn{0 \le F \le 1}).}
  \item{para}{The parameters from \code{\link{parpdq4}} or \code{\link{vec2par}}.}
  \item{paracheck}{A logical controlling whether the parameters are checked for validity. Overriding of this check might be extremely important and needed for use of the quantile function in the context of TL-moments with nonzero trimming.}
}
\value{
  Quantile value for nonexceedance probability \eqn{F}.
}
\details{
The PDQ4 was proposed by Hosking (2007) with the core justification of maximizing entropy and that \dQuote{maximizing entropy subject to a set of constraints can be regarded as deriving a distribution that is consistent with the information specified in the constraints while making minimal assumptions about the form of the distribution other than those embodied in the constraints.} The PDQ4 is that family constrained to the \eqn{\lambda_1}, \eqn{\lambda_2}, and \eqn{\tau_4} values of the L-moments. (See also the Polynomial Density-Quantile3 function for constraint on \eqn{\lambda_1}, \eqn{\lambda_2}, and \eqn{\tau_3} values of the L-moments, \code{\link{quapdq3}}.)

The PDQ4 is a symmetrical distribution (\eqn{\tau_3 = 0} everywhere) that has maximum entropy conditional on having specified values for the L-moments of \eqn{\lambda_1}, \eqn{\lambda_2}, and \eqn{\lambda_4 = \tau_4\lambda_2} with \eqn{\lambda_3 = \tau_3 = 0}.  The tails of the PDQ4 are exponentially decreasing and the distribution could be useful in distributional analysis with data showing similar tail characteristics. The attainable L-kurtosis range is \eqn{-1/4 < \tau_4 < 1} with the sign change from negative to positive of \eqn{\kappa} occurring at \eqn{\tau_4 = 1/6}. Finally, PDQ4 generalizes the logistic distribution, which is the special case \eqn{\kappa \rightarrow 0}, and contains distributions both lighter-tailed (\eqn{\kappa < 0}) and heavier-tailed (\eqn{\kappa > 0}) than the logistic.
}
\references{
Hosking, J.R.M., 2007, Distributions with maximum entropy subject to constraints on their L-moments or expected order statistics: Journal of Statistical Planning and Inference, v. 137, no. 9, pp. 2,870--2891, \doi{10.1016/j.jspi.2006.10.010}.
}
\author{W.H. Asquith}
\seealso{\code{\link{cdfpdq4}}, \code{\link{pdfpdq4}}, \code{\link{lmompdq4}}, \code{\link{parpdq4}}}
\examples{
  lmr <- lmoms(c(123, 34, 4, 654, 37, 78))
  quapdq4(0.5, parpdq4(lmr)) # [1] 155

\dontrun{
  FF <- seq(0.0001, 0.9999, by=0.001)
  para <- list(para=c(0, 0.4332, -0.7029), type="pdq4")
  plot( qnorm(FF, sd=1), quapdq4(FF, para), type="l", col=grey(0.8), lwd=4,
       xlab="Standard normal variate", ylab="Quantiles, Q(f)")
  lines(qnorm(FF, sd=1),   qnorm(FF, sd=1), lty=2)
  legend("topleft", c("Standard normal distribution",
                      "PDQ4 distribution with same L-moments as the standard normal"),
        lwd=c(1, 4), lty=c(2, 1), col=c(1, grey(0.8)), cex=0.8)
  mtext("Mimic Hosking (2007, fig. 3 [right])") # }

\dontrun{
  # A quick recipe to look at the shapes of quantile functions.
  FF <- seq(0.001, 0.999, by=0.001)
  plot(qnorm(FF), qnorm(FF), type="n", ylim=c(-7, 7),
       xlab="Standard normal variate", ylab="PDQ4 variate")
  abline(h=0, lty=2, lwd=0.9); abline(v=0, lty=2, lwd=0.9)

  lscale   <- 1 / sqrt(pi)
  tau4s    <- seq(-1/4, 0.7, by=.05)
  tau4s[1] <- tau4s[1] + 0.001
  for(i in 1:length(tau4s)) {
    lmr <- vec2lmom(c(0, lscale, 0, tau4s[i]))
    if(! are.lmom.valid(lmr)) next
    pdq4 <- parpdq4(lmr, snapt4uplimit=FALSE)
    lines(qnorm(FF), qlmomco(FF, pdq4), col=rgb(abs(tau4s[i]), 0, 1))
  }
  abline(0,1, col="darkgreen", lwd=3)
  txt <- "Standard normal distribution (Tau4=0.122602)"
  txt <- c(txt, paste0("PDQ4 distribution for varying Tau4 values",
                       " (color varies for accenting)"))
  legend("topleft", txt, col=c("darkgreen", rgb(0.2, 0, 1)),
                         cex=0.9, bty="n", lwd=c(3,1)) # }
}
\keyword{distribution}
\keyword{quantile function}
\keyword{Distribution: Polynomial Density-Quantile4}
