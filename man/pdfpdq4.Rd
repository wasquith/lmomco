\encoding{utf8}
\name{pdfpdq4}
\alias{pdfpdq4}
\title{Probability Density Function of the Polynomial Density-Quantile4 Distribution}
\description{
This function computes the probability density of the Polynomial Density-Quantile4 distribution given parameters (\eqn{\alpha} and \eqn{\beta}) computed by \code{\link{parpdq4}}. The probability density function has not explicit form. The implementation here simply uses a \emph{five-point stencil} to approciate the derivative of the cumulative distribution function \code{\link{cdfpdq4}} and hence an \code{eps} term is used and multipled to the scale parameter (\eqn{\alpha}) of the distribution. The distribution's canonical definition is in terms of the quantile function (\code{\link{quapdq4}}).
}
\usage{
pdfpdq4(x, para, paracheck=TRUE, h=NA, hfactor=0.2)
}
\arguments{
  \item{x}{A real value vector.}
  \item{para}{The parameters from \code{\link{parpdq4}} or \code{\link{vec2par}}.}
  \item{paracheck}{A logical switch as to whether the validity of the parameters should be checked. Default is \code{paracheck=TRUE}. This switch is made so that the root solution needed for \code{\link{cdfpdq4}} shows an extreme speed increase because of the repeated calls to \code{quapdq4}.}
  \item{h}{The differential element of the stencil, if provided, otherwise \code{hfactor} used.}
  \item{hfactor}{A term multiplied to the \eqn{\alpha} parameter to set the \eqn{h} in the numerical derivative. Not optimal, but seems to work for a variety of chosen parameters for plotting the density function.}
}
\value{
  Probability density (\eqn{f}) for \eqn{x}.
}
\references{
Hosking, J.R.M., 2007, Distributions with maximum entropy subject to constraints on their L-moments or expected order statistics: Journal of Statistical Planning and Inference, v. 137, no. 9, pp. 2870--2891, \doi{10.1016/j.jspi.2006.10.010}.
}
\author{W.H. Asquith}
\seealso{\code{\link{cdfpdq4}}, \code{\link{quapdq4}}, \code{\link{lmompdq4}}, \code{\link{parpdq4}}}
\examples{
\dontrun{
  para <- list(para=c(0, 0.4332, -0.7029), type="pdq4")
  X <- seq(-4, +4, by=(4 - -4) / 1000)
  plot( X, pdfpdq4(X, para), type="l", col=grey(0.8), lwd=4, ylim=c(0, 0.5))
  lines(X, dnorm(  X, sd=1), lty=2)
  legend("topleft", c("Standard normal distribution",
                      "PDQ4 distribution with same L-moments as the standard normal"),
                    lwd=c(1, 4), lty=c(2, 1), col=c(1, grey(0.8)), cex=0.8)
  mtext("Mimic Hosking (2007, fig. 3 [left])")
  check.pdf(pdfpdq4, para, hfactor=0.3) }

\dontrun{
  para <- list(para=c(100, 43.32, -0.7029), type="pdq4")
  minX <- quapdq4(0.0001, para)
  maxX <- quapdq4(0.9999, para)
  X <- seq(minX, maxX, by=(maxX - minX) / 1000)
  plot( X, pdfpdq4(X, para), type="l", col=grey(0.8), lwd=4)

  check.pdf(pdfpdq4, para, hfactor=0.3) }
}
\keyword{distribution}
\keyword{probability density function}
\keyword{Distribution: Polynomial Density-Quantile4}
