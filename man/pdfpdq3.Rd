\encoding{utf8}
\name{pdfpdq3}
\alias{pdfpdq3}
\title{Probability Density Function of the Polynomial Density-Quantile3 Distribution}
\description{
This function computes the probability density of the Polynomial Density-Quantile3 distribution given parameters (\eqn{\alpha} and \eqn{\beta})  computed by \code{\link{parpdq3}}. The probability density function has not explicit form. The implementation here simply uses a \emph{five-point stencil} to approciate the derivative of the cumulative distribution function \code{\link{cdfpdq3}} and hence an \code{eps} term is used and multipled to the scale parameter (\eqn{\alpha}) of the distribution. The distribution's canonical definition is in terms of the quantile function (\code{\link{quapdq3}}).
}
\usage{
pdfpdq3(x, para, paracheck=TRUE, h=NA, hfactor=0.2)
}
\arguments{
  \item{x}{A real value vector.}
  \item{para}{The parameters from \code{\link{parpdq4}} or \code{\link{vec2par}}.}
  \item{paracheck}{A logical switch as to whether the validity of the parameters should be checked. Default is \code{paracheck=TRUE}. This switch is made so that the root solution needed for \code{\link{cdfpdq3}} shows an extreme speed increase because of the repeated calls to \code{quapdq3}.}
    \item{h}{The differential element of the stencil, if provided, otherwise \code{hfactor} used.}
  \item{hfactor}{A term multiplied to the \eqn{\alpha} parameter to set the \eqn{h} in the numerical derivative. Not optimal, but seems to work for a variety of chosen parameters for plotting the density function.}
}
\value{
  Probability density (\eqn{f}) for \eqn{x}.
}
\references{
Hosking, J.R.M., 2007, Distributions with maximum entropy subject to
constraints on their L-moments or expected order statistics: Journal of
Statistical Planning and Inference, v. 137, no. 9, pp. 2870--2891, \doi{10.1016/j.jspi.2006.10.010}.
}
\author{W.H. Asquith}
\seealso{\code{\link{cdfpdq3}}, \code{\link{quapdq3}}, \code{\link{lmompdq3}}, \code{\link{parpdq3}}}
\examples{
\dontrun{
  para <- list(para=c(0.6933, 1.5495, 0.5488), type="pdq3")
  X <- seq(-5, +12, by=(12 - -5) / 1000)
  plot( X, pdfpdq3(X, para), type="l", col=grey(0.8), lwd=4, ylim=c(0, 0.3))
  lines(X, c(NA, diff(pf(exp(X), df1=7, df2=1))/((12 - -5) / 1000)), lty=2)
  legend("topleft", c("log F(7,1) distribution with same L-moments",
                      "PDQ3 distribution with same L-moments as the log F(7,1)"),
                    lwd=c(1, 4), lty=c(2, 1), col=c(1, grey(0.8)), cex=0.8)
  mtext("Mimic Hosking (2007, fig. 2 [left])")
  check.pdf(pdfpdq3, para) # }

\dontrun{
  para <- list(para=c(100, 43.32, -0.7029), type="pdq3")
  minX <- quapdq3(0.0001, para)
  maxX <- quapdq3(0.9999, para)
  X <- seq(minX, maxX, by=(maxX - minX) / 1000)
  plot( X, pdfpdq3(X, para), type="l", col=grey(0.8), lwd=4)
  check.pdf(pdfpdq3, para) # }

\dontrun{
  para <- vec2par(c(0.4729820, 3.0242067, 0.9880701), type="pdq3")
  print(lmom2par(par2lmom(para), type="pdq3"))
  # "|kappa| > 0.98, alpha (yes alpha) results could be unreliable"
  # So, we are entering into a problem for which the kappa parameter is
  # very large and instabilities in the algorithm will result, but
  # vec2par() has not mechanism for determining this type of situation.
  # Ultimately, things will manifest with a check.pdf() that fails.
  sup <- lmomco::supdist(para)$support
  xx <- seq(sup[1], sup[2], by=diff(range(sup)) / 2000)
  plot(xx, pdfpdq3(xx, para), type="l", col=grey(0.8))
  plot(xx, pdfpdq3(xx, para), type="l", col=grey(0.8), xlim=c(-1,10))
  # See hints of instability in the density shape in the second plot
  check.pdf(pdfpdq3, para) # non-finite function value }
}
\keyword{distribution}
\keyword{probability density function}
\keyword{Distribution: Polynomial Density-Quantile3}
