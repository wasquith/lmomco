\encoding{utf8}
\name{cdfpdq4}
\alias{cdfpdq4}
\title{Cumulative Distribution Function of the Polynomial Density-Quantile4 Distribution}
\description{
This function computes the cumulative probability or nonexceedance probability of the Polynomial Density-Quantile4 (PDQ4) distribution given parameters (\eqn{\xi}, \eqn{\alpha}, \eqn{\kappa}) computed by \code{\link{parpdq4}}. The cumulative distribution function has no explicit form and requires numerical methods. The \R function \code{uniroot()} is used to root the quantile function \code{\link{quapdq4}} to compute the nonexceedance probability. The distribution's canonical definition is in terms of the quantile function (\code{\link{quapdq4}}).
}
\usage{
cdfpdq4(x, para, paracheck=TRUE)
}
\arguments{
  \item{x}{A real value vector.}
  \item{para}{The parameters from \code{\link{parpdq4}} or  \code{\link{vec2par}}.}
  \item{paracheck}{A logical switch as to whether the validity of the parameters should be checked. Default is \code{paracheck=TRUE}. This switch is made so that the root solution needed for \code{\link{cdfpdq4}} shows an extreme speed increase because of the repeated calls to \code{quapdq4}.}
}
\value{
  Nonexceedance probability (\eqn{F}) for \eqn{x}.
}
\references{
Hosking, J.R.M., 2007, Distributions with maximum entropy subject to constraints on their L-moments or expected order statistics: Journal of Statistical Planning and Inference, v. 137, no. 9, pp. 2870--2891, \doi{10.1016/j.jspi.2006.10.010}.
}
\author{W.H. Asquith}
\seealso{\code{\link{pdfpdq4}}, \code{\link{quapdq4}}, \code{\link{lmompdq4}}, \code{\link{parpdq4}}}
\examples{
\dontrun{
  FF <- seq(0.001, 0.999, by=0.001)
  para  <- list(para=c(0, 0.4332, -0.7029), type="pdq4")
  Fpdq4 <- cdfpdq4(quapdq4(FF, para), para)
  plot(FF, Fpdq4, type="l", col=grey(0.8), lwd=4)
  # should be a 1:1 line, it is }

\dontrun{
  para <- list(para=c(0, 0.4332, -0.7029), type="pdq4")
  X <- seq(-5, +12, by=(12 - -5) / 500)
  plot( X, cdfpdq4(X, para), type="l", col=grey(0.8), lwd=4, ylim=c(0, 1))
  lines(X, pf( exp(X), df1=5, df2=4), lty=2)
  lines(X, c(NA, diff( cdfpdq4(X, para))          / ((12 - -5) / 500)))
  lines(X, c(NA, diff( pf( exp(X), df1=5, df2=4)) / ((12 - -5) / 500)), lty=2) # }
}
\keyword{distribution}
\keyword{cumulative distribution function}
\keyword{Distribution: Polynomial Density-Quantile4}
