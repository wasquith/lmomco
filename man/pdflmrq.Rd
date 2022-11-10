\name{pdflmrq}
\alias{pdflmrq}
\title{Probability Density Function of the Linear Mean Residual Quantile Function Distribution
}
\description{
This function computes the probability density function of the Linear Mean Residual Quantile Function distribution given parameters  computed by \code{\link{parlmrq}}.  The probability density function is
\deqn{f(x) = \frac{1 - F(x)}{2\alpha\,F(x) + (\mu - \alpha)}\mbox{,}}
where \eqn{f(x)} is the nonexceedance probability for quantile \eqn{x}, \eqn{F(x)} is the cumulative distribution function or nonexceedance probability at \eqn{x}, \eqn{\mu} is a location parameter, and \eqn{\alpha} is a scale parameter.
}
\usage{
pdflmrq(x, para)
}
\arguments{
  \item{x}{A real value vector.}
  \item{para}{The parameters from \code{\link{parlmrq}} or \code{\link{vec2par}}.}
}

\value{
  Probability density (\eqn{f}) for \eqn{x}.
}
\references{
Midhu, N.N., Sankaran, P.G., and Nair, N.U., 2013, A class of distributions with linear mean residual quantile function and it's generalizations: Statistical Methodology, v. 15, pp. 1--24.
}
\author{W.H. Asquith}
\seealso{\code{\link{cdflmrq}}, \code{\link{qualmrq}}, \code{\link{lmomlmrq}}, \code{\link{parlmrq}}
}
\examples{
lmr <- lmoms(c(3, 0.05, 1.6, 1.37, 0.57, 0.36, 2.2))
pdflmrq(3,parlmrq(lmr))
\dontrun{
para.lmrq <- list(para=c(2.1043, 0.4679), type="lmrq")
para.wei  <- vec2par(c(0,2,0.9), type="wei") # note switch from Midhu et al. ordering.
F <- seq(0.01,0.99,by=.01); x <- qualmrq(F, para.lmrq)
plot(x, pdflmrq(x, para.lmrq), type="l", ylab="", lwd=2, lty=2, col=2,
     xlab="The p.d.f. of Weibull and p.d.f. of LMRQD", xaxs="i", yaxs="i",
     xlim=c(0,9), ylim=c(0,0.8))
lines(x, pdfwei(x, para.wei))
mtext("Midhu et al. (2013, Statis. Meth.)")
}
}
\keyword{distribution}
\keyword{probability density function}
\keyword{Distribution: Linear Mean Residual Quantile Function}
