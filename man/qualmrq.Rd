\name{qualmrq}
\alias{qualmrq}
\title{Quantile Function of the Linear Mean Residual Quantile Function Distribution
}
\description{
This function computes the quantiles of the Linear Mean Residual Quantile Function distribution given parameters (\eqn{\mu} and \eqn{\alpha})  computed by \code{\link{parlmrq}}. The quantile function  is
\deqn{x(F) = -(\alpha + \mu)\times\log(1-F) - 2\alpha\times F\mbox{,}}
where \eqn{x(F)} is the quantile for nonexceedance probability \eqn{F}, \eqn{\mu} is a location parameter, and \eqn{\alpha} is a scale parameter.  The parameters must satisfy \eqn{\mu > 0} and \eqn{-\mu \le \alpha < \mu}.
}
\usage{
qualmrq(f, para, paracheck=TRUE)
}
\arguments{
  \item{f}{Nonexceedance probability (\eqn{0 \le F \le 1}).}
  \item{para}{The parameters from \code{\link{parlmrq}} or \code{\link{vec2par}}.}
  \item{paracheck}{A logical controlling whether the parameters are checked for validity. Overriding of this check might be extremely important and needed for use of the quantile function in the context of TL-moments with nonzero trimming.}
}
\value{
  Quantile value for nonexceedance probability \eqn{F}.
}
\references{
Midhu, N.N., Sankaran, P.G., and Nair, N.U., 2013, A class of distributions with linear mean residual quantile function and it's generalizations: Statistical Methodology, v. 15, pp. 1--24.
}
\author{W.H. Asquith}
\seealso{\code{\link{cdflmrq}}, \code{\link{pdflmrq}}, \code{\link{lmomlmrq}}, \code{\link{parlmrq}}
}
\examples{
lmr <- lmoms(c(3, 0.05, 1.6, 1.37, 0.57, 0.36, 2.2));
par <- parlmrq(lmr)
qualmrq(0.75,par)  
\dontrun{
# The distribution is said to have a linear mean residual quantile function.
# Let us have a look.
F <- nonexceeds(); par <- vec2par(c(101,21), type="lmrq")
plot(F, qlmomco(F,par), type="l", lwd=3, xlab="NONEXCEEDANCE PROBABILITY",
     ylab="LIFE TIME, RESIDUAL LIFE, OR REVERSED RESIDUAL LIFE")
lines(F, rmlmomco(F,par),  col=2, lwd=4) # heavy red line (residual life)
lines(F, rrmlmomco(F,par), col=2, lty=2) # dashed red (reversed res. life)
lines(F, cmlmomco(F,par),  col=4)        # conditional mean (blue)
# Notice that the rmlmomco() is a straight line as the name of the parent
# distribution: Linear Mean Residual Quantile Distribution suggests.
# Curiously, the reversed mean residual is not linear.
}
}
\keyword{distribution}
\keyword{quantile function}
\keyword{Distribution: Linear Mean Residual Quantile Function}
