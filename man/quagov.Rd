\name{quagov}
\alias{quagov}
\title{Quantile Function of the Govindarajulu Distribution}
\description{
This function computes the quantiles of the Govindarajulu distribution given parameters (\eqn{\xi}, \eqn{\alpha}, and \eqn{\beta})  computed by
\code{\link{pargov}}. The quantile function is
\deqn{x(F) = \xi + \alpha[(\beta+1)F^\beta - \beta F^{\beta+1}] \mbox{,}}
where \eqn{x(F)} is the quantile for nonexceedance probability \eqn{F}, \eqn{\xi} is location parameter, \eqn{\alpha} is a scale parameter, and \eqn{\beta} is a shape parameter.
}
\usage{
quagov(f, para, paracheck=TRUE)
}
\arguments{
  \item{f}{Nonexceedance probability (\eqn{0 \le F \le 1}).}
  \item{para}{The parameters from \code{\link{pargov}} or similar.}
  \item{paracheck}{A logical controlling whether the parameters are checked for validity. Overriding of this check might be extremely important and needed for use of the quantile function in the context of TL-moments with nonzero trimming.}
}
\value{
  Quantile value for for nonexceedance probability \eqn{F}.
}
\references{
Gilchrist, W.G., 2000, Statistical modelling with quantile functions: Chapman and Hall/CRC, Boca Raton.

Nair, N.U., Sankaran, P.G., Balakrishnan, N., 2013, Quantile-based reliability analysis: Springer, New York.

Nair, N.U., Sankaran, P.G., and Vineshkumar, B., 2012, The Govindarajulu distribution---Some Properties and applications: Communications in Statistics, Theory and Methods, 41(24), 4391--4406.
}
\author{ W.H. Asquith}
\seealso{\code{\link{cdfgov}}, \code{\link{pdfgov}}, \code{\link{lmomgov}}, \code{\link{pargov}}}
\examples{
lmr <- lmoms(c(123,34,4,654,37,78))
quagov(0.5,pargov(lmr))
\dontrun{
lmr <- lmoms(c(3, 0.05, 1.6, 1.37, 0.57, 0.36, 2.2));
par <- pargov(lmr)# LMRQ said to have a linear mean residual quantile function.
# Let us have a look.
F <- c(0,nonexceeds(),1)
plot(F, qlmomco(F,par), type="l", lwd=3, xlab="NONEXCEEDANCE PROBABILITY",
     ylab="LIFE TIME, RESIDUAL LIFE, OR REVERSED RESIDUAL LIFE")
lines(F, rmlmomco(F,par),  col=2, lwd=4)  # heavy red line (residual life)
lines(F, rrmlmomco(F,par), col=2, lty=2)  # dashed red (reversed res. life)
lines(F, cmlmomco(F,par),  col=4)         # conditional mean (blue)
# Notice how the conditional mean attaches to the parent at F=1, but it does not
# attached at F=0 because of the none zero origin.
cmlmomco(0,par)           # 1.307143 # expected life given birth only
lmomgov(par)$lambdas[1]   # 1.307143 # expected life of the parent distribution
rmlmomco(0, par)          # 1.288989 # residual life given birth only
qlmomco(0, par)           # 0.018153 # instantaneous life given birth
# Note: qlmomco(0,par) + rmlmomco(0,par) is the E[lifetime], but rmlmomco()
# is the RESIDUAL MEAN LIFE.
}
}
\keyword{distribution}
\keyword{quantile function}
\keyword{Distribution: Govindarajulu}
