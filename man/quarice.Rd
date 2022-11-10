\name{quarice}
\alias{quarice}
\title{Quantile Function of the Rice Distribution}
\description{
This function computes the quantiles of the Rice distribution given parameters (\eqn{\nu} and \eqn{\alpha})  computed by \code{\link{parrice}}. The quantile function  is complex and numerical rooting of the cumulative distribution function \code{\link{cdfrice}} is used.
}
\usage{
quarice(f, para, xmax=NULL, paracheck=TRUE)
}
\arguments{
  \item{f}{Nonexceedance probability (\eqn{0 \le F \le 1}).}
  \item{para}{The parameters from \code{\link{parrice}} or \code{\link{vec2par}}.}
  \item{xmax}{The maximum x value used for integeration.}
  \item{paracheck}{A logical controlling whether the parameters are checked for validity. Overriding of this check might be extremely important and needed for use of the quantile function in the context of TL-moments with nonzero trimming.}
}
\value{
  Quantile value for nonexceedance probability \eqn{F}.
}
\author{W.H. Asquith}
\references{
Asquith, W.H., 2011, Distributional analysis with L-moment statistics using the R environment for statistical computing: Createspace Independent Publishing Platform, ISBN 978--146350841--8.
}
\seealso{\code{\link{cdfrice}}, \code{\link{pdfrice}}, \code{\link{lmomrice}}, \code{\link{parrice}}}
\examples{
lmr <- vec2lmom(c(125,0.20), lscale=FALSE)
quarice(0.75,parrice(lmr))
# The quantile function of the Rice as implemented in lmomco
# is slow because of rooting the CDF, which is created by
# integration of the PDF. Rician random variates are easily created.
# Thus, in speed applications the rlmomco() with a Rice parameter
# object could be bypassed by the following function, rrice().
\dontrun{
"rrice" = function(n, nu, alpha) { # from the VGAM package
    theta = 1 # any number
    X = rnorm(n, mean=nu * cos(theta), sd=alpha)
    Y = rnorm(n, mean=nu * sin(theta), sd=alpha)
    return(sqrt(X^2 + Y^2))
}
n <- 5000; # suggest making it about 10,000
nu <- 100; alpha <- 10
set.seed(501); lmoms(rrice(n, nu, alpha))
set.seed(501); lmoms(rlmomco(n, vec2par(c(nu,alpha), type='rice')))
# There are slight numerical differences between the two?
}
}
\keyword{distribution}
\keyword{quantile function}
\keyword{Distribution: Rice}
