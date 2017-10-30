\name{cdfrice}
\alias{cdfrice}
\title{Cumulative Distribution Function of the Rice Distribution}
\description{
This function computes the cumulative probability or nonexceedance probability of the Rice distribution given parameters (\eqn{\nu} and \eqn{\mathrm{SNR}}) computed by \code{\link{parrice}}. The cumulative distribution function is complex and numerical integration of the probability density function \code{\link{pdfrice}} is used.
\deqn{
F(x) = 1 - Q\biggl(\frac{\nu}{\alpha}, \frac{x}{\alpha}\biggr)\mbox{,}
}
where \eqn{F(x)} is the nonexceedance probability for quantile \eqn{x}, \eqn{Q(a,b)} is the Marcum Q-function, and \eqn{\nu/\alpha} is a form of signal-to-noise ratio \eqn{\mathrm{SNR}}. If \eqn{\nu=0}, then the Rayleigh distribution results and \code{\link{pdfray}} is used. The Marcum Q-function is difficult to work with and the \pkg{lmomco} uses the \code{integrate} function on \code{\link{pdfrice}} (however, see the Note).
}
\usage{
cdfrice(x, para)
}
\arguments{
  \item{x}{A real value vector.}
  \item{para}{The parameters from \code{\link{parrice}} or  \code{\link{vec2par}}.}
}
\value{
  Nonexceedance probability (\eqn{F}) for \eqn{x}.
}
\note{
A user of \pkg{lmomco} reported that the Marcum Q function can be computed using \R functions. An implementation is shown in this note.
\preformatted{
See NEWS file and entries for version 2.0.1 for this "R Marcum"
"marcumq" <- function(a, b, nu=1) {
      pchisq(b^2, df=2*nu, ncp=a^2, lower.tail=FALSE) }
}
}
\author{ W.H. Asquith}
\references{
Asquith, W.H., 2011, Distributional analysis with L-moment statistics using the R environment for statistical computing: Createspace Independent Publishing Platform, ISBN 978--146350841--8.
}
\seealso{\code{\link{pdfrice}}, \code{\link{quarice}}, \code{\link{lmomrice}}, \code{\link{parrice}}}
\examples{
lmr <- vec2lmom(c(45,0.27), lscale=FALSE)
cdfrice(35,parrice(lmr))
}
\keyword{distribution}
\keyword{cumulative distribution function}
\keyword{Distribution: Rice}
