\name{quaemu}
\alias{quaemu}
\title{Quantile Function of the Eta-Mu Distribution}
\description{
This function computes the quantiles of the Eta-Mu (\eqn{\eta:\mu}) distribution given \eqn{\eta} and \eqn{\mu}) computed by \code{\link{paremu}}. The quantile function is complex and numerical rooting of the cumulative distribution function (\code{\link{cdfemu}}) is used.
}
\usage{
quaemu(f, para, paracheck=TRUE, yacoubsintegral=TRUE, eps=1e-7)
}
\arguments{
  \item{f}{Nonexceedance probability (\eqn{0 \le F \le 1}).}
  \item{para}{The parameters from \code{\link{paremu}} or \code{\link{vec2par}}.}
  \item{paracheck}{A logical controlling whether the parameters are checked for validity. Overriding of this check might be extremely important and needed for use of the quantile function in the context of TL-moments with nonzero trimming.}
  \item{yacoubsintegral}{A logical controlling whether the integral by Yacoub (2007) is used for the cumulative distribution function instead of numerical integration of \code{\link{pdfemu}}.}
  \item{eps}{A close-enough error term for the recursion process.}
}
\value{
  Quantile value for nonexceedance probability \eqn{F}.
}
\references{
Yacoub, M.D., 2007, The kappa-mu distribution and the eta-mu distribution: IEEE Antennas and Propagation Magazine, v. 49, no. 1, pp. 68--81
}
\author{ W.H. Asquith}
\seealso{\code{\link{cdfemu}}, \code{\link{pdfemu}}, \code{\link{lmomemu}},  \code{\link{paremu}}}
\examples{
\dontrun{
quaemu(0.75,vec2par(c(0.9, 1.5), type="emu")) #}
}
\keyword{distribution}
\keyword{quantile function}
\keyword{Distribution: Eta-Mu}
