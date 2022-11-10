\name{quakmu}
\alias{quakmu}
\title{Quantile Function of the Kappa-Mu Distribution}
\description{
This function computes the quantiles of the Kappa-Mu (\eqn{\kappa:\mu}) distribution given parameters (\eqn{\kappa} and \eqn{\alpha})  computed by \code{\link{parkmu}}. The quantile function is complex and numerical rooting of the cumulative distribution function (\code{\link{cdfkmu}}) is used.
}
\usage{
quakmu(f, para, paracheck=TRUE, getmed=FALSE, qualo=NA, quahi=NA, verbose=FALSE,
                marcumQ=TRUE, marcumQmethod=c("chisq", "delta", "integral"))
}
\arguments{
  \item{f}{Nonexceedance probability (\eqn{0 \le F \le 1}).}
  \item{para}{The parameters from \code{\link{parkmu}} or \code{\link{vec2par}}.}
  \item{paracheck}{A logical controlling whether the parameters are checked for validity. Overriding of this check might be extremely important and needed for use of the quantile function in the context of TL-moments with nonzero trimming.}
  \item{getmed}{Same argument for \code{\link{cdfkmu}}. Because of nesting a \code{quakmu} call in \code{\link{cdfkmu}}, this argument and the next two are shown here are to avoid confusion in use of \code{...} instead. This argument should not overrided by the user.}
    \item{qualo}{A lower limit of the range of \eqn{x} to look for a \code{uniroot} of \eqn{F(x)}.}
    \item{quahi}{An upper limit of the range of \eqn{x} to look for a \code{uniroot} of \eqn{F(x)}.}
  \item{verbose}{Should alert messages be shown by \code{message()}?}
  \item{marcumQ}{Same argument for \code{\link{cdfkmu}}, which the user can set change.}
  \item{marcumQmethod}{Same argument for \code{\link{cdfkmu}}, which the user can set change.}
}
\value{
  Quantile value for nonexceedance probability \eqn{F}.
}
\references{
Yacoub, M.D., 2007, The kappa-mu distribution and the eta-mu distribution: IEEE Antennas and Propagation Magazine, v. 49, no. 1, pp. 68--81
}
\author{W.H. Asquith}
\seealso{\code{\link{cdfkmu}}, \code{\link{pdfkmu}}, \code{\link{lmomkmu}}, \code{\link{parkmu}}}
\examples{
quakmu(0.75,vec2par(c(0.9, 1.5), type="kmu"))
}
\keyword{distribution}
\keyword{quantile function}
\keyword{Distribution: Kappa-Mu}
