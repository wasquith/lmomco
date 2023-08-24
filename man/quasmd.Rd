\name{quasmd}
\alias{quasmd}
\title{Quantile Function of the Singh--Maddala Distribution}
\description{
This function computes the quantiles of the Singh--Maddala (Burr Type XII) distribution given parameters (\eqn{\xi}, \eqn{a}, \eqn{b}, and \eqn{q}) computed by \code{\link{parsmd}}. The quantile function  is

\deqn{x(F) = \xi + a\biggl((1-F)^{-1/q} - 1 \biggr)^{1/b}\mbox{,}}

where \eqn{x(F)} with \eqn{0 \le x \le \infty} is the quantile for nonexceedance probability \eqn{F}, \eqn{\xi} is a location parameter, \eqn{a} is a scale parameter (\eqn{a > 0}), \eqn{b} is a shape parameter (\eqn{b > 0}), and \eqn{q} is another shape parameter (\eqn{q > 0}).
}
\usage{
quasmd(f, para, paracheck=TRUE)
}
\arguments{
  \item{f}{Nonexceedance probability (\eqn{0 \le F \le 1}).}
  \item{para}{The parameters from \code{\link{parsmd}} or \code{\link{vec2par}}.}
  \item{paracheck}{A logical controlling whether the parameters are checked for validity. Overriding of this check might be extremely important and needed for use of the quantile function in the context of TL-moments with nonzero trimming.}
}
\value{
  Quantile value for nonexceedance probability \eqn{F}.
}
\references{
Kumar, D., 2017, The Singh--Maddala distribution---Properties and estimation: International Journal of System Assurance Engineering and Management, v. 8, no. S2, 15 p., \doi{10.1007/s13198-017-0600-1}.

Shahzad, M.N., and Zahid, A., 2013, Parameter estimation of Singh Maddala distribution by moments: International Journal of Advanced Statistics and Probability, v. 1, no. 3, pp. 121--131, \doi{10.14419/ijasp.v1i3.1206}.
}
\author{W.H. Asquith}
\seealso{\code{\link{cdfsmd}}, \code{\link{pdfsmd}}, \code{\link{lmomsmd}}, \code{\link{parsmd}}}
\examples{
quasmd(0.99, parsmd(vec2lmom(c(155, 118.6, 0.6, 0.45)))) # 1547.337 99th percentile
}
\keyword{distribution}
\keyword{quantile function}
\keyword{Distribution: Singh-Maddala}
\keyword{Distribution: Burr Type XII}
