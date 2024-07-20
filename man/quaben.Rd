\encoding{utf8}
\name{quaben}
\alias{quaben}
\title{Quantile Function of the Benford Distribution}
\description{
This function computes the quantiles of the Benford distribution (Benford's Law) given parameter defining the number of first M-significant figures and the numeric base. The quantile function has no analytical form and summation of the probability mass function (to form the cumulative distribution function, see also \code{\link{cdfben}}) is used with clever use of the \code{cut()} function.
}
\usage{
quaben(f, para=list(para=c(1, 10)), ...)
}
\arguments{
  \item{f}{Nonexceedance probability (\eqn{0 \le F \le 1}).}
  \item{para}{The number of the first M-significant digits followed by the numerical base (only base10 supported) and the list structure mimics similar uses of the \pkg{lmomco} list structure. Default are the first significant digits and hence the digits 1 through 9.}
  \item{...}{Additional arguments to pass (not likely to be needed but changes in base handling might need this).}
}
\value{
  Quantile value for nonexceedance probability \eqn{F}.
}
\references{
Benford, F., 1938, The law of anomalous numbers: Proceedings of the American Philosophical Society, v. 78, no. 4, pp. 551--572, \url{https://www.jstor.org/stable/984802}.

Goodman, W., 2016, The promises and pitfalls of Benford’s law: Significance (Magazine), June 2015, pp. 38--41, \doi{10.1111/j.1740-9713.2016.00919.x}.
}
\author{W.H. Asquith}
\seealso{\code{\link{cdfben}}, \code{\link{pmfben}}}
\examples{
para <- list(para=c(1, 10))
quaben(    cdfben(  5, para=para) , para=para) # 5
quaben(sum(pmfben(1:5, para=para)), para=para) # 5
}
\keyword{distribution (discrete)}
\keyword{quantile function}
\keyword{Distribution: Benford}
\concept{Benford Law}
