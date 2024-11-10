\encoding{utf8}
\name{cdfben}
\alias{cdfben}
\title{Cumulative Distribution Function of the Benford Distribution}
\description{
This function computes the cumulative probability or nonexceedance probability of the Benford distribution (Benford's Law) given parameters defining the number of first M-significant digits and the numeric base. The cumulative distribution function has a somewhat simple analytical form by direct summation of the probability mass function (\code{\link{pmfben}}).
}
\usage{
cdfben(d, para=list(para=c(1, 10)), ...)
}
\arguments{
  \item{d}{A integer value vector of M-significant digits.}
  \item{para}{The number of the first M-significant digits followed by the numerical base (only base10 supported) and the list structure mimics similar uses of the \pkg{lmomco} list structure. Default are the first significant digits and hence the digits 1 through 9.}
  \item{...}{Additional arguments to pass (not likely to be needed but changes in base handling might need this).}
}

\value{
  Nonexceedance probability (\eqn{F}) for \eqn{x}.
}
\references{
Benford, F., 1938, The law of anomalous numbers: Proceedings of the American Philosophical Society, v. 78, no. 4, pp. 551--572, \url{https://www.jstor.org/stable/984802}.

Goodman, W., 2016, The promises and pitfalls of Benford’s law: Significance (Magazine), June 2015, pp. 38--41, \doi{10.1111/j.1740-9713.2016.00919.x}.
}
\author{W.H. Asquith}
\seealso{\code{\link{pmfben}}, \code{\link{quaben}}}
\examples{
para <- list(para=c(2, 10))
cdfben(c(15, 25), para=para) # 0.2041200 0.4149733

sum(diff(cdfben(seq(10,99,0.1), para=para))) + cdfben(10, para=para) # 1
}
\keyword{distribution (discrete)}
\keyword{cumulative distribution function}
\keyword{Distribution: Benford}
\concept{Benford Law}

