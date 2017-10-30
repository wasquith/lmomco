\name{cdfgov}
\alias{cdfgov}
\title{Cumulative Distribution Function of the Govindarajulu Distribution}
\description{
This function computes the cumulative probability or nonexceedance probability of the Govindarajulu distribution given parameters (\eqn{\xi}, \eqn{\alpha}, and \eqn{\beta}) computed by \code{\link{pargov}}. The cumulative distribution function has no explicit form and requires numerical methods. The \R function \code{uniroot} is used to root the quantile function \code{\link{quagov}} to compute the nonexceedance probability. The function returns 0 or 1 if the \code{x} argument is at or beyond the limits of the distribution as specified by the parameters.
}
\usage{
cdfgov(x, para)
}
\arguments{
  \item{x}{A real value vector.}
  \item{para}{The parameters from \code{\link{pargov}} or  \code{\link{vec2par}}.}
}

\value{
  Nonexceedance probability (\eqn{F}) for \eqn{x}.
}
\references{
Gilchrist, W.G., 2000, Statistical modelling with quantile functions: Chapman and Hall/CRC, Boca Raton.

Nair, N.U., Sankaran, P.G., and Balakrishnan, N., 2013, Quantile-based reliability analysis: Springer, New York.

Nair, N.U., Sankaran, P.G., and Vineshkumar, B., 2012, The Govindarajulu distribution---Some Properties and applications: Communications in Statistics, Theory and Methods, v. 41, no. 24, pp. 4391--4406.
}
\author{ W.H. Asquith}
\seealso{\code{\link{pdfgov}}, \code{\link{quagov}}, \code{\link{lmomgov}}, \code{\link{pargov}}}
\examples{
  lmr <- lmoms(c(123,34,4,654,37,78))
  cdfgov(50,pargov(lmr))
}
\keyword{distribution}
\keyword{cumulative distribution function}
\keyword{Distribution: Govindarajulu}
