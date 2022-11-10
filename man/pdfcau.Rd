\name{pdfcau}
\alias{pdfcau}
\title{Probability Density Function of the Cauchy Distribution}
\description{
This function computes the probability density
of the Cauchy distribution given parameters (\eqn{\xi} and \eqn{\alpha})  provided by \code{\link{parcau}}.  The probability density function is
\deqn{f(x) = \left(\pi \alpha \left[1 + \left({\frac{x-\xi}{\alpha}}\right)^2\right] \right)^{-1}\mbox{,}}
where \eqn{f(x)} is the probability density for quantile \eqn{x},
\eqn{\xi} is a location parameter, and \eqn{\alpha} is a scale parameter.
}
\usage{
pdfcau(x, para)
}
\arguments{
  \item{x}{A real value vector.}
  \item{para}{The parameters from \code{\link{parcau}} or \code{\link{vec2par}}.}
}

\value{
  Probability density (\eqn{f}) for \eqn{x}.
}
\references{
Elamir, E.A.H., and Seheult, A.H., 2003, Trimmed L-moments: Computational Statistics
and Data Analysis, vol. 43, pp. 299--314.

Evans, Merran, Hastings, Nicholas, Peacock, J.B., 2000, Statistical distributions: 3rd ed., Wiley, New York.

Gilchrist, W.G., 2000, Statistical modeling with quantile functions:
Chapman and Hall/CRC, Boca Raton, FL.
}
\author{W.H. Asquith}
\seealso{\code{\link{cdfcau}}, \code{\link{quacau}}, \code{\link{lmomcau}}, \code{\link{parcau}}, \code{\link{vec2par}}}
\examples{
  cau <- vec2par(c(12,12),type='cau')
  x <- quacau(0.5,cau)
  pdfcau(x,cau)
}
\keyword{distribution}
\keyword{probability density function}
\keyword{Distribution: Cauchy}
