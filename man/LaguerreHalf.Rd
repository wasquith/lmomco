\name{LaguerreHalf}
\alias{LaguerreHalf}
\title{Laguerre Polynomial (Half)}
\description{
This function computes the Laguerre polynomial, which is useful in applications involving the variance of the Rice distribution (see \code{\link{parrice}}). The Laguerre polynomial is
\deqn{
L_{1/2}(x) = \exp^{x/2}\times[(1-x)I_0(-x/2) - xI_1(-x/2)]\mbox{,}
}
where the modified Bessel function of the first kind is \eqn{I_k(x)}, which has an \R implementation in \code{besselI}, and for strictly integer \eqn{k} is defined as
\deqn{I_k(x) = \frac{1}{\pi} \int_0^\pi \exp(x\cos(\theta)) \cos(k \theta)\; \mathrm{d}\theta\mbox{.}}
}
\usage{
LaguerreHalf(x)
}
\arguments{
  \item{x}{A value.}
}
\value{
The value for the Laguerre polynomial is returned.
}
\author{ W.H. Asquith}
\seealso{\code{\link{pdfrice}}}
\examples{
LaguerreHalf(-100^2/(2*10^2))
}
\keyword{special function}
