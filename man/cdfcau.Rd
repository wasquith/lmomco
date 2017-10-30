\name{cdfcau}
\alias{cdfcau}
\title{Cumulative Distribution Function of the Cauchy Distribution}
\description{
This function computes the cumulative probability or nonexceedance probability of the Cauchy distribution given parameters (\eqn{\xi} and \eqn{\alpha}) computed by \code{\link{parcau}}.  The cumulative distribution function is
\deqn{F(x) = \frac{\arctan(Y)}{\pi}+0.5 \mbox{,}}
where \eqn{Y} is
\deqn{Y = \frac{x - \xi}{\alpha}\mbox{, and}}
where \eqn{F(x)} is the nonexceedance probability for quantile \eqn{x},  \eqn{\xi} is a location parameter, and \eqn{\alpha} is a scale parameter.
}
\usage{
cdfcau(x, para)
}
\arguments{
  \item{x}{A real value vector.}
  \item{para}{The parameters from \code{\link{parcau}} or \code{\link{vec2par}}.}
}

\value{
  Nonexceedance probability (\eqn{F}) for \eqn{x}.
}
\references{
Elamir, E.A.H., and Seheult, A.H., 2003, Trimmed L-moments: Computational Statistics
and Data Analysis, v. 43, pp. 299--314.

Gilchrist, W.G., 2000, Statistical modeling with quantile functions:
Chapman and Hall/CRC, Boca Raton, FL.
}
\author{ W.H. Asquith}
\seealso{\code{\link{pdfcau}}, \code{\link{quacau}}, \code{\link{lmomcau}}, \code{\link{parcau}}
}
\examples{
  para <- c(12,12)
  cdfcau(50,vec2par(para,type='cau'))
}
\keyword{distribution}
\keyword{cumulative distribution function}
\keyword{Distribution: Cauchy}
