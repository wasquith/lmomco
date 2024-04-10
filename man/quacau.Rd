\name{quacau}
\alias{quacau}
\title{Quantile Function of the Cauchy Distribution}
\description{
This function computes the quantiles of the Cauchy distribution given parameters (\eqn{\xi} and \eqn{\alpha}) of the distribution provided by \code{\link{parcau}}. The quantile function of the distribution is
\deqn{x(F) = \xi + \alpha \times \tan\bigl(\pi(F-0.5)\bigr) \mbox{,}}
where \eqn{x(F)} is the quantile for nonexceedance probability \eqn{F}, \eqn{\xi} is a location parameter and \eqn{\alpha} is a scale parameter. The quantile function of the Cauchy distribution is supported by \R function \code{qcauchy}. This function does not use \code{qcauchy} because \code{qcauchy} does not return \code{Inf} for \eqn{F = 1} although it
returns \code{-Inf} for \eqn{F = 0}.
}
\usage{
quacau(f, para, paracheck=TRUE)
}
\arguments{
  \item{f}{Nonexceedance probability (\eqn{0 \le F \le 1}).}
  \item{para}{The parameters from \code{\link{parcau}} or \code{\link{vec2par}}.}
  \item{paracheck}{A logical controlling whether the parameters are checked for validity. Overriding of this check might be extremely important and needed for use of the distribution quantile function in the context of TL-moments with nonzero trimming.}
}
\value{
  Quantile value for for nonexceedance probability \eqn{F}.
}
\references{
Elamir, E.A.H., and Seheult, A.H., 2003, Trimmed L-moments: Computational Statistics
and Data Analysis, v. 43, pp. 299--314.

Gilchirst, W.G., 2000, Statistical modeling with quantile functions:
Chapman and Hall/CRC, Boca Raton, FL.
}
\author{W.H. Asquith}
\seealso{\code{\link{cdfcau}}, \code{\link{pdfcau}}, \code{\link{lmomcau}}, \code{\link{parcau}}
}
\examples{
  para <- c(12,12)
  quacau(.5,vec2par(para,type='cau'))
}
\keyword{distribution}
\keyword{quantile function}
\keyword{Distribution: Cauchy}
