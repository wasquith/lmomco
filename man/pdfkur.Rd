\name{pdfkur}
\alias{pdfkur}
\title{Probability Density Function of the Kumaraswamy Distribution}
\description{
This function computes the probability density
of the Kumaraswamy distribution given parameters (\eqn{\alpha} and \eqn{\beta})  computed by \code{\link{parkur}}. The probability density function is
\deqn{f(x) = \alpha\beta x^{\alpha - 1}(1-x^\alpha)^{\beta-1} \mbox{,}}
where \eqn{f(x)} is the nonexceedance probability for quantile \eqn{x},
\eqn{\alpha} is a shape parameter, and \eqn{\beta} is a shape parameter.
}
\usage{
pdfkur(x, para)
}
\arguments{
  \item{x}{A real value vector.}
  \item{para}{The parameters from \code{\link{parkur}} or \code{\link{vec2par}}.}
}
\value{
  Probability density (\eqn{f}) for \eqn{x}.
}
\references{
Jones, M.C., 2009, Kumaraswamy's distribution---A beta-type distribution with
some tractability advantages: Statistical Methodology, v. 6, pp. 70--81.
}
\author{W.H. Asquith}
\seealso{\code{\link{cdfkur}}, \code{\link{quakur}}, \code{\link{lmomkur}}, \code{\link{parkur}}}
\examples{
  lmr <- lmoms(c(0.25, 0.4, 0.6, 0.65, 0.67, 0.9))
  kur <- parkur(lmr)
  x <- quakur(0.5,kur)
  pdfkur(x,kur)
}
\keyword{distribution}
\keyword{probability density function}
\keyword{Distribution: Kumaraswamy}
