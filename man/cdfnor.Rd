\name{cdfnor}
\alias{cdfnor}
\title{Cumulative Distribution Function of the Normal Distribution}
\description{
This function computes the cumulative probability or nonexceedance probability of the Normal distribution given parameters of the distribution computed by \code{\link{parnor}}.  The cumulative distribution function is
\deqn{F(x) = \Phi((x-\mu)/\sigma) \mbox{,}}
where \eqn{F(x)} is the nonexceedance probability for quantile \eqn{x}, \eqn{\mu} is the arithmetic mean, and \eqn{\sigma} is the standard deviation, and \eqn{\Phi} is the cumulative distribution function of the Standard Normal
distribution, and thus the \R function \code{pnorm} is used.
}
\usage{
cdfnor(x, para)
}
\arguments{
  \item{x}{A real value vector.}
  \item{para}{The parameters from \code{\link{parnor}} or  \code{\link{vec2par}}.}
}

\value{
  Nonexceedance probability (\eqn{F}) for \eqn{x}.
}
\references{
Hosking, J.R.M., 1990, L-moments---Analysis and estimation of
distributions using linear combinations of order statistics: Journal
of the Royal Statistical Society, Series B, v. 52, pp. 105--124.

Hosking, J.R.M., and Wallis, J.R., 1997, Regional frequency analysis---An
approach based on L-moments: Cambridge University Press.
}
\author{W.H. Asquith}
\seealso{\code{\link{pdfnor}}, \code{\link{quanor}}, \code{\link{lmomnor}}, \code{\link{parnor}}}
\examples{
  lmr <- lmoms(c(123,34,4,654,37,78))
  cdfnor(50,parnor(lmr))
}
\keyword{distribution}
\keyword{cumulative distribution function}
\keyword{Distribution: Normal}
