\name{cdfkur}
\alias{cdfkur}
\title{Cumulative Distribution Function of the Kumaraswamy Distribution}
\description{
This function computes the cumulative probability or nonexceedance probability of the Kumaraswamy distribution given parameters (\eqn{\alpha} and \eqn{\beta}) computed by \code{\link{parkur}}. The cumulative distribution function is
\deqn{F(x) = 1 - (1-x^\alpha)^\beta \mbox{,}}
where \eqn{F(x)} is the nonexceedance probability for quantile \eqn{x},
\eqn{\alpha} is a shape parameter, and \eqn{\beta} is a shape parameter.
}
\usage{
cdfkur(x, para)
}
\arguments{
  \item{x}{A real value vector.}
  \item{para}{The parameters from \code{\link{parkur}} or  \code{\link{vec2par}}.}
}

\value{
  Nonexceedance probability (\eqn{F}) for \eqn{x}.
}
\references{
Jones, M.C., 2009, Kumaraswamy's distribution---A beta-type distribution with
some tractability advantages: Statistical Methodology, v. 6, pp. 70--81.
}
\author{ W.H. Asquith}
\seealso{\code{\link{pdfkur}}, \code{\link{quakur}}, \code{\link{lmomkur}}, \code{\link{parkur}}}
\examples{
  lmr <- lmoms(c(0.25, 0.4, 0.6, 0.65, 0.67, 0.9))
  cdfkur(0.5,parkur(lmr))
}
\keyword{distribution}
\keyword{cumulative distribution function}
\keyword{Distribution: Kumaraswamy}
