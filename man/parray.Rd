\name{parray}
\alias{parray}
\title{Estimate the Parameters of the Rayleigh Distribution}
\description{
This function estimates the parameters of the Rayleigh distribution given the L-moments of the data in an L-moment object such as that returned by \code{\link{lmoms}}. The relations between distribution parameters and L-moments are
\deqn{\alpha = \frac{2\lambda_2\sqrt{\pi}}{\sqrt{2} - 1}\mbox{,}}
and
\deqn{\xi = \lambda_1 - \alpha\sqrt{\pi/2}\mbox{.}}
}
\usage{
parray(lmom, xi=NULL, checklmom=TRUE, ...)
}
\arguments{
  \item{lmom}{An L-moment object created by \code{\link{lmoms}} or \code{\link{vec2lmom}}.}
  \item{xi}{The lower limit of the distribution. If \eqn{\xi} is known then alternative algorithms are triggered and only the first L-moment is required for fitting.}
  \item{checklmom}{Should the \code{lmom} be checked for validity using the \code{\link{are.lmom.valid}
} function. Normally this should be left as the default and it is very unlikely that the L-moments will not be viable (particularly in the \eqn{\tau_4} and \eqn{\tau_3} inequality). However, for some circumstances or large simulation exercises then one might want to bypass this check.}
  \item{...}{Other arguments to pass.}
}

\value{
  An \R \code{list} is returned.

  \item{type}{The type of distribution: \code{ray}.}
  \item{para}{The parameters of the distribution.}
  \item{source}{The source of the parameters: \dQuote{parray}.}
}
\references{
Hosking, J.R.M., 1986, The theory of probability weighted moments: Research Report RC12210, IBM Research Division, Yorkton Heights, N.Y.
}
\author{W.H. Asquith}
\seealso{\code{\link{lmomray}},
\code{\link{cdfray}}, \code{\link{pdfray}}, \code{\link{quaray}}
}
\examples{
lmr <- lmoms(rnorm(20))
parray(lmr)
}
\keyword{distribution (parameters)}
\keyword{Distribution: Rayleigh}
