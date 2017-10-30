\name{parkap}
\alias{parkap}
\title{Estimate the Parameters of the Kappa Distribution}
\description{
This function estimates the parameters of the Kappa distribution given the L-moments of the data in an L-moment object such as that returned by \code{\link{lmoms}}. The relations between distribution parameters and L-moments are seen under \code{\link{lmomkap}}.
}
\usage{
parkap(lmom, checklmom=TRUE, ...)
}
\arguments{
  \item{lmom}{An L-moment object created by \code{\link{lmoms}} or \code{\link{vec2lmom}}.}
  \item{checklmom}{Should the \code{lmom} be checked for validity using the \code{\link{are.lmom.valid}} function. Normally this should be left as the default and it is very unlikely that the L-moments will not be viable (particularly in the \eqn{\tau_4} and \eqn{\tau_3} inequality). However, for some circumstances or large simulation exercises then one might want to bypass this check.}
  \item{...}{Other arguments to pass.}
}
\value{
  An \R \code{list} is returned.

  \item{type}{The type of distribution: \code{kap}.}
  \item{para}{The parameters of the distribution.}
  \item{source}{The source of the parameters: \dQuote{parkap}.}
  \item{support}{The support (or range) of the fitted distribution.}
  \item{ifail}{A numeric failure code.}
  \item{ifailtext}{A text message for the failure code.}
}
\references{
Hosking, J.R.M., 1994, The four-parameter kappa distribution: IBM Journal of Reserach and Development, v. 38, no. 3, pp. 251--258.

Hosking, J.R.M., and Wallis, J.R., 1997, Regional frequency analysis---An approach based on L-moments: Cambridge University Press.
}
\author{ W.H. Asquith}
\seealso{\code{\link{lmomkap}},
\code{\link{cdfkap}}, \code{\link{pdfkap}}, \code{\link{quakap}}
}
\examples{
lmr <- lmoms(rnorm(20))
parkap(lmr)
}
\keyword{distribution (parameters)}
\keyword{Distribution: Kappa}
