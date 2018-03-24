\name{parkap}
\alias{parkap}
\title{Estimate the Parameters of the Kappa Distribution}
\description{
This function estimates the parameters of the Kappa distribution given the L-moments of the data in an L-moment object such as that returned by \code{\link{lmoms}}. The relations between distribution parameters and L-moments are seen under \code{\link{lmomkap}}, but of relevance to this documentation, the upper bounds of L-kurtosis (\eqn{\tau_4}) and a function of L-skew (\eqn{\tau_3}) is given by
\deqn{\tau_4 < \frac{5\tau_3^2+1}{6}}
This bounds is equal to the Generalized Logistic distribution (\code{\link{parglo}}) and failure occurs if this upper bounds is exceeded. However, the argument \code{snap.tau4}, if set, will set \eqn{\tau_4} equal to the upper bounds of \eqn{\tau_4} of the distribution to the relation above. This value of \eqn{\tau_4} should be close enough numerically The argument \code{nudge.tau4} is provided to offset \eqn{\tau_4} downward just a little. This keeps the relation operator as \dQuote{\eqn{<}} in the bounds above to match Hosking's tradition as his sources declare \dQuote{\eqn{\ge}} as above the GLO. The nudge here hence is not zero, which is a little different compared to the conceptually similar snapping in \code{\link{paraep4}}.
}
\usage{
parkap(lmom, checklmom=TRUE,
             snap.tau4=FALSE, nudge.tau4=sqrt(.Machine$double.eps), ...)
}
\arguments{
  \item{lmom}{An L-moment object created by \code{\link{lmoms}} or \code{\link{vec2lmom}}.}
  \item{checklmom}{Should the \code{lmom} be checked for validity using the \code{\link{are.lmom.valid}} function. Normally this should be left as the default and it is very unlikely that the L-moments will not be viable (particularly in the \eqn{\tau_4} and \eqn{\tau_3} inequality). However, for some circumstances or large simulation exercises then one might want to bypass this check.}
  \item{snap.tau4}{A logical to \dQuote{snap} the \eqn{\tau_4} downwards to the lower boundary if the given \eqn{\tau_4} is greater than the boundary described as above.}
  \item{nudge.tau4}{An offset to the snapping of \eqn{\tau_4} intended to move \eqn{\tau_4} just below the upper bounds. (The absolute value of the nudge is made internally to ensure only downward adjustment by a subtraction operation.)}
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

\dontrun{
parkap(vec2lmom(c(0,1,.3,.8)), snap.tau4=TRUE) # Tau=0.8 is way above the GLO.}
}
\keyword{distribution (parameters)}
\keyword{Distribution: Kappa}
