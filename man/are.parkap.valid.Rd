\name{are.parkap.valid}
\alias{are.parkap.valid}
\title{Are the Distribution Parameters Consistent with the Kappa Distribution}
\description{
Is the distribution parameter object consistent with the corresponding distribution?  The distribution functions  (\code{\link{cdfkap}}, \code{\link{pdfkap}}, \code{\link{quakap}}, and
\code{\link{lmomkap}}) require consistent parameters to return the cumulative
probability (nonexceedance), density, quantile, and L-moments of the distribution,
respectively.  These functions internally use the \code{\link{are.parkap.valid}}
function.
}
\usage{
are.parkap.valid(para, nowarn=FALSE)
}
\arguments{
  \item{para}{A distribution parameter list returned by \code{\link{parkap}} or \code{\link{vec2par}}.}
  \item{nowarn}{A logical switch on warning suppression. If \code{TRUE} then \code{options(warn=-1)} is made and restored on return. This switch is to permit calls in which warnings are not desired as the user knows how to handle the returned value---say in an optimization algorithm.}
}
\value{
  \item{TRUE}{If the parameters are \code{kap} consistent.}
  \item{FALSE}{If the parameters are not \code{kap} consistent.}
}
\references{
Hosking, J.R.M., 1994, The four-parameter kappa distribution: IBM Journal of Reserach and Development, v. 38, no. 3, pp. 251--258.

Hosking, J.R.M., and Wallis, J.R., 1997, Regional frequency analysis---An approach based on L-moments: Cambridge University Press.
}
\note{
This function calls \code{\link{is.kap}} to verify consistency between
the distribution parameter object and the intent of the user.
}
\author{ W.H. Asquith}
\seealso{ \code{\link{is.kap}}, \code{\link{parkap}} }
\examples{
para <- parkap(lmoms(c(123,34,4,654,37,78)))
if(are.parkap.valid(para)) Q <- quakap(0.5,para)
}
\keyword{utility (distribution)}
\keyword{Distribution: Kappa}
\keyword{utility (distribution/parameter validation)}

