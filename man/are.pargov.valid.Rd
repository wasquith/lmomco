\name{are.pargov.valid}
\alias{are.pargov.valid}
\title{Are the Distribution Parameters Consistent with the Govindarajulu Distribution}
\description{
Is the distribution parameter object consistent with the corresponding distribution?  The distribution functions  (\code{\link{cdfgov}}, \code{\link{pdfgov}}, \code{\link{quagov}}, and
\code{\link{lmomgov}}) require consistent parameters to return the cumulative
probability (nonexceedance), density, quantile, and L-moments of the distribution,
respectively.  These functions internally use the \code{\link{are.pargov.valid}}
function.
}
\usage{
are.pargov.valid(para, nowarn=FALSE)
}
\arguments{
  \item{para}{A distribution parameter list returned by \code{\link{pargov}} or \code{\link{vec2par}}.}
  \item{nowarn}{A logical switch on warning suppression. If \code{TRUE} then \code{options(warn=-1)} is made and restored on return. This switch is to permit calls in which warnings are not desired as the user knows how to handle the returned value---say in an optimization algorithm.}
}
\value{
  \item{TRUE}{If the parameters are \code{gov} consistent.}
  \item{FALSE}{If the parameters are not \code{gov} consistent.}
}
\references{
Gilchrist, W.G., 2000, Statistical modelling with quantile functions: Chapman and Hall/CRC, Boca Raton.

Nair, N.U., Sankaran, P.G., and Balakrishnan, N., 2013, Quantile-based reliability analysis: Springer, New York.
}
\note{
This function calls \code{\link{is.gov}} to verify consistency between
the distribution parameter object and the intent of the user.
}
\author{ W.H. Asquith}
\seealso{ \code{\link{is.gov}}, \code{\link{pargov}} }
\examples{
para <- pargov(lmoms(c(123,34,4,654,37,78)))
if(are.pargov.valid(para)) Q <- quagov(0.5,para)
}
\keyword{utility (distribution)}
\keyword{Distribution: Govindarajulu}
\keyword{utility (distribution/parameter validation)}


