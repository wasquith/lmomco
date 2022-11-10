\name{are.pargno.valid}
\alias{are.pargno.valid}
\title{Are the Distribution Parameters Consistent with the Generalized Normal Distribution}
\description{
Is the distribution parameter object consistent with the corresponding distribution?  The distribution functions  (\code{\link{cdfgno}}, \code{\link{pdfgno}}, \code{\link{quagno}}, and
\code{\link{lmomgno}}) require consistent parameters to return the cumulative
probability (nonexceedance), density, quantile, and L-moments of the distribution,
respectively.  These functions internally use the \code{\link{are.pargno.valid}}
function.
}
\usage{
are.pargno.valid(para, nowarn=FALSE)
}
\arguments{
  \item{para}{A distribution parameter list returned by \code{\link{pargno}} or \code{\link{vec2par}}.}
  \item{nowarn}{A logical switch on warning suppression. If \code{TRUE} then \code{options(warn=-1)} is made and restored on return. This switch is to permit calls in which warnings are not desired as the user knows how to handle the returned value---say in an optimization algorithm.}
}
\value{
  \item{TRUE}{If the parameters are \code{gno} consistent.}
  \item{FALSE}{If the parameters are not \code{gno} consistent.}
}
\references{
Hosking, J.R.M., 1996, FORTRAN routines for use with the method of L-moments:
Version 3, IBM Research Report RC20525, T.J. Watson Research Center,
Yorktown Heights, New York.

Hosking, J.R.M., and Wallis, J.R., 1997, Regional frequency analysis---An
approach based on L-moments: Cambridge University Press.
}
\note{
This function calls \code{\link{is.gno}} to verify consistency between
the distribution parameter object and the intent of the user.
}
\author{W.H. Asquith}
\seealso{\code{\link{is.gno}}, \code{\link{pargno}}
}
\examples{
para <- pargno(lmoms(c(123,34,4,654,37,78)))
if(are.pargno.valid(para)) Q <- quagno(0.5,para)
}
\keyword{utility (distribution)}
\keyword{Distribution: Generalized Normal}
\keyword{utility (distribution/parameter validation)}

