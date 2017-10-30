\name{are.parpe3.valid}
\alias{are.parpe3.valid}
\title{Are the Distribution Parameters Consistent with the Pearson Type III Distribution}
\description{
Is the distribution parameter object consistent with the corresponding distribution?  The distribution functions  (\code{\link{cdfpe3}}, \code{\link{pdfpe3}}, \code{\link{quape3}}, and
\code{\link{lmompe3}}) require consistent parameters to return the cumulative
probability (nonexceedance), density, quantile, and L-moments of the distribution,
respectively.  These functions internally use the \code{\link{are.parpe3.valid}}
function.
}
\usage{
are.parpe3.valid(para, nowarn=FALSE)
}
\arguments{
  \item{para}{A distribution parameter list returned by \code{\link{parpe3}} or \code{\link{vec2par}}.}
  \item{nowarn}{A logical switch on warning suppression. If \code{TRUE} then \code{options(warn=-1)} is made and restored on return. This switch is to permit calls in which warnings are not desired as the user knows how to handle the returned value---say in an optimization algorithm.}
}
\value{
  \item{TRUE}{If the parameters are \code{pe3} consistent.}
  \item{FALSE}{If the parameters are not \code{pe3} consistent.}
}
\references{
Hosking, J.R.M., 1996, FORTRAN routines for use with the method of L-moments:
Version 3, IBM Research Report RC20525, T.J. Watson Research Center,
Yorktown Heights, New York.

Hosking, J.R.M., and Wallis, J.R., 1997, Regional frequency analysis---An
approach based on L-moments: Cambridge University Press.
}
\note{
This function calls \code{\link{is.pe3}} to verify consistency between
the distribution parameter object and the intent of the user.
}
\author{ W.H. Asquith}
\seealso{ \code{\link{is.pe3}}, \code{\link{parpe3}} }
\examples{
para <- parpe3(lmoms(c(123,34,4,654,37,78)))
if(are.parpe3.valid(para)) Q <- quape3(0.5,para)
}
\keyword{utility (distribution)}
\keyword{Distribution: Pearson Type III}
\keyword{utility (distribution, parameter validation)}

