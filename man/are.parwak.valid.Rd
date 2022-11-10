\name{are.parwak.valid}
\alias{are.parwak.valid}
\title{Are the Distribution Parameters Consistent with the Wakeby Distribution}
\description{
Is the distribution parameter object consistent with the corresponding distribution?  The distribution functions  (\code{\link{cdfwak}}, \code{\link{pdfwak}}, \code{\link{quawak}}, and
\code{\link{lmomwak}}) require consistent parameters to return the cumulative
probability (nonexceedance), density, quantile, and L-moments of the distribution,
respectively.  These functions internally use the \code{\link{are.parwak.valid}}
function.
}
\usage{
are.parwak.valid(para, nowarn=FALSE)
}
\arguments{
  \item{para}{A distribution parameter list returned by \code{\link{parwak}} or \code{\link{vec2par}}.}
  \item{nowarn}{A logical switch on warning suppression. If \code{TRUE} then \code{options(warn=-1)} is made and restored on return. This switch is to permit calls in which warnings are not desired as the user knows how to handle the returned value---say in an optimization algorithm.}
}
\value{
  \item{TRUE}{If the parameters are \code{wak} consistent.}
  \item{FALSE}{If the parameters are not \code{wak} consistent.}
}
\references{
Hosking, J.R.M., 1996, FORTRAN routines for use with the method of L-moments:
Version 3, IBM Research Report RC20525, T.J. Watson Research Center,
Yorktown Heights, New York.

Hosking, J.R.M., and Wallis, J.R., 1997, Regional frequency analysis---An
approach based on L-moments: Cambridge University Press.
}
\note{
This function calls \code{\link{is.wak}} to verify consistency between
the distribution parameter object and the intent of the user.
}
\author{W.H. Asquith}
\seealso{\code{\link{is.wak}}, \code{\link{parwak}} }
\examples{
para <- parwak(lmoms(c(123,34,4,654,37,78)))
if(are.parwak.valid(para)) Q <- quawak(0.5,para)
}
\keyword{utility (distribution)}
\keyword{Distribution: Wakeby}
\keyword{utility (distribution/parameter validation)}

