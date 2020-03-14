\name{are.parlap.valid}
\alias{are.parlap.valid}
\title{Are the Distribution Parameters Consistent with the Laplace Distribution}
\description{
Is the distribution parameter object consistent with the corresponding distribution?  The distribution functions  (\code{\link{cdflap}}, \code{\link{pdflap}}, \code{\link{qualap}}, and \code{\link{lmomlap}}) require consistent parameters to return the cumulative probability (nonexceedance), density, quantile, and L-moments of the distribution, respectively.  These functions internally use the \code{\link{are.parlap.valid}} function.
}
\usage{
are.parlap.valid(para, nowarn=FALSE)
}
\arguments{
  \item{para}{A distribution parameter list returned by \code{\link{parlap}} or \code{\link{vec2par}}.}
  \item{nowarn}{A logical switch on warning suppression. If \code{TRUE} then \code{options(warn=-1)} is made and restored on return. This switch is to permit calls in which warnings are not desired as the user knows how to handle the returned value---say in an optimization algorithm.}
}
\value{
  \item{TRUE}{If the parameters are \code{lap} consistent.}
  \item{FALSE}{If the parameters are not \code{lap} consistent.}
}
\references{
Hosking, J.R.M., 1986, The theory of probability weighted moments: IBM Research Report RC12210, T.J. Watson Research Center, Yorktown Heights, New York.
}
\note{
This function calls \code{\link{is.lap}} to verify consistency between
the distribution parameter object and the intent of the user.
}
\author{ W.H. Asquith}
\seealso{ \code{\link{is.lap}}, \code{\link{parlap}} }
\examples{
para <- parlap(lmoms(c(123,34,4,654,37,78)))
if(are.parlap.valid(para)) Q <- qualap(0.5,para)
}
\keyword{utility (distribution)}
\keyword{Distribution: Laplace}
\keyword{utility (distribution/parameter validation)}

