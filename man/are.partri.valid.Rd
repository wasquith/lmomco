\name{are.partri.valid}
\alias{are.partri.valid}
\title{Are the Distribution Parameters Consistent with the Asymmetric Triangular Distribution}
\description{
Is the distribution parameter object consistent with the corresponding distribution?  The distribution functions  (\code{\link{cdftri}}, \code{\link{pdftri}}, \code{\link{quatri}}, and
\code{\link{lmomtri}}) require consistent parameters to return the cumulative
probability (nonexceedance), density, quantile, and L-moments of the distribution,
respectively.  These functions internally use the \code{\link{are.partri.valid}}
function.
}
\usage{
are.partri.valid(para, nowarn=FALSE)
}
\arguments{
  \item{para}{A distribution parameter list returned by \code{\link{partri}} or \code{\link{vec2par}}.}
  \item{nowarn}{A logical switch on warning suppression. If \code{TRUE} then \code{options(warn=-1)} is made and restored on return. This switch is to permit calls in which warnings are not desired as the user knows how to handle the returned value---say in an optimization algorithm.}
}
\value{
  \item{TRUE}{If the parameters are \code{tri} consistent.}
  \item{FALSE}{If the parameters are not \code{tri} consistent.}
}
\note{
This function calls \code{\link{is.tri}} to verify consistency between
the distribution parameter object and the intent of the user.
}
\author{ W.H. Asquith}
\seealso{ \code{\link{is.tri}}, \code{\link{partri}} }
\examples{
para <- partri(lmoms(c(46, 70, 59, 36, 71, 48, 46, 63, 35, 52)))
if(are.partri.valid(para)) Q <- quatri(0.5,para)
}
\keyword{utility (distribution)}
\keyword{Distribution: Asymmetric Triangular}
\keyword{Distribution: Triangular}
\keyword{utility (distribution, parameter validation)}

