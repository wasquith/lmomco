\name{are.parsmd.valid}
\alias{are.parsmd.valid}
\title{Are the Distribution Parameters Consistent with the Singh--Maddala Distribution}
\description{
Is the distribution parameter object consistent with the corresponding distribution?  The distribution functions (\code{\link{cdfsmd}}, \code{\link{pdfsmd}}, \code{\link{quasmd}}, and \code{\link{lmomsmd}}) require consistent parameters to return the cumulative probability (nonexceedance), density, quantile, and L-moments of the distribution, respectively. These functions internally use the \code{are.parsmd.valid} function. The parameter constraints are simple \eqn{a > 0} (scale), \eqn{b > 0} (shape), and \eqn{q > 0} (shape).
}
\usage{
are.parsmd.valid(para, nowarn=FALSE)
}
\arguments{
  \item{para}{A distribution parameter list returned by \code{\link{parsmd}} or \code{\link{vec2par}}.}
  \item{nowarn}{A logical switch on warning suppression. If \code{TRUE} then \code{options(warn=-1)} is made and restored on return. This switch is to permit calls in which warnings are not desired as the user knows how to handle the returned value---say in an optimization algorithm.}
}
\value{
  \item{TRUE}{If the parameters are \code{smd} consistent.}
  \item{FALSE}{If the parameters are not \code{smd} consistent.}
}
\references{
Shahzad, M.N., and Zahid, A., 2013, Parameter estimation of Singh Maddala distribution by moments: International Journal of Advanced Statistics and Probability, v. 1, no. 3, pp. 121--131, \doi{10.14419/ijasp.v1i3.1206}.
}
\note{
This function calls \code{\link{is.smd}} to verify consistency between
the distribution parameter object and the intent of the user.
}
\author{W.H. Asquith}
\seealso{\code{\link{is.smd}}, \code{\link{parsmd}} }
\examples{
#para <- parsmd(lmoms(c(123, 34, 4, 654, 37, 78)))
#if(are.parsmd.valid(para)) Q <- quasmd(0.5, para)
}
\keyword{utility (distribution)}
\keyword{Distribution: Singh-Maddala}
\keyword{Distribution: Burr Type XII}
\keyword{utility (distribution/parameter validation)}

