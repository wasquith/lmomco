\name{are.paraep4.valid}
\alias{are.paraep4.valid}
\title{Are the Distribution Parameters Consistent with the 4-Parameter Asymmetric Exponential Power Distribution}
\description{
Is the distribution parameter object consistent with the corresponding distribution?  The distribution functions  (\code{\link{cdfaep4}}, \code{\link{pdfaep4}}, \code{\link{quaaep4}}, and
\code{\link{lmomaep4}}) require consistent parameters to return the cumulative
probability (nonexceedance), density, quantile, and L-moments of the distribution,
respectively.  These functions internally use the \code{\link{are.paraep4.valid}}
function.
}
\usage{
are.paraep4.valid(para, nowarn=FALSE)
}
\arguments{
  \item{para}{A distribution parameter list returned by \code{\link{paraep4}} or \code{\link{vec2par}}.}
  \item{nowarn}{A logical switch on warning suppression. If \code{TRUE} then \code{options(warn=-1)} is made and restored on return. This switch is to permit calls in which warnings are not desired as the user knows how to handle the returned value---say in an optimization algorithm.}
}
\value{
  \item{TRUE}{If the parameters are \code{aep4} consistent.}
  \item{FALSE}{If the parameters are not \code{aep4} consistent.}
}
\note{
This function calls \code{\link{is.aep4}} to verify consistency between
the distribution parameter object and the intent of the user.
}
\author{W.H. Asquith}
\references{
Asquith, W.H., 2014, Parameter estimation for the 4-parameter asymmetric exponential power distribution by the method of L-moments using R: Computational Statistics and Data Analysis, v. 71, pp. 955--970.

Delicado, P., and Goria, M.N., 2008, A small sample comparison of maximum likelihood,
moments and L-moments methods for the asymmetric exponential power distribution:
Computational Statistics and Data Analysis, v. 52, no. 3, pp. 1661--1673.
}
\seealso{\code{\link{is.aep4}}, \code{\link{paraep4}} }
\examples{
para <- vec2par(c(0,1, 0.5, 4), type="aep4")
if(are.paraep4.valid(para)) Q <- quaaep4(0.5,para)
}
\keyword{utility (distribution)}
\keyword{Distribution: Asymmetric Exponential Power}
\keyword{utility (distribution/parameter validation)}
\keyword{Distribution: Exponential Power}

