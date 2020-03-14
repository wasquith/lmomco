\name{are.parlmrq.valid}
\alias{are.parlmrq.valid}
\title{Are the Distribution Parameters Consistent with the Linear Mean Residual Quantile Function Distribution}
\description{
Is the distribution parameter object consistent with the corresponding distribution?  The distribution functions (\code{\link{cdflmrq}}, \code{\link{pdflmrq}}, \code{\link{qualmrq}}, and \code{\link{lmomlmrq}}) require consistent parameters to return the cumulative probability (nonexceedance), density, quantile, and L-moments of the distribution, respectively.  These functions internally use the \code{\link{are.parlmrq.valid}} function. The constraints on the parameters are listed under \code{\link{qualmrq}}. The documentation for \code{\link{qualmrq}} provides the conditions for valid parameters.
}
\usage{
are.parlmrq.valid(para, nowarn=FALSE)
}
\arguments{
  \item{para}{A distribution parameter list returned by \code{\link{parlmrq}} or \code{\link{vec2par}}.}
  \item{nowarn}{A logical switch on warning suppression. If \code{TRUE} then \code{options(warn=-1)} is made and restored on return. This switch is to permit calls in which warnings are not desired as the user knows how to handle the returned value---say in an optimization algorithm.}
}
\value{
  \item{TRUE}{If the parameters are \code{lmrq} consistent.}
  \item{FALSE}{If the parameters are not \code{lmrq} consistent.}
}
\references{
Midhu, N.N., Sankaran, P.G., and Nair, N.U., 2013, A class of distributions with linear mean residual quantile function and it's generalizations: Statistical Methodology, v. 15, pp. 1--24.
}
\note{
This function calls \code{\link{is.lmrq}} to verify consistency between the distribution parameter object and the intent of the user.
}
\author{ W.H. Asquith}
\seealso{ \code{\link{is.lmrq}}, \code{\link{parlmrq}} }
\examples{
para <- parlmrq(lmoms(c(3, 0.05, 1.6, 1.37, 0.57, 0.36, 2.2)))
if(are.parlmrq.valid(para)) Q <- qualmrq(0.5,para)
}
\keyword{utility (distribution)}
\keyword{Distribution: Linear Mean Residual Quantile Function}
\keyword{utility (distribution/parameter validation)}

