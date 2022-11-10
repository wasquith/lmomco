\name{are.parpdq4.valid}
\alias{are.parpdq4.valid}
\title{Are the Distribution Parameters Consistent with the Polynomial Density-Quantile4}
\description{
Is the distribution parameter object consistent with the corresponding distribution?  The distribution functions  (\code{\link{cdfpdq4}}, \code{\link{pdfpdq4}}, \code{\link{quapdq4}}, and
\code{\link{lmompdq4}}) require consistent parameters to return the cumulative
probability (nonexceedance), density, quantile, and L-moments of the distribution,
respectively.  These functions internally use the \code{\link{are.parpdq4.valid}}
function.
}
\usage{
are.parpdq4.valid(para, nowarn=FALSE)
}
\arguments{
  \item{para}{A distribution parameter list returned by \code{\link{parpdq4}} or \code{\link{vec2par}}.}
  \item{nowarn}{A logical switch on warning suppression. If \code{TRUE} then \code{options(warn=-1)} is made and restored on return. This switch is to permit calls in which warnings are not desired as the user knows how to handle the returned value---say in an optimization algorithm.}
}
\value{
  \item{TRUE}{If the parameters are \code{pdq4} consistent.}
  \item{FALSE}{If the parameters are not \code{pdq4} consistent.}
}
\note{
This function calls \code{\link{is.pdq4}} to verify consistency between
the distribution parameter object and the intent of the user.
}
\author{W.H. Asquith}
\seealso{\code{\link{is.pdq4}}, \code{\link{parpdq4}} }
\examples{
para <- parpdq4(lmoms(c(46, 70, 59, 36, 71, 48, 46, 63, 35, 52)))
if(are.parpdq4.valid(para)) Q <- quapdq4(0.5, para)
}
\keyword{utility (distribution)}
\keyword{Distribution: Polynomial Density-Quantile4}
\keyword{utility (distribution/parameter validation)}

